from __future__ import annotations

import math
from typing import Mapping, Sequence
import numpy as np

from .models import Material, Laminate, LaminateSection, ReducedBending


def qbar_orthotropic(mat: Material, theta_rad: float) -> np.ndarray:
    E1, E2, G, v1 = mat.E1, mat.E2, mat.G12, mat.nu12
    v2 = E2 / E1 * v1
    den = 1.0 - v1 * v2
    Q11 = E1 / den
    Q22 = E2 / den
    Q12 = v1 * E2 / den
    Q66 = G
    c, s = math.cos(theta_rad), math.sin(theta_rad)
    Q11t = c**4 * Q11 + s**4 * Q22 + 2*c*c*s*s*(Q12 + 2*Q66)
    Q22t = s**4 * Q11 + c**4 * Q22 + 2*c*c*s*s*(Q12 + 2*Q66)
    Q12t = c*c*s*s*(Q11 + Q22 - 4*Q66) + (c**4 + s**4)*Q12
    Q66t = c*c*s*s*(Q11 + Q22 - 2*Q12) + (c*c - s*s)**2 * Q66
    Q16t = c*s*(c*c*Q11 - s*s*Q22 - (c*c-s*s)*(Q12+2*Q66))
    Q26t = c*s*(s*s*Q11 - c*c*Q22 + (c*c-s*s)*(Q12+2*Q66))
    return np.array([[Q11t,Q12t,Q16t],[Q12t,Q22t,Q26t],[Q16t,Q26t,Q66t]], float)


def material_invariants(mat: Material) -> np.ndarray:
    nu21 = mat.nu12 * mat.E2 / mat.E1
    den = 1.0 - mat.nu12 * nu21
    Q11 = mat.E1 / den
    Q22 = mat.E2 / den
    Q12 = mat.nu12 * mat.E2 / den
    Q66 = mat.G12
    return np.array([
        (3*Q11+3*Q22+2*Q12+4*Q66)/8,
        (Q11-Q22)/2,
        (Q11+Q22-2*Q12-4*Q66)/8,
        (Q11+Q22+6*Q12-4*Q66)/8,
        (Q11+Q22-2*Q12+4*Q66)/8,
    ], float)


_GL3_X = np.array([-math.sqrt(0.6), 0.0, math.sqrt(0.6)])
_GL3_W = np.array([5.0, 8.0, 5.0]) / 9.0


def laminate_section(lam: Laminate, materials: Mapping[str, Material]) -> LaminateSection:
    if not lam.plies:
        raise ValueError("Laminate must contain at least one ply")
    tt = lam.thickness
    z = np.concatenate(([-tt/2], -tt/2 + np.cumsum([p.thickness for p in lam.plies])))
    n = len(lam.plies)
    qbars = np.zeros((n,3,3), float)
    gth = np.zeros(n, float)
    for k,p in enumerate(lam.plies):
        if p.material not in materials:
            raise KeyError(f"Unknown material {p.material!r} in ply {k+1}")
        m = materials[p.material]
        qbars[k] = qbar_orthotropic(m, math.radians(p.angle_deg))
        gth[k] = float(m.Gth if m.Gth is not None else m.G12)
    A=np.zeros((3,3)); B=np.zeros((3,3)); D=np.zeros((3,3)); A44=0.0
    for k in range(n):
        z0,z1=z[k],z[k+1]; qb=qbars[k]
        A += qb*(z1-z0)
        B += 0.5*qb*(z1*z1-z0*z0)
        D += (1/3)*qb*(z1**3-z0**3)
        A44 += gth[k]*(z1-z0)
    KABD=np.block([[A,B],[B,D]])
    CABD=np.linalg.solve(KABD,np.eye(6))
    d=-CABD[0,3]/CABD[3,3]
    Dn=D-2*d*B+d*d*A
    Bn=B-d*A
    # Integral of g(z)^2 / G_k through each ply, g(z) the first moment of Q11 about the neutral axis. g is
    # quadratic in z, so the integrand is a quartic: three-point Gauss-Legendre integration is exact
    # (identical to the adaptive quadrature of earlier versions, about 30 times faster).
    Gssum=0.0; gcum=0.0
    for k in range(n):
        z0,z1=z[k],z[k+1]; q11=qbars[k,0,0]; gk=gth[k]
        half=0.5*(z1-z0); mid=0.5*(z1+z0)
        zz=mid+half*_GL3_X
        g=gcum+0.5*q11*((zz-d)**2-(z0-d)**2)
        Gssum += half*float(_GL3_W@(g*g))/gk
        gcum += 0.5*q11*((z1-d)**2-(z0-d)**2)
    if Gssum <= 0:
        raise ValueError("Unable to compute positive equivalent transverse shear stiffness")
    Seq=Dn[0,0]**2/Gssum
    kappa=Seq/A44
    scaleB=max(1.0,np.linalg.norm(A,'fro')*tt)
    sym=float(np.linalg.norm(B,'fro')/scaleB)
    return LaminateSection(A=A,B=B,Bneutral=Bn,Doriginal=D,D=Dn,Seq=float(Seq),A44=float(A44),
                           kappa_shear=float(kappa),d_neutral=float(d),thickness=float(tt),z=z,Qbar=qbars,
                           symmetry_indicator=sym)


# Calibrated twist release (theory guide, Section 5.2): eta = (A + B psi^2) lambda^N with the twist length
# lambda = sqrt(D66 / A44) / L of the wall; fitted to periodic shell FE of 320 random symmetric walls
# (0.010 <= lambda <= 0.100). eta is the fraction of the bend-twist reduction psi^2 of D11 that the wall releases.
TWIST_RELEASE_A = 2.0
TWIST_RELEASE_B = 2.15
TWIST_RELEASE_N = 0.73
TWIST_LENGTH_RANGE = (0.010, 0.100)


def twist_length(sec, length: float) -> float:
    """lambda = sqrt(D66 / A44) / L: length over which the wall's transverse shear resists twist, relative to L."""
    D = np.asarray(sec.D, float)
    if not (length and length > 0 and sec.A44 > 0 and D[2, 2] > 0):
        return 0.0
    return float(math.sqrt(D[2, 2] / sec.A44) / length)


TWIST_RELEASE_DEFAULT = (TWIST_RELEASE_A, TWIST_RELEASE_B, TWIST_RELEASE_N)


def twist_release(psi: float, lam: float, coeffs=None) -> float:
    """eta = (a + b psi^2) lambda^n, capped at 1; coeffs = (a, b, n), default the calibrated values."""
    a, b, n = TWIST_RELEASE_DEFAULT if coeffs is None else (float(c) for c in coeffs)
    if lam <= 0:
        return 0.0
    return min(1.0, max(0.0, (a + b * psi * psi) * lam ** n))


def reduced_bending_stiffness(sec: LaminateSection, mode: str="freeMoment", length: float | None = None, twist_coeffs=None) -> ReducedBending:
    A=sec.A; D0=sec.Doriginal; B0=sec.B; t=sec.thickness
    bscale=max(np.linalg.norm(A,'fro')*max(t,np.finfo(float).eps),np.finfo(float).eps)
    sym=float(np.linalg.norm(B0,'fro')/bscale)
    if sym>1e-10:
        Dsec=D0-B0@np.linalg.solve(A,B0)
        used=True
    else:
        Dsec=D0.copy(); used=False
    Dsec=0.5*(Dsec+Dsec.T)
    key=mode.lower().replace('_','')
    if key in ("freemoment","free","anisotropic"):
        Drr=Dsec[1:3,1:3]; Dr1=Dsec[1:3,0]
        if np.linalg.cond(Drr)>1e13:
            raise np.linalg.LinAlgError("Transverse anisotropic bending submatrix is singular/ill-conditioned")
        ratio=-np.linalg.solve(Drr,Dr1)
        dstrip=float(Dsec[0,0]-Dsec[0,1:3]@np.linalg.solve(Drr,Dr1))
    elif key in ("cylindrical","thesis"):
        ratio=np.zeros(2); dstrip=float(Dsec[0,0])
    elif key in ("twistrelease","twistcorrected","twist"):
        # Calibrated partial twist release (theory guide, Section 5.2):
        # D_b = D11 (1 - eta psi^2),  eta = (A + B psi^2) lambda^N, capped at 1
        psi=float(Dsec[0,2]/math.sqrt(Dsec[0,0]*Dsec[2,2])) if Dsec[0,0]>0 and Dsec[2,2]>0 else 0.0
        eta=twist_release(psi,twist_length(sec,length),twist_coeffs)
        ratio=np.array([0.0,0.0]); dstrip=float(Dsec[0,0]*(1.0-eta*psi*psi))
    else:
        raise ValueError(f"Unknown strip bending mode: {mode}")
    if not math.isfinite(dstrip) or dstrip<=0:
        raise ValueError(f"Reduced strip bending stiffness must be positive; got {dstrip}")
    return ReducedBending(Dfull=D0,Dsection=Dsec,Dstrip=dstrip,curvature_ratio=ratio,
                          used_B_reduction=used,symmetry_measure=sym,mode=mode)


def lamination_parameters(lam: Laminate, materials: Mapping[str, Material]) -> dict:
    if not lam.plies:
        raise ValueError("Laminate must contain plies")
    T=lam.thickness
    z=np.concatenate(([-T/2],-T/2+np.cumsum([p.thickness for p in lam.plies])))
    xiA=np.zeros(4); xiD=np.zeros(4)
    used=[]
    for p in lam.plies:
        if p.material not in used: used.append(p.material)
    by={m:{"material":m,"rhoA":0.0,"xiA":np.zeros(4),"rhoD":0.0,"xiD":np.zeros(4)} for m in used}
    for k,p in enumerate(lam.plies):
        th=math.radians(p.angle_deg)
        f=np.array([math.cos(2*th),math.cos(4*th),math.sin(2*th),math.sin(4*th)])
        z0,z1=z[k],z[k+1]; dz=z1-z0; d3=z1**3-z0**3
        xiA+=(dz/T)*f; xiD+=(4*d3/T**3)*f
        b=by[p.material]; b["rhoA"]+=dz/T; b["xiA"]+=(dz/T)*f; b["rhoD"]+=4*d3/T**3; b["xiD"]+=(4*d3/T**3)*f
    return {"thickness":T,"z":z,"xiA":xiA,"xiD":xiD,"by_material":list(by.values()),
            "sum_rhoA":sum(v["rhoA"] for v in by.values()),"sum_rhoD":sum(v["rhoD"] for v in by.values())}


def _make_matrix(U: np.ndarray, xi: Sequence[float], scale: float, rho: float=1.0) -> np.ndarray:
    x1,x2,x3,x4=map(float,xi)
    M11=scale*(U[0]*rho+U[1]*x1+U[2]*x2)
    M22=scale*(U[0]*rho-U[1]*x1+U[2]*x2)
    M12=scale*(U[3]*rho-U[2]*x2)
    M66=scale*(U[4]*rho-U[2]*x2)
    M16=scale*(0.5*U[1]*x3+U[2]*x4)
    M26=scale*(0.5*U[1]*x3-U[2]*x4)
    return np.array([[M11,M12,M16],[M12,M22,M26],[M16,M26,M66]],float)


def section_from_invariants_lp(mat: Material,t:float,xiA:Sequence[float],xiD:Sequence[float],Seq:float) -> LaminateSection:
    U=material_invariants(mat); A=_make_matrix(U,xiA,t); D=_make_matrix(U,xiD,t**3/12)
    z=np.array([-t/2,t/2]); q=np.zeros((0,3,3))
    return LaminateSection(A=A,B=np.zeros((3,3)),Bneutral=np.zeros((3,3)),Doriginal=D,D=D,Seq=float(Seq),A44=math.nan,
                           kappa_shear=math.nan,d_neutral=0.0,thickness=t,z=z,Qbar=q,symmetry_indicator=0.0)


def section_from_hybrid_lp(materials: Mapping[str,Material], blocks: Sequence[dict], t:float, Seq:float) -> LaminateSection:
    A=np.zeros((3,3)); D=np.zeros((3,3))
    for block in blocks:
        m=materials[block["material"]]; U=material_invariants(m)
        A+=_make_matrix(U,block["xiA"],t,float(block["rhoA"]))
        D+=_make_matrix(U,block["xiD"],t**3/12,float(block["rhoD"]))
    return LaminateSection(A=A,B=np.zeros((3,3)),Bneutral=np.zeros((3,3)),Doriginal=D,D=D,Seq=float(Seq),A44=math.nan,
                           kappa_shear=math.nan,d_neutral=0.0,thickness=t,z=np.array([-t/2,t/2]),Qbar=np.zeros((0,3,3)),symmetry_indicator=0.0)


def analyze_laminate(lam: Laminate, materials: Mapping[str,Material]) -> dict:
    sec=lam.section_override or laminate_section(lam,materials)
    lp=lamination_parameters(lam,materials)
    inv=[]
    seen=[]
    for p in lam.plies:
        if p.material in seen: continue
        seen.append(p.material); m=materials[p.material]
        inv.append({"material":m.name,"E1":m.E1,"E2":m.E2,"G12":m.G12,"nu12":m.nu12,"U":material_invariants(m)})
    return {"section":sec,"A":sec.A,"B":sec.B,"D":sec.Doriginal,"Dneutral":sec.D,
            "invariants":inv,"lp":lp,"strip_cyl":reduced_bending_stiffness(sec,"cylindrical"),
            "strip_free":reduced_bending_stiffness(sec,"freeMoment")}
