from __future__ import annotations
import numpy as np
from .models import LaminateSection
from .laminate import reduced_bending_stiffness


WALL_SHEAR_OPTIONS = ("euler-bernoulli", "timoshenko")


def is_timoshenko(wall_shear: str) -> bool:
    """True for the shear-deformable wall; False for the Euler-Bernoulli wall of the original model."""
    w = str(wall_shear or "euler-bernoulli").lower().replace("_", "-").replace(" ", "-")
    if w.startswith("timo"):
        return True
    if w in ("euler-bernoulli", "eulerbernoulli", "eb", "bernoulli", "euler"):
        return False
    raise ValueError(f"Unknown wall_shear: {wall_shear!r} (use 'euler-bernoulli' or 'timoshenko')")


def _wall_matrix(A: np.ndarray, Db: float, Seq: float, L: float, h: float, shear: bool):
    """Wall stiffness: membrane from A, out-of-plane bending from Db.

    Timoshenko (shear=True): beta = 12 Db / (L^2 S_eq); Euler-Bernoulli (shear=False): beta = 0, i.e.
    the classical 12 Db h / L^3 beam terms of the original model.
    """
    beta = 12 * Db / (L * L * Seq) if shear else 0.0
    s = Seq * L * L * beta * beta if shear else 0.0          # = 144 Db^2 / (L^2 S_eq)
    X = s + 12 * Db
    kb = h * X / (L ** 3 * (beta + 1) ** 2)
    kbr = h * X / (2 * L ** 2 * (beta + 1) ** 2)
    krD = h * (s + 4 * Db * beta * beta + 8 * Db * beta + 16 * Db) / (4 * L * (beta + 1) ** 2)
    krO = h * (s - 4 * Db * beta * beta - 8 * Db * beta + 8 * Db) / (4 * L * (beta + 1) ** 2)
    K = np.zeros((8, 8))
    M = h / L * np.array([[A[0, 0], A[0, 2]], [A[0, 2], A[2, 2]]])
    ii = [0, 2]; jj = [4, 6]
    K[np.ix_(ii, ii)] += M; K[np.ix_(ii, jj)] -= M; K[np.ix_(jj, ii)] -= M; K[np.ix_(jj, jj)] += M
    ib = [1, 3, 5, 7]
    K[np.ix_(ib, ib)] += np.array([[kb, kbr, -kb, kbr], [kbr, krD, -kbr, krO], [-kb, -kbr, kb, -kbr], [kbr, krO, -kbr, krD]])
    g = h * np.array([-A[0, 1], 0, -A[1, 2], 0, A[0, 1], 0, A[1, 2], 0], float)
    r = float(A[1, 1] * L * h)
    return 0.5 * (K + K.T), g, r, beta


def wall_element_paper(A:np.ndarray,D:np.ndarray,Seq:float,L:float,h:float,shear:bool=True):
    """Cylindrical-bending wall (D11). ``shear=False`` is the Euler-Bernoulli wall of the original model."""
    if Seq<=0 or D[0,0]<=0 or L<=0 or h<=0:
        raise ValueError("D11, Seq, L and h must be positive")
    K,g,r,beta=_wall_matrix(A,float(D[0,0]),Seq,L,h,shear)
    return K,g,r,{"mode":"paper","Dstrip":float(D[0,0]),"curvature_ratio":np.zeros(2),"beta":beta,
                  "beam":"timoshenko" if shear else "euler-bernoulli"}


def wall_element_anisotropic_strip(sec:LaminateSection,L:float,h:float,mode:str="freeMoment",shear:bool=True,twist_coeffs=None):
    if L<=0 or h<=0 or sec.Seq<=0: raise ValueError("L, h and Seq must be positive")
    red=reduced_bending_stiffness(sec,mode,L,twist_coeffs)
    K,g,r,beta=_wall_matrix(sec.A,red.Dstrip,sec.Seq,L,h,shear)
    meta={"mode":red.mode,"Dstrip":red.Dstrip,"curvature_ratio":red.curvature_ratio,"used_B_reduction":red.used_B_reduction,
          "symmetry_measure":red.symmetry_measure,"beta":beta,"beam":"timoshenko" if shear else "euler-bernoulli",
          "twist_coeffs":None if twist_coeffs is None else [float(c) for c in twist_coeffs]}
    return K,g,r,meta


def wall_transform(alpha:float)->np.ndarray:
    c,s=np.cos(alpha),np.sin(alpha); R=np.array([[c,-s],[s,c]])
    T=np.eye(8); T[0:2,0:2]=R; T[4:6,4:6]=R
    return T

# Backward-compatible alias for projects built with earlier releases.
wall_element_thesis = wall_element_paper
