from __future__ import annotations

"""Post-processing helpers for batch homogenization results.

The module deliberately works only with quantities already produced by RVE Studio.
It does not introduce strength estimates or assume a parent-solid modulus.  The
Gibson-Ashby style plots are therefore scaling plots of effective stiffness versus
relative density rather than fitted constitutive laws for a particular material.
"""

from dataclasses import dataclass
import re
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PowerLawFit:
    coefficient: float
    exponent: float
    r_squared: float
    n_points: int


def _num(s: pd.Series | None, index) -> pd.Series:
    if s is None:
        return pd.Series(np.nan, index=index, dtype=float)
    return pd.to_numeric(s, errors="coerce")


def _safe_div(a: pd.Series, b: pd.Series) -> pd.Series:
    aa = pd.to_numeric(a, errors="coerce").astype(float)
    bb = pd.to_numeric(b, errors="coerce").astype(float)
    out = pd.Series(np.nan, index=aa.index, dtype=float)
    mask = np.isfinite(aa) & np.isfinite(bb) & (np.abs(bb) > 0.0)
    out.loc[mask] = aa.loc[mask] / bb.loc[mask]
    return out


def _ratio_ge_one(a: pd.Series, b: pd.Series) -> pd.Series:
    aa = np.abs(pd.to_numeric(a, errors="coerce").astype(float))
    bb = np.abs(pd.to_numeric(b, errors="coerce").astype(float))
    mx = pd.concat([aa, bb], axis=1).max(axis=1)
    mn = pd.concat([aa, bb], axis=1).min(axis=1)
    return _safe_div(mx, mn)


def enrich_results(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with derived efficiency and anisotropy metrics.

    Specific moduli use the numerically convenient unit GPa/(Mg/m^3).  For an
    effective modulus supplied in MPa and density in kg/m^3, the numerical value
    is simply ``E_MPa / rho_kg_m3`` because the 10^3 conversions cancel.
    """
    out = df.copy()
    idx = out.index
    rho = _num(out.get("AverageDensityKgM3"), idx)
    rrel = _num(out.get("RelativeDensity"), idx)
    out["RelativeDensityPct"] = 100.0 * rrel
    out["AverageDensityMgM3"] = rho / 1000.0

    stiffness_names = ["E1", "E2", "E3", "G12", "G13", "G23", "Q11", "Q22", "Q12", "Q66", "Q44", "Q55", "C11", "C22", "C33", "C44", "C55", "C66"]
    for name in stiffness_names:
        if name in out.columns:
            out[f"Specific{name}_GPa_per_Mg_m3"] = _safe_div(_num(out[name], idx), rho)

    if "E1" in out and "E2" in out:
        out["E1_over_E2"] = _safe_div(_num(out["E1"], idx), _num(out["E2"], idx))
        out["E12_Anisotropy"] = _ratio_ge_one(out["E1"], out["E2"])
    if "Q11" in out and "Q22" in out:
        out["Q11_over_Q22"] = _safe_div(_num(out["Q11"], idx), _num(out["Q22"], idx))
        out["Q12_Anisotropy"] = _ratio_ge_one(out["Q11"], out["Q22"])
    if "G13" in out and "G23" in out:
        out["G13_over_G23"] = _safe_div(_num(out["G13"], idx), _num(out["G23"], idx))
        out["OutOfPlaneShear_Anisotropy"] = _ratio_ge_one(out["G13"], out["G23"])
    if all(c in out for c in ("Q11", "Q22", "Q12")):
        den = np.sqrt(np.abs(_num(out["Q11"], idx) * _num(out["Q22"], idx)))
        out["Q12_CouplingIndex"] = _safe_div(np.abs(_num(out["Q12"], idx)), den)

    # Geometry anisotropy: preserve directional ratio and also supply a >=1 aspect ratio.
    if "LengthA_mm" in out and "LengthB_mm" in out:
        la, lb = _num(out["LengthA_mm"], idx), _num(out["LengthB_mm"], idx)
        out["GeometryLengthRatio"] = _safe_div(la, lb)
        out["GeometryAspectRatio"] = _ratio_ge_one(la, lb)
    elif "LengthRatio" in out:
        lr = np.abs(_num(out["LengthRatio"], idx))
        out["GeometryLengthRatio"] = _num(out["LengthRatio"], idx)
        inv = _safe_div(pd.Series(1.0, index=idx), lr)
        out["GeometryAspectRatio"] = pd.concat([lr, inv], axis=1).max(axis=1)

    # A compact in-plane efficiency score useful for ranking; it is not a strength metric.
    if all(c in out for c in ("E1", "E2")):
        emin = pd.concat([_num(out["E1"], idx), _num(out["E2"], idx)], axis=1).min(axis=1)
        out["MinInPlaneSpecificE_GPa_per_Mg_m3"] = _safe_div(emin, rho)

    # Cost-efficiency metrics use the user-defined material cost units.  They are
    # intentionally unit-agnostic: if the library stores GBP/kg, the results are
    # in GBP; if it stores a normalized dummy index/kg, the same index is used.
    cost = _num(out.get("TotalMaterialCost"), idx)
    if "TotalMaterialCost" in out.columns:
        for name in ["E1","E2","E3","G12","G13","G23","Q11","Q22","Q66","C11","C22","C33","C44","C55","C66"]:
            if name in out.columns:
                out[f"{name}_per_Cost"] = _safe_div(_num(out[name], idx), cost)
        if "MassKg" in out.columns:
            out["CostPerKgFromCase"] = _safe_div(cost, _num(out["MassKg"], idx))

    # Manufacturability descriptors.  These are screening quantities only; they
    # do not modify the homogenisation solution.  New v1.3.4 batch files carry
    # the T1_mm, T2_mm, ... design columns so these can be evaluated directly.
    tcols=[c for c in out.columns if re.fullmatch(r"T\d+_mm",str(c))]
    if tcols:
        T=out[tcols].apply(pd.to_numeric,errors="coerce")
        out["MinWallThickness_mm"]=T.min(axis=1,skipna=True)
        out["MaxWallThickness_mm"]=T.max(axis=1,skipna=True)
        out["MeanWallThickness_mm"]=T.mean(axis=1,skipna=True)
        out["WallThicknessRatio"]=_safe_div(out["MaxWallThickness_mm"],out["MinWallThickness_mm"])
        out["WallThicknessSpreadPct"]=100.0*_safe_div(out["MaxWallThickness_mm"]-out["MinWallThickness_mm"],out["MeanWallThickness_mm"])

    if "LayupAngles" in out.columns:
        def parse_angles(value):
            if pd.isna(value): return []
            parts=re.split(r"[/,;\s]+",str(value).strip())
            vals=[]
            for x in parts:
                if not x: continue
                try: vals.append(float(x))
                except Exception: pass
            return vals
        def balanced(a):
            if not a: return False
            vals=[round(((x+180.0)%180.0),8) for x in a]
            # 0 and 90 degree plies are self-balanced; every other +theta must
            # have the corresponding -theta modulo 180.
            for x in a:
                xm=((x+180.0)%180.0)
                if min(abs(xm),abs(xm-90.0),abs(xm-180.0))<1e-7: continue
                target=round(((-x)+180.0)%180.0,8)
                if vals.count(round(xm,8))!=vals.count(target): return False
            return True
        parsed=out["LayupAngles"].map(parse_angles)
        out["PlyCount"]=parsed.map(len).astype(float)
        out["SymmetricLayup"]=parsed.map(lambda a: bool(a) and all(abs(x-y)<1e-8 for x,y in zip(a,reversed(a))))
        out["BalancedLayup"]=parsed.map(balanced)
    return out


def fit_power_law(x, y) -> PowerLawFit | None:
    """Fit y = C x^n in log space for positive finite observations."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if int(mask.sum()) < 3:
        return None
    lx = np.log(x[mask]); ly = np.log(y[mask])
    n, logc = np.polyfit(lx, ly, 1)
    pred = logc + n * lx
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return PowerLawFit(float(np.exp(logc)), float(n), float(r2), int(mask.sum()))


def pareto_mask_min_x_max_y(x, y) -> np.ndarray:
    """Boolean Pareto mask for minimizing x and maximizing y."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    result = np.zeros(len(x), dtype=bool)
    order = np.where(valid)[0][np.argsort(x[valid], kind="mergesort")]
    best = -np.inf
    for i in order:
        if y[i] > best:
            result[i] = True
            best = y[i]
    return result


def directional_youngs_modulus(E1: float, E2: float, G12: float, nu12: float, angles_deg=None):
    """Orthotropic plane-stress E(theta) derived from E1,E2,G12,nu12.

    Returns angles in degrees and E(theta) in the same units as E1/E2/G12.
    """
    if angles_deg is None:
        angles_deg = np.linspace(0.0, 360.0, 361)
    theta = np.deg2rad(np.asarray(angles_deg, dtype=float))
    E1=float(E1); E2=float(E2); G12=float(G12); nu12=float(nu12)
    if min(E1, E2, G12) <= 0:
        return np.asarray(angles_deg, dtype=float), np.full_like(theta, np.nan)
    s11=1.0/E1; s22=1.0/E2; s12=-nu12/E1; s66=1.0/G12
    c=np.cos(theta); s=np.sin(theta)
    invE=s11*c**4+s22*s**4+(2.0*s12+s66)*(s*c)**2
    E=np.where(invE>0,1.0/invE,np.nan)
    return np.asarray(angles_deg,dtype=float), E


def merge_input_metadata(results: pd.DataFrame, input_df: pd.DataFrame) -> pd.DataFrame:
    """Attach design metadata to an older result file using CaseID."""
    if "CaseID" not in results or "CaseID" not in input_df:
        raise ValueError("Both files must contain CaseID")
    meta_cols = [c for c in input_df.columns if c != "CaseID" and c not in results.columns]
    return results.merge(input_df[["CaseID", *meta_cols]], on="CaseID", how="left", validate="one_to_one")


@dataclass(frozen=True)
class ObjectiveSpec:
    """One database-based design objective.

    goal is ``max``, ``min`` or ``target``.  ``target`` minimizes the absolute
    relative deviation from ``target`` (absolute deviation when target=0).
    """
    metric: str
    goal: str = "max"
    target: float | None = None
    weight: float = 1.0


@dataclass(frozen=True)
class ConstraintSpec:
    metric: str
    relation: str
    value: float
    upper: float | None = None


def _objective_loss(values: np.ndarray, spec: ObjectiveSpec) -> np.ndarray:
    x=np.asarray(values,dtype=float); goal=str(spec.goal).strip().lower()
    if goal.startswith('max'): return -x
    if goal.startswith('min'): return x
    if goal.startswith('target'):
        if spec.target is None or not np.isfinite(float(spec.target)):
            raise ValueError(f"Target objective {spec.metric!r} requires a finite target")
        t=float(spec.target); scale=abs(t) if abs(t)>1e-15 else 1.0
        return np.abs(x-t)/scale
    raise ValueError(f"Unknown objective goal {spec.goal!r}; use max, min or target")


def apply_constraints(df: pd.DataFrame, constraints: list[ConstraintSpec] | None) -> np.ndarray:
    """Return a Boolean feasible mask for numeric design constraints."""
    mask=np.ones(len(df),dtype=bool)
    for c in constraints or []:
        if c.metric not in df.columns:
            raise KeyError(f"Constraint metric {c.metric!r} is not in the result table")
        x=pd.to_numeric(df[c.metric],errors='coerce').to_numpy(float); rel=str(c.relation).strip().lower()
        m=np.isfinite(x)
        if rel in ('<=','≤','max'): m &= x <= float(c.value)
        elif rel in ('>=','≥','min'): m &= x >= float(c.value)
        elif rel in ('=','==','target'):
            tol=abs(float(c.upper)) if c.upper is not None else 0.0; m &= np.abs(x-float(c.value)) <= tol
        elif rel in ('between','range'):
            if c.upper is None: raise ValueError(f"Between constraint {c.metric!r} requires an upper value")
            lo,hi=sorted((float(c.value),float(c.upper))); m &= (x>=lo)&(x<=hi)
        else: raise ValueError(f"Unknown constraint relation {c.relation!r}")
        mask &= m
    return mask


def pareto_mask_multiobjective(df: pd.DataFrame, objectives: list[ObjectiveSpec], feasible_mask=None) -> np.ndarray:
    """Return non-dominated cases for arbitrary min/max/target objectives.

    The search is over the finite rows of the supplied batch database; this is
    database-based Pareto optimization, not a surrogate or continuous optimizer.
    """
    if not objectives: raise ValueError('At least one objective is required')
    n=len(df); feasible=np.ones(n,dtype=bool) if feasible_mask is None else np.asarray(feasible_mask,dtype=bool).copy()
    losses=[]
    for o in objectives:
        if o.metric not in df.columns: raise KeyError(f"Objective metric {o.metric!r} is not in the result table")
        x=pd.to_numeric(df[o.metric],errors='coerce').to_numpy(float); loss=_objective_loss(x,o); feasible &= np.isfinite(loss); losses.append(loss)
    out=np.zeros(n,dtype=bool); idx=np.where(feasible)[0]
    if not len(idx): return out
    L=np.column_stack(losses)[idx]
    keep=np.ones(len(idx),dtype=bool)
    # A point is dominated when another feasible design is no worse in every
    # objective and strictly better in at least one.
    for i in range(len(idx)):
        if not keep[i]: continue
        dominated_by=np.all(L <= L[i],axis=1) & np.any(L < L[i],axis=1)
        dominated_by[i]=False
        if np.any(dominated_by): keep[i]=False
    out[idx[keep]]=True
    return out


def weighted_tradeoff_score(df: pd.DataFrame, objectives: list[ObjectiveSpec], feasible_mask=None) -> np.ndarray:
    """Return a 0..1 normalized weighted utility for feasible database designs."""
    if not objectives: raise ValueError('At least one objective is required')
    n=len(df); feasible=np.ones(n,dtype=bool) if feasible_mask is None else np.asarray(feasible_mask,dtype=bool).copy()
    raw=[]; weights=[]
    for o in objectives:
        if o.metric not in df.columns: raise KeyError(f"Objective metric {o.metric!r} is not in the result table")
        x=pd.to_numeric(df[o.metric],errors='coerce').to_numpy(float); loss=_objective_loss(x,o); feasible &= np.isfinite(loss); raw.append(loss); weights.append(max(0.0,float(o.weight)))
    score=np.full(n,np.nan,dtype=float); idx=np.where(feasible)[0]
    if not len(idx): return score
    W=np.asarray(weights,dtype=float)
    if not np.any(W>0): W[:]=1.0
    W=W/W.sum(); utilities=[]
    for loss in raw:
        v=loss[idx]; lo=float(np.min(v)); hi=float(np.max(v))
        if hi-lo <= 1e-15: u=np.ones_like(v)
        else: u=1.0-(v-lo)/(hi-lo)
        utilities.append(np.clip(u,0.0,1.0))
    score[idx]=np.column_stack(utilities)@W
    return score


def optimize_batch_designs(df: pd.DataFrame, objectives: list[ObjectiveSpec], constraints: list[ConstraintSpec] | None=None, manufacturing_rules=None) -> pd.DataFrame:
    """Add feasibility, Pareto and weighted trade-off columns to batch results.

    ``manufacturing_rules`` is optional and retained as a late-bound argument so
    older scripts calling the v1.3.3 API remain compatible.
    """
    out=enrich_results(df.copy()); feasible=apply_constraints(out,constraints)
    if manufacturing_rules is not None:
        # Function is defined later in the module; resolution occurs at call time.
        mfg=manufacturing_feasible_mask(out,manufacturing_rules); out['ManufacturingFeasible']=mfg; feasible &= mfg
    pareto=pareto_mask_multiobjective(out,objectives,feasible); score=weighted_tradeoff_score(out,objectives,feasible)
    out['OptimizationFeasible']=feasible
    out['OptimizationPareto']=pareto
    out['OptimizationScore']=score
    # Store target deviations explicitly so target-matching studies are auditable.
    for o in objectives:
        if str(o.goal).lower().startswith('target'):
            x=pd.to_numeric(out[o.metric],errors='coerce').to_numpy(float); t=float(o.target); scale=abs(t) if abs(t)>1e-15 else 1.0
            out[f'TargetDeviation_{o.metric}']=np.abs(x-t)/scale
    return out


@dataclass(frozen=True)
class ManufacturingRules:
    """Optional manufacturability screening rules for solved database cases.

    A ``None`` numerical value disables that rule.  These checks operate on the
    design metadata saved with each case; they do not alter the mechanics model.
    """
    min_wall_thickness_mm: float | None = None
    max_wall_thickness_mm: float | None = None
    max_wall_thickness_ratio: float | None = None
    max_geometry_aspect_ratio: float | None = None
    max_relative_density: float | None = None
    min_plies: int | None = None
    max_plies: int | None = None
    require_balanced_layup: bool = False
    require_symmetric_layup: bool = False
    allowed_ply_angles: tuple[float, ...] | None = None
    min_buckling_factor: float | None = None


def manufacturing_feasible_mask(df: pd.DataFrame, rules: ManufacturingRules | None) -> np.ndarray:
    """Return a manufacturability mask using derived design descriptors."""
    if rules is None:
        return np.ones(len(df),dtype=bool)
    d=enrich_results(df)
    mask=np.ones(len(d),dtype=bool)
    def apply_numeric(col,op,value):
        nonlocal mask
        if value is None: return
        if col not in d.columns:
            mask &= False; return
        x=pd.to_numeric(d[col],errors='coerce').to_numpy(float); m=np.isfinite(x)
        if op=='>=': m &= x>=float(value)
        else: m &= x<=float(value)
        mask &= m
    apply_numeric('MinWallThickness_mm','>=',rules.min_wall_thickness_mm)
    apply_numeric('MaxWallThickness_mm','<=',rules.max_wall_thickness_mm)
    apply_numeric('WallThicknessRatio','<=',rules.max_wall_thickness_ratio)
    apply_numeric('GeometryAspectRatio','<=',rules.max_geometry_aspect_ratio)
    apply_numeric('RelativeDensity','<=',rules.max_relative_density)
    apply_numeric('PlyCount','>=',rules.min_plies)
    apply_numeric('PlyCount','<=',rules.max_plies)
    apply_numeric('BucklingLambdaMin','>=',rules.min_buckling_factor)
    if rules.require_balanced_layup:
        if 'BalancedLayup' not in d.columns: mask &= False
        else: mask &= d['BalancedLayup'].fillna(False).astype(bool).to_numpy()
    if rules.require_symmetric_layup:
        if 'SymmetricLayup' not in d.columns: mask &= False
        else: mask &= d['SymmetricLayup'].fillna(False).astype(bool).to_numpy()
    if rules.allowed_ply_angles:
        if 'LayupAngles' not in d.columns: mask &= False
        else:
            allowed=np.asarray([((float(a)+180.0)%180.0) for a in rules.allowed_ply_angles],dtype=float)
            def ok_angles(value):
                if pd.isna(value): return False
                vals=[]
                for token in re.split(r'[/,;\s]+',str(value).strip()):
                    if token:
                        try: vals.append((float(token)+180.0)%180.0)
                        except Exception: return False
                return bool(vals) and all(np.min(np.abs(allowed-v))<1e-7 for v in vals)
            mask &= d['LayupAngles'].map(ok_angles).to_numpy(bool)
    return mask


def design_variable_columns(df: pd.DataFrame) -> tuple[list[str],list[str]]:
    """Infer numerical and categorical design variables from a batch table."""
    numeric=[]
    preferred={'Height','LengthA_mm','LengthB_mm','LengthRatio','TargetRelativeDensity','MaterialCostPerKg'}
    for c in df.columns:
        sc=str(c)
        if (sc in preferred or re.fullmatch(r'Param\d+Value',sc) or re.fullmatch(r'T\d+_mm',sc)) and pd.api.types.is_numeric_dtype(df[c]):
            numeric.append(sc)
    categorical=[c for c in ['GeometryID','GeometryType','Material','LayupID','LayupAngles'] if c in df.columns]
    return numeric,categorical


def sensitivity_analysis(df: pd.DataFrame, response: str, variables: list[str] | None=None, include_categories: bool=True) -> pd.DataFrame:
    """Screen sensitivity of a solved database response to design variables.

    Numerical variables report Pearson r, Spearman rho and standardized linear
    regression beta.  Categorical variables report eta-squared, the fraction of
    response variance explained by group means.  The analysis is descriptive and
    should not be interpreted as proof of causality.
    """
    if response not in df.columns: raise KeyError(f'Response {response!r} is not in the result table')
    yall=pd.to_numeric(df[response],errors='coerce')
    inferred_num,inferred_cat=design_variable_columns(df)
    numvars=inferred_num if variables is None else [v for v in variables if v in df.columns and pd.api.types.is_numeric_dtype(df[v])]
    catvars=inferred_cat if variables is None else [v for v in variables if v in df.columns and not pd.api.types.is_numeric_dtype(df[v])]
    rows=[]
    # Standardized regression is fitted jointly to finite, non-constant numerical variables.
    beta_map={}
    usable=[]
    for v in numvars:
        x=pd.to_numeric(df[v],errors='coerce')
        m=np.isfinite(x)&np.isfinite(yall)
        if int(m.sum())>=3 and float(x[m].std(ddof=0))>1e-15: usable.append(v)
    if usable:
        block=pd.DataFrame({v:pd.to_numeric(df[v],errors='coerce') for v in usable}); block['_y']=yall
        block=block.replace([np.inf,-np.inf],np.nan).dropna()
        if len(block)>=max(3,len(usable)+1):
            X=block[usable].to_numpy(float); y=block['_y'].to_numpy(float)
            sx=X.std(axis=0); sy=y.std()
            keep=sx>1e-15
            if sy>1e-15 and np.any(keep):
                Z=(X[:,keep]-X[:,keep].mean(axis=0))/sx[keep]; zy=(y-y.mean())/sy
                beta=np.linalg.lstsq(Z,zy,rcond=None)[0]
                for name,b in zip(np.asarray(usable)[keep],beta): beta_map[str(name)]=float(b)
    for v in numvars:
        x=pd.to_numeric(df[v],errors='coerce'); m=np.isfinite(x)&np.isfinite(yall)
        n=int(m.sum())
        if n<3: continue
        xv=x[m]; yv=yall[m]
        pear=float(xv.corr(yv,method='pearson')) if xv.nunique()>1 and yv.nunique()>1 else np.nan
        spear=float(xv.corr(yv,method='spearman')) if xv.nunique()>1 and yv.nunique()>1 else np.nan
        beta=beta_map.get(v,np.nan)
        strength=float(np.nanmax(np.abs([spear,beta]))) if np.any(np.isfinite([spear,beta])) else (abs(pear) if np.isfinite(pear) else np.nan)
        rows.append({'Variable':v,'VariableType':'Numeric','PearsonR':pear,'SpearmanR':spear,'StandardizedBeta':beta,'EtaSquared':np.nan,'EffectStrength':strength,'Samples':n})
    if include_categories:
        ymean=float(yall[np.isfinite(yall)].mean()) if np.isfinite(yall).any() else np.nan
        for v in catvars:
            tmp=pd.DataFrame({'g':df[v],'y':yall}).dropna(); n=len(tmp)
            if n<3 or tmp['g'].nunique()<2: continue
            total=float(((tmp.y-ymean)**2).sum())
            between=float(sum(len(g)*(float(g.y.mean())-ymean)**2 for _,g in tmp.groupby('g',dropna=False)))
            eta=between/total if total>1e-15 else np.nan
            rows.append({'Variable':v,'VariableType':'Categorical','PearsonR':np.nan,'SpearmanR':np.nan,'StandardizedBeta':np.nan,'EtaSquared':eta,'EffectStrength':eta,'Samples':n})
    if not rows: return pd.DataFrame(columns=['Variable','VariableType','PearsonR','SpearmanR','StandardizedBeta','EtaSquared','EffectStrength','Samples'])
    return pd.DataFrame(rows).sort_values('EffectStrength',ascending=False,na_position='last').reset_index(drop=True)


def _robust_losses(df: pd.DataFrame, objectives: list[ObjectiveSpec], relative_uncertainty_pct: float, k_sigma: float):
    rel=max(0.0,float(relative_uncertainty_pct))/100.0; k=max(0.0,float(k_sigma)); losses=[]; robust_values={}; valid=np.ones(len(df),dtype=bool)
    for o in objectives:
        if o.metric not in df.columns: raise KeyError(f"Objective metric {o.metric!r} is not in the result table")
        x=pd.to_numeric(df[o.metric],errors='coerce').to_numpy(float); sigma=np.abs(x)*rel; goal=str(o.goal).lower()
        valid &= np.isfinite(x)
        if goal.startswith('max'):
            rv=x-k*sigma; loss=-rv
        elif goal.startswith('min'):
            rv=x+k*sigma; loss=rv
        elif goal.startswith('target'):
            if o.target is None: raise ValueError(f'Target objective {o.metric!r} requires a target')
            t=float(o.target); scale=abs(t) if abs(t)>1e-15 else 1.0; rv=x
            loss=(np.abs(x-t)+k*sigma)/scale
        else: raise ValueError(f'Unknown objective goal {o.goal!r}')
        robust_values[o.metric]=rv; losses.append(loss)
    return losses,robust_values,valid


def robust_optimize_batch_designs(df: pd.DataFrame, objectives: list[ObjectiveSpec], constraints: list[ConstraintSpec] | None=None,
                                  relative_uncertainty_pct: float=5.0, k_sigma: float=1.0,
                                  manufacturing_rules: ManufacturingRules | None=None,
                                  robust_constraints: bool=True) -> pd.DataFrame:
    """Conservative k-sigma optimisation over the solved batch database.

    The user supplies an assumed relative response uncertainty.  For a maximise
    objective the robust value is x-k*sigma; for minimise it is x+k*sigma.  This
    is a screening model, not a stochastic re-analysis of the RVE solver.
    """
    if not objectives: raise ValueError('At least one objective is required')
    out=enrich_results(df.copy()); rel=max(0.0,float(relative_uncertainty_pct))/100.0; k=max(0.0,float(k_sigma))
    feasible=np.ones(len(out),dtype=bool)
    for c in constraints or []:
        if c.metric not in out.columns: raise KeyError(f"Constraint metric {c.metric!r} is not in the result table")
        x=pd.to_numeric(out[c.metric],errors='coerce').to_numpy(float); sigma=np.abs(x)*rel if robust_constraints else np.zeros_like(x); r=str(c.relation).lower(); m=np.isfinite(x)
        if r in ('<=','≤','max'): m &= (x+k*sigma)<=float(c.value)
        elif r in ('>=','≥','min'): m &= (x-k*sigma)>=float(c.value)
        elif r in ('between','range'):
            if c.upper is None: raise ValueError(f'Between constraint {c.metric!r} requires upper')
            lo,hi=sorted((float(c.value),float(c.upper))); m &= (x-k*sigma)>=lo; m &= (x+k*sigma)<=hi
        elif r in ('=','==','target'):
            tol=abs(float(c.upper)) if c.upper is not None else 0.0; m &= (np.abs(x-float(c.value))+k*sigma)<=tol
        else: raise ValueError(f'Unknown constraint relation {c.relation!r}')
        feasible &= m
    mfg=manufacturing_feasible_mask(out,manufacturing_rules); out['ManufacturingFeasible']=mfg; feasible &= mfg
    losses,robust_values,finite=_robust_losses(out,objectives,relative_uncertainty_pct,k_sigma); feasible &= finite
    idx=np.where(feasible)[0]; pareto=np.zeros(len(out),dtype=bool); score=np.full(len(out),np.nan)
    if len(idx):
        L=np.column_stack(losses)[idx]; keep=np.ones(len(idx),dtype=bool)
        for i in range(len(idx)):
            dominated=np.all(L<=L[i],axis=1)&np.any(L<L[i],axis=1); dominated[i]=False
            if np.any(dominated): keep[i]=False
        pareto[idx[keep]]=True
        weights=np.array([max(0.0,float(o.weight)) for o in objectives]);
        if not np.any(weights>0): weights[:]=1.0
        weights=weights/weights.sum(); util=[]
        for loss in losses:
            v=loss[idx]; lo=float(np.min(v)); hi=float(np.max(v)); u=np.ones_like(v) if hi-lo<=1e-15 else 1-(v-lo)/(hi-lo); util.append(np.clip(u,0,1))
        score[idx]=np.column_stack(util)@weights
    out['RobustOptimizationFeasible']=feasible; out['RobustOptimizationPareto']=pareto; out['RobustOptimizationScore']=score
    out['AssumedResponseUncertaintyPct']=float(relative_uncertainty_pct); out['RobustKSigma']=float(k_sigma)
    for metric,rv in robust_values.items(): out[f'RobustValue_{metric}']=rv
    return out
