from __future__ import annotations
import math,os
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

from .geometry import BUILTIN_GEOMETRIES
from .plugins import discover_geometry_plugins
from .library import load_material_library,load_geometry_library,load_layup_library,material_mapping,PACKAGE_ROOT,ensure_library_dirs
from .models import Ply,Laminate,Element,RVEModel
from .homogenization import homogenize_rve
from .reliability import assess as assess_reliability
from .validation import extract_six
from .postprocess import enrich_results
from .load_recovery import recover_wall_loads
from .buckling import evaluate_wall_buckling


def lhs(n:int,d:int,rng:np.random.Generator)->np.ndarray:
    X=np.empty((n,d))
    for j in range(d):
        perm=rng.permutation(n); X[:,j]=(perm+rng.random(n))/n
    return X


def _read_table(path:str|Path)->pd.DataFrame:
    path=Path(path); suf=path.suffix.lower()
    if suf=='.csv': return pd.read_csv(path)
    if suf in ('.xlsx','.xls'):
        xf=pd.ExcelFile(path); sheet='ParametricCases' if 'ParametricCases' in xf.sheet_names else ('Results' if 'Results' in xf.sheet_names else xf.sheet_names[0]); return pd.read_excel(path,sheet_name=sheet)
    if suf=='.parquet': return pd.read_parquet(path)
    raise ValueError('Supported database formats: .csv, .xlsx, .parquet')


def _write_table(df:pd.DataFrame,path:str|Path,sheet='ParametricCases')->None:
    path=Path(path); suf=path.suffix.lower()
    if suf=='.csv': df.to_csv(path,index=False)
    elif suf in ('.xlsx','.xls'): df.to_excel(path,index=False,sheet_name=sheet)
    elif suf=='.parquet':
        try: df.to_parquet(path,index=False)
        except ImportError as e: raise ImportError('Parquet requires pyarrow or fastparquet. Install RVE-Studio[parquet].') from e
    else: raise ValueError('Supported output formats: .csv, .xlsx, .parquet')


def generate_parametric_database(output_file:str|Path,num_cases:int=5000,seed:int=20260906,
        relative_density_range=(0.03,0.15),height:float=150.0,max_walls:int=20,
        materials_rows:list[dict]|None=None,geometry_rows:list[dict]|None=None,layup_rows:list[dict]|None=None,
        plugin_dir:str|Path|None=None,wall_shear:str='euler-bernoulli',bending:str|None=None,twist_release=None)->pd.DataFrame:
    mats=[r for r in (materials_rows or load_material_library()) if r.get('active',True)]
    geoms=[r for r in (geometry_rows or load_geometry_library()) if r.get('active',True)]
    lays=[r for r in (layup_rows or load_layup_library()) if r.get('active',True)]
    if not mats or not geoms or not lays: raise ValueError('At least one active material, geometry and layup are required')
    from .walls import is_timoshenko
    _beam_name='timoshenko' if is_timoshenko(wall_shear) else 'euler-bernoulli'
    from .validation import bending_option
    _bend=bending_option(bending) if bending else None      # None: the solver of each geometry definition
    _twist=None if twist_release is None else [float(x) for x in twist_release]
    registry=discover_geometry_plugins(Path(plugin_dir) if plugin_dir else PACKAGE_ROOT/'plugins')
    rng=np.random.default_rng(seed); weights=np.array([float(g.get('weight',1.0)) for g in geoms]); ideal=num_cases*weights/weights.sum(); counts=np.floor(ideal).astype(int)
    rem=num_cases-counts.sum(); order=np.argsort(-(ideal-counts)); counts[order[:rem]]+=1
    rows=[]; case_no=0
    for gd,n in zip(geoms,counts):
        if n<=0: continue
        gid=str(gd['id']).upper()
        if gid not in registry: raise ValueError(f'Geometry plugin {gid} not found')
        params=gd.get('parameters',[]); d=max(1,len(params))+1+max_walls; X=lhs(n,d,rng)
        for i in range(n):
            case_no+=1; pvals={}
            for k,p in enumerate(params): pvals[p['name']]=float(p['min'])+(float(p['max'])-float(p['min']))*X[i,k]
            geom=registry[gid](pvals); nw=len(geom.connectivity)
            if nw>max_walls: raise ValueError(f'{gid} has {nw} walls > max_walls={max_walls}')
            rho=relative_density_range[0]+(relative_density_range[1]-relative_density_range[0])*X[i,max(1,len(params))]
            factors=0.65+(1.35-0.65)*X[i,max(1,len(params))+1:max(1,len(params))+1+max_walls]
            lengths=np.array([np.linalg.norm(geom.nodes[j]-geom.nodes[ii]) for ii,j in geom.connectivity]); area=abs(np.linalg.det(np.column_stack((geom.a1,geom.a2))))
            scale=rho*area/np.sum(geom.wall_weights*lengths*factors[:nw]); T=np.full(max_walls,np.nan); T[:nw]=scale*factors[:nw]
            m=mats[(case_no-1)%len(mats)]; lay=lays[(case_no-1)%len(lays)]; angles=lay['angles']
            row={'CaseID':f'C{case_no:06d}','Active':True,'GeometryID':gid,'GeometryType':gd.get('name',geom.name),'Height':height,
                 'Material':m['name'],'MaterialCostPerKg':m.get('cost_per_kg',np.nan),'LayupID':lay['id'],'LayupAngles':'/'.join(str(a) for a in angles),'TargetRelativeDensity':rho,'NumWalls':nw,
                 'a1x':geom.a1[0],'a1y':geom.a1[1],'a2x':geom.a2[0],'a2y':geom.a2[1],'SolverType':_bend[0] if _bend else gd.get('solver','paper'),'WallShear':_beam_name,
                 'StripMode':_bend[1] if _bend else gd.get('strip_mode','cylindrical'),
                 'TwistA':_twist[0] if _twist else np.nan,'TwistB':_twist[1] if _twist else np.nan,'TwistN':_twist[2] if _twist else np.nan,'PeriodicTol':1e-8,'AnchorNode':1,'LengthUnit':'mm','Gravity_m_s2':9.80665,'GeometryNote':geom.note}
            for k in range(4):
                if k<len(params): row[f'Param{k+1}Name']=params[k]['name']; row[f'Param{k+1}Value']=pvals[params[k]['name']]
                else: row[f'Param{k+1}Name']=''; row[f'Param{k+1}Value']=np.nan
            row['LengthA_mm']=geom.length_a; row['LengthB_mm']=geom.length_b; row['LengthRatio']=geom.length_a/geom.length_b if np.isfinite(geom.length_b) and geom.length_b!=0 else np.nan
            for k in range(max_walls): row[f'T{k+1}_mm']=T[k]
            rows.append(row)
    df=pd.DataFrame(rows)
    out_path=Path(output_file)
    if out_path.suffix.lower() in ('.xlsx','.xls'):
        with pd.ExcelWriter(out_path) as writer:
            df.to_excel(writer,index=False,sheet_name='ParametricCases')
            pd.DataFrame(mats).to_excel(writer,index=False,sheet_name='Materials')
            pd.DataFrame(geoms).to_excel(writer,index=False,sheet_name='GeometryLibrary')
            pd.DataFrame(lays).to_excel(writer,index=False,sheet_name='LayupLibrary')
    else:
        _write_table(df,out_path)
    return df


def _row_twist(row):
    try:
        c=[float(row.get(k)) for k in ('TwistA','TwistB','TwistN')]
    except (TypeError,ValueError):
        return None
    return tuple(c) if all(math.isfinite(x) for x in c) else None


_LIB_CACHE: dict = {}


def _library_context():
    """Material, geometry and layup libraries, read once per process and re-read only when a file changes."""
    from .library import DATA_DIR
    files=[DATA_DIR/'materials.json',DATA_DIR/'geometries.json',DATA_DIR/'layups.json']
    sig=tuple((f.stat().st_mtime_ns if f.exists() else 0) for f in files)
    if _LIB_CACHE.get('sig')!=sig:
        mats_rows=load_material_library(); _LIB_CACHE.update(sig=sig,value=(mats_rows,material_mapping(mats_rows),load_geometry_library(),load_layup_library()))
    return _LIB_CACHE['value']


def build_case_from_row(row:pd.Series|dict,materials_rows=None,geometry_rows=None,layup_rows=None,plugin_dir=None):
    if isinstance(row,pd.Series): row=row.to_dict()
    if materials_rows is None and geometry_rows is None and layup_rows is None:
        mats_rows,mats,geoms,lays=_library_context()
    else:
        mats_rows=materials_rows or load_material_library(); mats=material_mapping(mats_rows); geoms=geometry_rows or load_geometry_library(); lays=layup_rows or load_layup_library()
    gd=next((g for g in geoms if str(g['id']).upper()==str(row['GeometryID']).upper()),None)
    if gd is None: raise KeyError(f"GeometryID {row['GeometryID']} not in library")
    registry=discover_geometry_plugins(Path(plugin_dir) if plugin_dir else PACKAGE_ROOT/'plugins'); gid=str(gd['id']).upper(); pvals={}
    for k,p in enumerate(gd.get('parameters',[]),1):
        val=row.get(f'Param{k}Value',np.nan); pvals[p['name']]=float(val) if pd.notna(val) else 0.5*(float(p['min'])+float(p['max']))
    geom=registry[gid](pvals); global_layid=int(row['LayupID']); global_material=str(row['Material']); nw=len(geom.connectivity)
    laminates=[]; elements=[]

    def _nonblank(value, default=None):
        if value is None: return default
        try:
            if pd.isna(value): return default
        except Exception:
            pass
        text=str(value).strip()
        return text if text else default

    def _angles_from(value, fallback):
        text=_nonblank(value)
        if text is None: return [float(a) for a in fallback]
        vals=[x for x in str(text).replace(',', '/').replace(';','/').split('/') if x.strip()]
        return [float(x) for x in vals] if vals else [float(a) for a in fallback]

    for e,(ii,jj) in enumerate(geom.connectivity):
        wi=e+1
        # Backward-compatible global Material/LayupID plus optional per-wall overrides.
        # Accepted columns: W1Material, W1LayupID, W1LayupAngles, ...
        material=_nonblank(row.get(f'W{wi}Material'), global_material)
        layid_raw=_nonblank(row.get(f'W{wi}LayupID'), global_layid)
        layid=int(float(layid_raw)); lay=next(l for l in lays if int(l['id'])==layid)
        angles=_angles_from(row.get(f'W{wi}LayupAngles'), lay['angles'])
        t=float(row[f'T{wi}_mm'])
        if not math.isfinite(t) or t<=0: raise ValueError(f'Wall W{wi} thickness must be positive')
        if not angles: raise ValueError(f'Wall W{wi} has no laminate ply angles')
        plies=[Ply(material,float(a),t/len(angles)) for a in angles]
        laminates.append(Laminate(plies,name=f'{material} wall {wi}')); elements.append(Element(int(ii),int(jj),e,float(geom.wall_weights[e])))
    model=RVEModel(geom.nodes,elements,geom.a1,geom.a2,float(row.get('Height',150.0)),periodic_tol=float(row.get('PeriodicTol',1e-8)),anchor_node=int(row.get('AnchorNode',1))-1,
                   wall_element_type=str(row.get('SolverType','paper')),strip_bending_mode=str(row.get('StripMode','cylindrical')),length_unit=str(row.get('LengthUnit','mm')),gravity=float(row.get('Gravity_m_s2',9.80665)),name=str(row.get('CaseID','Case')),wall_groups=list(geom.wall_groups),wall_shear=str(row.get('WallShear','timoshenko') if str(row.get('WallShear','')).strip() not in ('', 'nan') else 'timoshenko'),twist_release=_row_twist(row))
    return model,mats,laminates


def _result_metadata(row_dict:dict)->dict:
    # Keep design descriptors beside solved properties so Results Explorer can
    # filter/colour/rank designs without reopening the input database.  v1.3.4
    # also preserves wall thickness columns for manufacturability/sensitivity.
    wanted=['GeometryID','GeometryType','Material','MaterialCostPerKg','LayupID','LayupAngles','TargetRelativeDensity','Height','NumWalls',
            'LengthA_mm','LengthB_mm','LengthRatio','SolverType','StripMode','WallShear','TwistA','TwistB','TwistN',
            'Param1Name','Param1Value','Param2Name','Param2Value','Param3Name','Param3Value','Param4Name','Param4Value']
    out={k:row_dict.get(k,np.nan) for k in wanted if k in row_dict}
    for k,v in row_dict.items():
        ks=str(k)
        if ks.startswith('T') and ks.endswith('_mm') and ks[1:-3].isdigit(): out[ks]=v
        elif ks.startswith('W') and any(tag in ks for tag in ('Material','LayupID','LayupAngles')): out[ks]=v
    return out


def solve_database_case(row_dict:dict,buckling_settings:dict|None=None)->dict:
    cid=str(row_dict.get('CaseID','')); meta=_result_metadata(row_dict)
    try:
        model,mats,lams=build_case_from_row(row_dict); out=homogenize_rve(model,mats,lams); six=extract_six(out)
        p=out.props
        ch={f'C{i+1}{j+1}':float(out.CH[i,j]) for i in range(6) for j in range(i,6)}
        solved={'CaseID':cid,**meta,'Status':'OK','Error':'',**ch,'Q11':six[0],'Q22':six[1],'Q12':six[2],'Q66':six[3],'Q44':six[4],'Q55':six[5],
                'E1':p.E1,'E2':p.E2,'E3':p.E3,'G12':p.G12,'G13':p.G13,'G23':p.G23,'nu12':p.nu12,'nu13':p.nu13,'nu23':p.nu23,
                'RelativeDensity':out.physical.relative_density,'MassKg':out.physical.total_mass_kg,'WeightN':out.physical.weight_n,'AverageDensityKgM3':out.physical.average_density_kg_m3,
                'TotalMaterialCost':out.physical.total_material_cost,'EffectiveCostPerKg':out.physical.effective_cost_per_kg,'CostPerRVEVolume':out.physical.cost_per_rve_m3}
        rel=assess_reliability(out)
        solved.update({'PsiMax':max((abs(w.psi) for w in rel.walls),default=0.0),'BendTwistIndexMax':max((w.bend_twist_index for w in rel.walls),default=0.0),'ChiBMax':max((w.chi_B for w in rel.walls),default=0.0),
                       'BetaMax':max((w.beta_max for w in rel.walls),default=0.0),'ReliabilityInPlane':rel.in_plane,'ReliabilityOutOfPlane':rel.out_of_plane})
        if buckling_settings:
            mode=str(buckling_settings.get('control_mode','stress')).lower()
            values=buckling_settings.get('macro_values',[-1.0,-1.0,0.0,0.0,0.0,0.25])
            rec=recover_wall_loads(model,out,macro_stress=values) if mode.startswith('stress') else recover_wall_loads(model,out,macro_strain=values)
            buck=evaluate_wall_buckling(model,out,rec,m_max=int(buckling_settings.get('m_max',8)),n_max=int(buckling_settings.get('n_max',8)),
                safety_factor=float(buckling_settings.get('safety_factor',1.5)),use_abd_condensation=bool(buckling_settings.get('use_abd_condensation',True)),
                convergence_check=bool(buckling_settings.get('convergence_check',True)))
            critical=next((w for w in buck.walls if int(w['wall'])==int(buck.critical_wall or -1)),None)
            lam=float(buck.minimum_lambda)
            solved.update({'BucklingLambdaMin':lam,'BucklingCriticalWall':buck.critical_wall or np.nan,'BucklingPass':buck.passed,
                'BucklingSafetyFactor':buck.safety_factor,'BucklingUtilization':(1.0/lam if np.isfinite(lam) and lam>0 else 0.0),
                'BucklingDesignUtilization':(buck.safety_factor/lam if np.isfinite(lam) and lam>0 else 0.0),
                'BucklingMmax':int(buck.settings['m_max']),'BucklingNmax':int(buck.settings['n_max']),
                'BucklingBoundaryCondition':'SSSS','BucklingTheory':'CLPT multi-mode Ritz/Navier','BucklingABDCondensed':bool(buck.settings['use_abd_condensation']),
                'BucklingMaxConvergenceError':max([float(w['mode_convergence_relative_error']) for w in buck.walls if np.isfinite(w['mode_convergence_relative_error'])],default=np.nan)})
            if critical:
                solved.update({'BucklingCriticalN11':critical['N11_N_per_mm'],'BucklingCriticalN33':critical['N33_N_per_mm'],'BucklingCriticalN13':critical['N13_N_per_mm'],
                    'BucklingCriticalModeM':critical['dominant_m'],'BucklingCriticalModeN':critical['dominant_n'],
                    'BucklingCriticalN11AtBifurcation':critical['N11_N_per_mm']*lam if np.isfinite(lam) else np.nan,
                    'BucklingCriticalN33AtBifurcation':critical['N33_N_per_mm']*lam if np.isfinite(lam) else np.nan,
                    'BucklingCriticalN13AtBifurcation':critical['N13_N_per_mm']*lam if np.isfinite(lam) else np.nan})
            stress_labels=('11','22','33','23','13','12')
            for lab,v in zip(stress_labels,rec.macro_stress): solved[f'BucklingSigma{lab}_MPa' if lab in ('11','22','33') else f'BucklingTau{lab}_MPa']=float(v)
            for lab,v in zip(stress_labels,rec.macro_strain): solved[f'BucklingEps{lab}' if lab in ('11','22','33') else f'BucklingGamma{lab}']=float(v)
        return solved
    except Exception as exc:
        return {'CaseID':cid,**meta,'Status':'FAILED','Error':str(exc)}


def run_batch(input_file:str|Path,output_file:str|Path,workers:int=1,case_limit:int|None=None,progress=None,buckling_settings:dict|None=None)->pd.DataFrame:
    df=_read_table(input_file); active=df[df.get('Active',True).astype(bool)] if 'Active' in df else df
    if case_limit: active=active.head(case_limit)
    rows=active.to_dict('records'); results=[]
    if workers<=1:
        for i,r in enumerate(rows,1):
            results.append(solve_database_case(r,buckling_settings));
            if progress: progress(i,len(rows),results[-1])
    else:
        # Cases are sent to the workers in chunks (far less inter-process traffic than one task per case);
        # results come back in input order.
        chunk=max(1,min(256,len(rows)//(workers*16) or 1))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for i,res in enumerate(ex.map(solve_database_case,rows,[buckling_settings]*len(rows),chunksize=chunk),1):
                results.append(res)
                if progress: progress(i,len(rows),res)
    out=enrich_results(pd.DataFrame(results)); _write_table(out,output_file,sheet='Results'); return out
