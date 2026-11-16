import numpy as np
import pandas as pd
from rvestudio.postprocess import enrich_results, fit_power_law, pareto_mask_min_x_max_y, directional_youngs_modulus, merge_input_metadata


def test_enrich_specific_and_anisotropy():
    df=pd.DataFrame({
        'CaseID':['A'], 'E1':[100.0], 'E2':[25.0], 'E3':[50.0], 'G12':[10.0], 'G13':[8.0], 'G23':[4.0],
        'Q11':[120.0], 'Q22':[30.0], 'Q12':[15.0], 'Q66':[10.0], 'Q44':[8.0], 'Q55':[4.0],
        'RelativeDensity':[0.1], 'AverageDensityKgM3':[200.0], 'LengthA_mm':[60.0], 'LengthB_mm':[30.0]
    })
    out=enrich_results(df)
    assert np.isclose(out.loc[0,'SpecificE1_GPa_per_Mg_m3'],0.5)
    assert np.isclose(out.loc[0,'E12_Anisotropy'],4.0)
    assert np.isclose(out.loc[0,'GeometryAspectRatio'],2.0)
    assert np.isclose(out.loc[0,'Q12_CouplingIndex'],15/np.sqrt(120*30))


def test_power_law_recovers_exponent():
    x=np.array([.02,.04,.08,.16])
    y=3.5*x**2.2
    fit=fit_power_law(x,y)
    assert fit is not None
    assert np.isclose(fit.exponent,2.2,atol=1e-10)
    assert np.isclose(fit.coefficient,3.5,atol=1e-10)
    assert fit.r_squared > .999999


def test_pareto_and_directional_modulus():
    m=pareto_mask_min_x_max_y([1,2,3,4],[1,3,2,5])
    assert m.tolist()==[True,True,False,True]
    ang,E=directional_youngs_modulus(100,25,10,.2,[0,90])
    assert np.isclose(E[0],100)
    assert np.isclose(E[1],25)


def test_merge_old_results_with_input_metadata():
    res=pd.DataFrame({'CaseID':['C1'],'E1':[1.]})
    inp=pd.DataFrame({'CaseID':['C1'],'GeometryID':['HEX'],'LengthA_mm':[50.]})
    out=merge_input_metadata(res,inp)
    assert out.loc[0,'GeometryID']=='HEX'
    assert out.loc[0,'LengthA_mm']==50.

from rvestudio.postprocess import ObjectiveSpec, ConstraintSpec, optimize_batch_designs


def test_multiobjective_pareto_and_target_selection():
    df=pd.DataFrame({
        'CaseID':['A','B','C','D'],
        'E1':[100.,130.,115.,90.],
        'MassKg':[1.0,1.8,1.1,0.8],
        'TotalMaterialCost':[10.,30.,12.,8.],
        'Q11':[50.,70.,60.,40.],
    })
    obj=[ObjectiveSpec('E1','max',weight=1.0),ObjectiveSpec('MassKg','min',weight=1.0),ObjectiveSpec('TotalMaterialCost','min',weight=.5)]
    out=optimize_batch_designs(df,obj,[ConstraintSpec('Q11','>=',45.)])
    assert not bool(out.loc[out.CaseID=='D','OptimizationFeasible'].iloc[0])
    assert out.OptimizationPareto.sum() >= 1
    assert np.isfinite(out.loc[out.OptimizationFeasible,'OptimizationScore']).all()

    target=optimize_batch_designs(df,[ObjectiveSpec('E1','target',target=115.,weight=1.)])
    best=target.loc[target.OptimizationScore.idxmax(),'CaseID']
    assert best=='C'
    assert np.isclose(target.loc[target.CaseID=='C','TargetDeviation_E1'].iloc[0],0.0)

from rvestudio.postprocess import (
    ManufacturingRules, manufacturing_feasible_mask, sensitivity_analysis,
    robust_optimize_batch_designs, design_variable_columns,
)


def test_manufacturing_metrics_and_rules():
    df=pd.DataFrame({
        'CaseID':['A','B'], 'T1_mm':[0.5,0.1], 'T2_mm':[0.75,0.8],
        'LengthA_mm':[50.,100.], 'LengthB_mm':[50.,25.], 'RelativeDensity':[0.08,0.2],
        'LayupAngles':['45/-45/-45/45','0/45/90'], 'E1':[100.,120.]
    })
    e=enrich_results(df)
    assert np.isclose(e.loc[0,'MinWallThickness_mm'],0.5)
    assert np.isclose(e.loc[0,'WallThicknessRatio'],1.5)
    assert bool(e.loc[0,'SymmetricLayup'])
    assert bool(e.loc[0,'BalancedLayup'])
    rules=ManufacturingRules(min_wall_thickness_mm=.25,max_wall_thickness_ratio=2.0,max_geometry_aspect_ratio=2.0,max_relative_density=.1,require_symmetric_layup=True)
    m=manufacturing_feasible_mask(e,rules)
    assert m.tolist()==[True,False]


def test_sensitivity_recovers_dominant_variable_and_category_effect():
    n=30
    x=np.linspace(0,1,n); z=np.sin(np.linspace(0,8,n)); group=np.where(np.arange(n)<15,'A','B')
    y=10*x+0.05*z+(group=='B')*2.0
    df=pd.DataFrame({'Param1Value':x,'T1_mm':z+2,'GeometryID':group,'E1':y})
    num,cat=design_variable_columns(df)
    assert 'Param1Value' in num and 'T1_mm' in num and 'GeometryID' in cat
    s=sensitivity_analysis(df,'E1')
    top=s.iloc[0]['Variable']
    assert top in ('Param1Value','GeometryID')
    row=s[s.Variable=='Param1Value'].iloc[0]
    assert row.SpearmanR > .9


def test_robust_optimisation_penalises_uncertainty_and_manufacturing():
    df=pd.DataFrame({
        'CaseID':['A','B','C'], 'E1':[100.,110.,90.], 'MassKg':[1.0,1.2,.8],
        'T1_mm':[.5,.15,.4], 'T2_mm':[.5,.15,.4], 'LengthA_mm':[50.,50.,50.], 'LengthB_mm':[50.,50.,50.],
        'RelativeDensity':[.08,.08,.08], 'LayupAngles':['0/90/90/0']*3,
    })
    obj=[ObjectiveSpec('E1','max',weight=1),ObjectiveSpec('MassKg','min',weight=.5)]
    rules=ManufacturingRules(min_wall_thickness_mm=.25)
    out=robust_optimize_batch_designs(df,obj,relative_uncertainty_pct=10,k_sigma=1,manufacturing_rules=rules)
    assert not bool(out.loc[out.CaseID=='B','RobustOptimizationFeasible'].iloc[0])
    assert np.isclose(out.loc[out.CaseID=='A','RobustValue_E1'].iloc[0],90.0)
    assert np.isfinite(out.loc[out.RobustOptimizationFeasible,'RobustOptimizationScore']).all()

def test_buckling_factor_can_be_manufacturing_constraint():
    df=pd.DataFrame({'CaseID':['A','B'],'E1':[100.,120.],'BucklingLambdaMin':[1.8,1.2]})
    mask=manufacturing_feasible_mask(df,ManufacturingRules(min_buckling_factor=1.5))
    assert mask.tolist()==[True,False]
