from pathlib import Path
from rvestudio.database import generate_parametric_database,run_batch


def test_database_generation_and_batch(tmp_path:Path):
    inp=tmp_path/'db.csv'; out=tmp_path/'res.csv'
    db=generate_parametric_database(inp,num_cases=8,seed=123)
    assert len(db)==8
    res=run_batch(inp,out,workers=1,case_limit=4)
    assert len(res)==4
    assert (res.Status=='OK').all()


def test_batch_contains_full_ch_and_cost(tmp_path:Path):
    inp=tmp_path/'db_cost.csv'; out=tmp_path/'res_cost.csv'
    generate_parametric_database(inp,num_cases=2,seed=7)
    res=run_batch(inp,out,workers=1,case_limit=1)
    row=res.iloc[0]
    for c in ['C11','C12','C66','TotalMaterialCost','EffectiveCostPerKg','CostPerRVEVolume']:
        assert c in res.columns
    assert row['TotalMaterialCost'] > 0
