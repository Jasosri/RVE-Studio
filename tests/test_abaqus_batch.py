import json
from pathlib import Path

import pandas as pd

from rvestudio.library import PACKAGE_ROOT
from rvestudio.database import _read_table, build_case_from_row
from rvestudio.examples_data import thesis_hexagonal_case
from rvestudio.abaqus_homogenization import PeriodicHomogenizationSettings, prepare_periodic_homogenization
from rvestudio.abaqus_batch import run_abaqus_batch, BACKEND_CH21


def test_database_row_accepts_per_wall_material_layup_and_angles():
    row = _read_table(PACKAGE_ROOT / 'data' / 'RVE_5000_Python_Input.csv').iloc[0].to_dict()
    row['W1Material'] = 'Boron/Epoxy'
    row['W1LayupID'] = 2
    row['W1LayupAngles'] = '0/90'
    model, mats, lams = build_case_from_row(row)
    assert model.elements[0].laminate == 0
    assert [p.material for p in lams[0].plies] == ['Boron/Epoxy', 'Boron/Epoxy']
    assert [p.angle_deg for p in lams[0].plies] == [0.0, 90.0]
    assert abs(sum(p.thickness for p in lams[0].plies) - float(row['T1_mm'])) < 1e-12


def test_periodic_package_resume_preserves_completed_case(tmp_path):
    model, mats, lams = thesis_hexagonal_case(1)
    settings = PeriodicHomogenizationSettings(mesh_size=10.0, strain_amplitude=0.001)
    files = prepare_periodic_homogenization(tmp_path, model, mats, lams, settings, run_name='Case_C000001')
    completed = {'status': 'COMPLETED', 'strain_energy_ALLSE': 1.25, 'macro_reaction_u1': {}, 'macro_displacement_u1': {}}
    files.cases[0].results_file.write_text(json.dumps(completed), encoding='utf-8')
    files2 = prepare_periodic_homogenization(tmp_path, model, mats, lams, settings, run_name='Case_C000001', resume=True)
    assert files2.run_dir == files.run_dir
    assert json.loads(files2.cases[0].results_file.read_text())['status'] == 'COMPLETED'


def test_abaqus_batch_checkpoints_and_skips_completed_output(tmp_path, monkeypatch):
    src = _read_table(PACKAGE_ROOT / 'data' / 'RVE_5000_Python_Input.csv').head(2)
    inp = tmp_path / 'input.csv'; out = tmp_path / 'results.csv'
    src.to_csv(inp, index=False)
    calls = []

    def fake_case(row, **kwargs):
        cid = str(row['CaseID']); calls.append(cid)
        return {'CaseID': cid, 'Status': 'OK', 'Backend': BACKEND_CH21, 'Error': '', 'C11': 123.0, 'RunDir': str(tmp_path/cid), 'Runtime_s': 0.1}

    import rvestudio.abaqus_batch as mod
    monkeypatch.setattr(mod, '_run_periodic_case', fake_case)
    df = run_abaqus_batch(inp, out, backend=BACKEND_CH21, abaqus_command='fake', run_root=tmp_path/'runs', resume=True)
    assert list(df['Status']) == ['OK', 'OK']
    assert len(calls) == 2
    assert out.exists()

    calls.clear()
    df2 = run_abaqus_batch(inp, out, backend=BACKEND_CH21, abaqus_command='fake', run_root=tmp_path/'runs', resume=True)
    assert len(calls) == 0
    assert list(df2['CaseID']) == list(src['CaseID'].astype(str))
