from pathlib import Path
import numpy as np

from rvestudio.abaqus import AbaqusSettings, generate_abaqus_input, prepare_abaqus_run, _extractor_source
from rvestudio.models import Material, Ply, Laminate, Element, RVEModel


def _case():
    mats={"M": Material("M", 100000, 10000, 5000, 0.25, 4000, 1600)}
    lams=[Laminate([Ply("M",0,0.2),Ply("M",45,0.2)],"L")]
    model=RVEModel(np.array([[0.,0.],[10.,0.],[5.,8.]]),[Element(0,1,0),Element(1,2,0)],np.array([10.,0.]),np.array([0.,8.]),20.)
    return model,mats,lams


def test_abaqus_s4_input_contains_shells_layup_and_bc():
    model,mats,lams=_case(); s=AbaqusSettings(mesh_size=4,load_case="XZ shear")
    text=generate_abaqus_input(model,mats,lams,s,job_name="T")
    assert "*ELEMENT, TYPE=S4" in text
    assert "*SHELL SECTION, ELSET=WALL_1, COMPOSITE" in text
    assert "*ORIENTATION, NAME=ORI_W1" in text
    assert "*ELASTIC, TYPE=LAMINA" in text
    assert "TOP, 1, 1, 1" in text
    assert "BOTTOM, 1, 6, 0." in text


def test_abaqus_run_package_is_reproducible(tmp_path):
    model,mats,lams=_case(); s=AbaqusSettings(mesh_size=5,load_case="Z compression",displacement=0.5)
    files=prepare_abaqus_run(tmp_path,model,mats,lams,s)
    assert files.input_file.exists() and files.extractor_file.exists() and files.snapshot_file.exists()
    assert "TYPE=S4" in files.input_file.read_text()
    assert '"abaqus_settings"' in files.snapshot_file.read_text()


def test_abaqus_odb_extractor_is_2021_repository_safe():
    src = _extractor_source()
    assert ".nodeSets.get(" not in src
    assert "repo_get(assembly.nodeSets, 'TOP')" in src
    assert "'status': 'ERROR'" in src
    assert "traceback.format_exc()" in src
