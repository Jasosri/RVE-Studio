"""The Help menu opens the user and theory guides; the citation list names the papers the software builds on."""
from rvestudio.guides import CITATIONS, GUIDES, guide_path


def test_both_guides_are_found_as_pdf():
    for kind in GUIDES:
        p = guide_path(kind)
        assert p is not None and p.suffix == ".pdf" and p.read_bytes()[:4] == b"%PDF"


def test_citations_name_the_current_and_both_earlier_papers():
    text = " ".join(t + " " + (u or "") for t, u in CITATIONS)
    assert "RVE Studio" in text
    assert "10.1016/j.engstruct.2023.117106" in text and "10.1016/j.compstruct.2024.118629" in text


def test_installed_application_keeps_libraries_in_the_user_folder(tmp_path):
    from rvestudio.library import USER_ROOT, library_root
    (tmp_path / "data").mkdir(); (tmp_path / "data" / "materials.json").write_text("[]")
    assert library_root(tmp_path, frozen=False) == tmp_path          # source checkout edits its own data/
    assert library_root(tmp_path, frozen=True) == USER_ROOT          # Windows application: never the install folder
    assert library_root(tmp_path / "missing", frozen=False) == USER_ROOT


def test_windows_build_files_are_consistent():
    from pathlib import Path
    import rvestudio
    root = Path(__file__).resolve().parents[1]
    iss = (root / "installer" / "RVE_Studio.iss").read_text()
    spec = (root / "RVE_Studio.spec").read_text()
    assert "name='RVE_Studio'" in spec and "installer/rve_studio.ico" in spec and (root / "installer" / "rve_studio.ico").is_file()
    assert "RVE_Studio.exe" in iss and "RVE_Studio_User_Guide.pdf" in iss and "RVE_Studio_Theory_Guide.pdf" in iss
    assert f"filevers=({rvestudio.__version__.replace('.', ', ')}, 0)" in (root / "installer" / "version_info.txt").read_text()
    for src in ("rvestudio/defaults", "rvestudio/gui/assets", "rvestudio/docs"):
        assert (root / src).is_dir() and f"('{src}', '{src}')" in spec
