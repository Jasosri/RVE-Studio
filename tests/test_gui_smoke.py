import os
import pytest
import numpy as np


def test_gui_main_window_constructs_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow

    app = QApplication.instance() or QApplication([])
    window = MainWindow()
    assert window.windowTitle().startswith("RVE Studio")
    assert window.single_page.nodes.rowCount() > 0
    assert window.single_page.walls.rowCount() > 0
    window.close()


def test_sketch_scene_survives_gc_and_can_draw_lattice_vectors():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import gc
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.sketch import RVESketchView

    app = QApplication.instance() or QApplication([])
    view = RVESketchView()
    gc.collect()
    assert view.scene() is not None
    assert view.scene() is view._scene
    view.set_lattice_vectors([10.0, 0.0], [-5.0, 8.6602540378])
    assert len(view._vector_items) == 6
    assert all(item.scene() is view._scene for item in view._vector_items)
    view.close()


def test_gui_all_pages_and_examples_construct_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow

    app=QApplication.instance() or QApplication([])
    window=MainWindow()
    # Construct/show every page to catch lazy widget/scene lifetime problems.
    for row in range(window.nav.count()):
        window.nav.setCurrentRow(row); app.processEvents()
        assert window.stack.currentIndex()==row
    # Load every supplied benchmark configuration through the real GUI path.
    for config in range(1,15):
        window.single_page.load_example_configuration(config); app.processEvents()
        summary=window.single_page.current_geometry_summary()
        assert summary.is_valid
        assert summary.detected_geometry_id=="HEX"
        assert summary.parameter_sources.get("L1")
        assert summary.parameter_sources.get("L2")
    # Exercise sketch dimension overlays and a simple editing constraint.
    page=window.single_page
    page.sketch.set_show_dimensions(True)
    assert all("mm" in w.label.text() for w in page.sketch.walls)
    assert page.sketch.constrain_wall(1,"vertical")
    app.processEvents()
    window.close()


def test_geometry_definition_wizard_constructs_preview_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.geometry_analysis import analyse_geometry
    from rvestudio.gui.geometry_wizard import GeometryDefinitionWizard

    app=QApplication.instance() or QApplication([])
    model,_,_=paper_hexagonal_case(1)
    summary=analyse_geometry(model.nodes,model.elements,model.a1,model.a2)
    dlg=GeometryDefinitionWizard(summary,existing_ids=['HEX','SQUARE','MIXED','TRI'])
    assert hasattr(dlg,'preview')
    assert dlg.preview.text()
    # Creating is now the default action, not silently reusing HEX.
    assert dlg.mode.currentData()=='scale_xy'
    assert dlg.id_edit.text().startswith('CUSTOM_HEX')
    dlg.close()


def test_existing_geometry_dropdown_loads_library_geometry_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow

    app=QApplication.instance() or QApplication([])
    window=MainWindow(); page=window.single_page
    ids=[page.geometry_selector.itemData(i) for i in range(page.geometry_selector.count())]
    for gid in ('HEX','SQUARE','MIXED','TRI'):
        assert gid in ids
    idx=page.geometry_selector.findData('HEX'); page.geometry_selector.setCurrentIndex(idx); app.processEvents()
    # L1/L2 appear as editable physical parameters in the selector row.
    assert page.geometry_parameter_labels[0].text()=='L1'
    assert page.geometry_parameter_labels[1].text()=='L2'
    page.geometry_parameter_spins[0].setValue(50.0); page.geometry_parameter_spins[1].setValue(50.0)
    page.load_selected_geometry(); app.processEvents()
    assert page.nodes.rowCount()==5
    assert page.walls.rowCount()==3
    summary=page.current_geometry_summary()
    assert summary.detected_geometry_id=='HEX'
    window.close()


def test_wall_load_recovery_gui_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow

    app=QApplication.instance() or QApplication([])
    window=MainWindow(); page=window.single_page
    page.wall_load_control.setCurrentText('Stress controlled [MPa]')
    vals=[-2.0,-1.0,-0.2,0.15,-0.1,0.5]
    for box,val in zip(page.wall_load_spins,vals): box.setValue(val)
    page.recover_current_wall_loads(); app.processEvents()
    assert page.wall_load_table.rowCount()==page.walls.rowCount()
    assert page.wall_load_validation.rowCount()>=7
    assert page.wall_load_df is not None and not page.wall_load_df.empty
    assert 'N11_N_per_mm' in page.wall_load_df.columns
    window.close()


def test_wall_buckling_gui_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow

    app=QApplication.instance() or QApplication([])
    window=MainWindow(); page=window.single_page
    page.wall_load_control.setCurrentText('Stress controlled [MPa]')
    vals=[-1.0,0.0,-1.0,0.0,0.25,0.0]
    for box,val in zip(page.wall_load_spins,vals): box.setValue(val)
    page.recover_current_wall_loads(); app.processEvents()
    page.buck_m.setValue(4); page.buck_n.setValue(4); page.buck_sf.setValue(1.0)
    page.evaluate_current_buckling(); app.processEvents()
    assert page.buckling_df is not None and not page.buckling_df.empty
    assert page.buck_table.rowCount()==page.walls.rowCount()
    assert 'lambda_cr' in page.buckling_df.columns
    assert np.isfinite(page.buckling_df['lambda_cr'].replace([np.inf,-np.inf],np.nan).dropna()).all()
    window.close()


def test_fe_page_shows_prepared_jobs(tmp_path):
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio import fe_workflow as fw
    from rvestudio.gui.main_window import MainWindow

    app = QApplication.instance() or QApplication([])
    fw.prepare_study(tmp_path / "s", paper_configs="1", settings=fw.FEStudySettings(meshes=[10.0], checks=["dbc"]))
    window = MainWindow()
    window.go_to("FE Design & Validation")
    page = window.fe_page
    assert window.current_page() is page
    page.root_edit.setText(str(tmp_path / "s"))
    page.reload()
    app.processEvents()
    assert page.jobs_table.rowCount() == 10            # six shell-term and four U3-restrained DBC cases
    assert "pending 10" in page.counts_label.text()
    assert page.project_combo.count() == 1
    window.go_to("Help")
    assert window.nav.currentItem().text() == "Help"
    window.close()


def test_fe_page_checks_the_current_design_without_abaqus(tmp_path):
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow

    app = QApplication.instance() or QApplication([])
    window = MainWindow()
    window.single_page.load_example_configuration(1)
    page = window.fe_page
    page.root_edit.setText(str(tmp_path / "study"))
    page.reload()
    page.add_current()
    page.meshes.setText("10")
    page.check_now()
    page.worker.wait(300000) if page.worker is not None else None
    for _ in range(50):
        app.processEvents()
    assert page.project_combo.count() == 1
    assert page.shell_table.rowCount() == 6 and page.solid_table.rowCount() == 6
    assert page.shell_table.item(0, 0).text() == "Q11" and page.solid_table.item(0, 0).text() == "C1111"
    assert page.matrix_tables["C"][1].rowCount() == 6 and page.matrix_tables["Q"][1].rowCount() == 5
    assert page.term_tables["eng"].rowCount() == 12 and page.term_tables["C"].item(0, 0).text() == "C1111"
    assert (tmp_path / "study" / "FE_STUDY_REPORT.html").exists()
    for i in range(page.def_source.count()):
        page.def_source.setCurrentIndex(i)
        app.processEvents()
    window.close()


def test_define_fe_rve_dialog_draws_and_centres_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from types import SimpleNamespace
    from PySide6.QtWidgets import QApplication
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.fe_rve import select_fe_rve
    from rvestudio.gui.fe_rve_dialog import FERVEDialog

    app = QApplication.instance() or QApplication([])
    model, _m, _l = paper_hexagonal_case(1)
    dlg = FERVEDialog(model)
    assert dlg.mode() == "auto" and dlg.current_spec() is None
    lay = dlg.layout_data
    rw, rh = lay["rect_width"], lay["rect_height"]
    ax = dlg.fig.axes[0]
    # Drag a sloppy 2.7 x 1.2-cell rectangle: snapped to 3 x 1 cells and centred on the geometry.
    dlg._press(SimpleNamespace(inaxes=ax, button=1, xdata=13.0, ydata=7.0))
    dlg._motion(SimpleNamespace(inaxes=ax, xdata=13 + 2.7 * rw, ydata=7 + 1.2 * rh))
    dlg._release(SimpleNamespace(inaxes=ax, xdata=13 + 2.7 * rw, ydata=7 + 1.2 * rh))
    app.processEvents()
    spec = dlg.current_spec()
    assert dlg.mode() == "window" and spec["cells"] == [3, 1]
    model.fe_rve = spec
    assert not select_fe_rve(model)[1]["selection"]["warnings"]
    dlg.rb_cells.setChecked(True)
    dlg.n1.setValue(3); dlg.n2.setValue(2)
    assert dlg.current_spec() == {"mode": "cells", "cells": [3, 2], "centre_cell": 0}
    dlg.accept()
    assert dlg.spec()["cells"] == [3, 2]
    dlg.close()


def test_fe_page_keeps_a_rve_choice_per_design_when_qt_is_available(tmp_path, monkeypatch):
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.io import load_project
    from rvestudio.gui import fe_rve_dialog
    from rvestudio.gui.main_window import MainWindow

    app = QApplication.instance() or QApplication([])
    window = MainWindow()
    window.single_page.load_example_configuration(1)
    page = window.fe_page
    page.add_current()
    monkeypatch.setattr(fe_rve_dialog.FERVEDialog, "exec", lambda self: True)
    monkeypatch.setattr(fe_rve_dialog.FERVEDialog, "spec", lambda self: {"mode": "cells", "cells": [3, 2]})
    page.define_rve()
    item = page.project_list.item(0)
    assert "3 x 2 cells" in item.text()
    path = page._project_path(item)
    assert load_project(path)[0].fe_rve == {"mode": "cells", "cells": [3, 2]}
    page.reset_rve()
    assert "automatic" in item.text() and load_project(page._project_path(item))[0].fe_rve is None
    window.close()


def test_geometry_library_shows_parameter_diagrams_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow, GeometryLibraryPage

    app = QApplication.instance() or QApplication([])
    window = MainWindow()
    page = next(window.page_at(i) for i in range(window.stack.count())
                if isinstance(window.page_at(i), GeometryLibraryPage))
    for r in range(page.table.rowCount()):
        page.table.setCurrentCell(r, 1)
        app.processEvents()
        assert page.diagram_fig.axes or "No diagram" in page.diagram_caption.text()
    page.table.setCurrentCell(0, 1)
    assert "<b>L1</b>" in page.diagram_caption.text()
    window.close()


def test_wall_editor_reverses_and_returns_rows_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QDialogButtonBox
    from rvestudio.gui.wall_editor import WallEditorDialog

    app = QApplication.instance() or QApplication([])
    nodes = [[0, 0], [-43.3, 75], [86.6, 0], [43.3, 75], [0, 50]]
    dlg = WallEditorDialog(nodes, [[1, 5, 1.0, 1, ''], [2, 5, 1.0, 1, ''], [4, 5, 1.0, 1, 'x']], [86.6, 0], [-43.3, 75])
    assert dlg.table.rowCount() == 3 and dlg.table.item(0, 6).text() == "90"
    dlg.select_rows([0, 2])
    dlg.reverse_selected()
    assert dlg.rows() == [[5, 1, 1.0, 1, ''], [2, 5, 1.0, 1, ''], [5, 4, 1.0, 1, 'x']]
    assert dlg.table.item(0, 6).text() == "-90"
    dlg.table.item(1, 1).setText("9")                                  # no node 9: blocked
    app.processEvents()
    assert dlg.problems() and not dlg.buttons.button(QDialogButtonBox.StandardButton.Ok).isEnabled()
    dlg.table.item(1, 1).setText("2")
    app.processEvents()
    assert not dlg.problems()
    dlg.close()


def test_wall_table_reverse_and_editor_round_trip_when_qt_is_available(monkeypatch):
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui import wall_editor
    from rvestudio.gui.main_window import MainWindow

    app = QApplication.instance() or QApplication([])
    window = MainWindow()
    page = window.single_page
    page.load_example_configuration(1)
    assert page._wall_table_rows()[0][:2] == [1, 5]
    assert all(w.arrow.isVisible() for w in page.sketch.walls)
    page.walls.selectRow(0)
    page.reverse_selected_walls()
    app.processEvents()
    assert page._wall_table_rows()[0][:2] == [5, 1]
    assert page.parse_model()[0].elements[0].node_i == 4                 # the solver sees Node I = N5
    monkeypatch.setattr(wall_editor.WallEditorDialog, "exec", lambda self: True)
    monkeypatch.setattr(wall_editor.WallEditorDialog, "rows", lambda self: [[1, 5, 1.0, 1, ''], [2, 5, 1.0, 1, ''], [4, 5, 1.0, 1, '']])
    page.open_wall_editor()
    assert page._wall_table_rows()[0][:2] == [1, 5]
    page.show_dir.setChecked(False)
    assert not any(w.arrow.isVisible() for w in page.sketch.walls)
    window.close()


def test_dashboard_tiles_and_table_editor_dialog_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QLabel, QPushButton
    from rvestudio.gui.main_window import MainWindow, DashboardPage
    from rvestudio.gui.ui_helpers import asset

    app = QApplication.instance() or QApplication([])
    window = MainWindow(); window.show(); app.processEvents()
    dash = window.page_at(0)
    assert isinstance(dash, DashboardPage)
    texts = [w.text() for w in dash.findChildren(QLabel)]
    assert not any("Paper benchmark" in t or "14 configurations" in t for t in texts)
    for img, *_ in DashboardPage.TILES:
        assert asset(img).exists(), img
    # buttons are never narrower than their text
    for b in window.findChildren(QPushButton):
        if b.text() and b.isVisible():
            assert b.width() >= b.fontMetrics().horizontalAdvance(b.text().replace("&", "")), b.text()
    page = window.single_page
    page.load_example_configuration(1); app.processEvents()
    assert "5 nodes" in page.table_summary.text() and "3 walls" in page.table_summary.text()
    page.open_table_editor(1); app.processEvents()
    assert page.table_dialog.isVisible() and page.raw_tabs.currentIndex() == 1
    assert page.walls.window() is page.table_dialog
    page.add_wall_row(); app.processEvents()
    assert "4 walls" in page.table_summary.text()
    page.table_dialog.close()
    dash._open("Single RVE Designer", 5)
    assert window.current_page() is page and page.tabs.currentIndex() == 5
    # opening the designer again after the wall-buckling tile shows the design tab, not wall buckling
    window.go_to("Dashboard"); dash._open("Single RVE Designer")
    assert window.current_page() is page and page.tabs.currentIndex() == 0
    # the dashboard has no scroll bar
    from PySide6.QtWidgets import QScrollArea
    assert not dash.findChildren(QScrollArea) and not isinstance(window.stack.widget(0), QScrollArea)
    window.close()


def test_benchmark_designs_page_and_db_wall_model_options_when_qt_is_available():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow, ValidationPage, GenerateDBPage
    app = QApplication.instance() or QApplication([])
    window = MainWindow(); window.go_to("Benchmark Designs"); app.processEvents()
    page = window.current_page(); assert isinstance(page, ValidationPage)
    headers = [page.table.horizontalHeaderItem(i).text() for i in range(page.table.columnCount())]
    assert not any("FE" in h or "Published" in h for h in headers)
    page.wall_model.beam.setCurrentIndex(1); page.wall_model.bending.setCurrentIndex(1); page.run(); app.processEvents()
    assert page.last_benchmark["wall_shear"] == "timoshenko" and page.last_benchmark["bending"] == "twist-corrected"
    assert page.table.rowCount() == 84
    window.go_to("Generate DB"); db = window.current_page(); assert isinstance(db, GenerateDBPage)
    db.wall_model.bending.setCurrentIndex(1); db.wall_model.coeff[0].setValue(3.0)
    assert db.wall_model.twist_release() == (3.0, 2.15, 0.73)
    assert window.page_names[9] == "FE Design & Validation"
    window.close()


def test_help_page_only_opens_the_two_guides_when_qt_is_available(monkeypatch):
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    import rvestudio.gui.main_window as mw
    app = QApplication.instance() or QApplication([])
    opened = []

    class FakeDesktop:                       # never start a PDF viewer on the test machine
        @staticmethod
        def openUrl(url):
            opened.append(url.toLocalFile()); return True
    monkeypatch.setattr(mw, "QDesktopServices", FakeDesktop)
    window = mw.MainWindow()
    assert set(window.help_page.guide_buttons) == {"user", "theory"}
    window.help_page.guide_buttons["theory"].click()
    assert opened and opened[-1].endswith("RVE_Studio_Theory_Guide.pdf")
    menu = [m for m in window.menuBar().actions() if m.text() == "&Help"][0].menu()
    assert [a.text() for a in menu.actions()] == ["User Guide", "Theory Guide"]
    window.close()
