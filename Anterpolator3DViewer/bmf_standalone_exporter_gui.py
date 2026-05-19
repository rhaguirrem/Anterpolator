#!/usr/bin/env python3
"""Standalone GUI for BMF export and browsing."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from PyQt5 import QtCore, QtWidgets

import bmf_standalone_exporter as bmf_tools


def _display_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if pd.isna(value):
            return ""
        return f"{value:g}"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def _metadata_frame(result: dict[str, Any]) -> pd.DataFrame:
    report = dict(result.get("report") or {})
    metadata = dict(result.get("metadata") or {})
    merged = {
        "reader_mode": result.get("reader_mode", "unknown"),
        "rows_loaded": result.get("rows_loaded", 0),
        "row_count": result.get("row_count", 0),
        "file": report.get("file", ""),
        "size_bytes": report.get("size_bytes", 0),
        "starts_with_tbms2": report.get("starts_with_tbms2", False),
        "header_fields": report.get("header_fields", {}),
    }
    for key, value in metadata.items():
        merged[key] = value
    rows = [{"key": key, "value": _display_value(value)} for key, value in merged.items()]
    return pd.DataFrame(rows, columns=["key", "value"])


class DataFrameTableModel(QtCore.QAbstractTableModel):
    def __init__(self, frame: pd.DataFrame | None = None):
        super().__init__()
        self._frame = frame if frame is not None else pd.DataFrame()

    def set_frame(self, frame: pd.DataFrame | None) -> None:
        self.beginResetModel()
        self._frame = frame if frame is not None else pd.DataFrame()
        self.endResetModel()

    def rowCount(self, parent: QtCore.QModelIndex | None = None) -> int:
        if parent and parent.isValid():
            return 0
        return int(len(self._frame.index))

    def columnCount(self, parent: QtCore.QModelIndex | None = None) -> int:
        if parent and parent.isValid():
            return 0
        return int(len(self._frame.columns))

    def data(self, index: QtCore.QModelIndex, role: int = QtCore.Qt.DisplayRole) -> Any:
        if not index.isValid() or role not in (QtCore.Qt.DisplayRole, QtCore.Qt.ToolTipRole):
            return None
        value = self._frame.iat[index.row(), index.column()]
        return _display_value(value)

    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = QtCore.Qt.DisplayRole) -> Any:
        if role != QtCore.Qt.DisplayRole:
            return None
        if orientation == QtCore.Qt.Horizontal:
            if 0 <= section < len(self._frame.columns):
                return str(self._frame.columns[section])
            return ""
        return str(section)


class BmfLoadWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal(dict)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, path: str, row_limit: int | None):
        super().__init__()
        self.path = path
        self.row_limit = row_limit

    @QtCore.pyqtSlot()
    def run(self) -> None:
        try:
            result = bmf_tools.load_bmf_table(self.path, row_limit=self.row_limit)
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


class BmfStandaloneWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Anterpolator BMF Standalone")
        self.resize(1200, 800)

        self.metadata_model = DataFrameTableModel()
        self.fields_model = DataFrameTableModel()
        self.rows_model = DataFrameTableModel()
        self._bmf_load_thread: QtCore.QThread | None = None
        self._bmf_load_worker: BmfLoadWorker | None = None
        self._bmf_load_progress: QtWidgets.QProgressDialog | None = None

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs)

        self.export_tab = QtWidgets.QWidget()
        self.browse_tab = QtWidgets.QWidget()
        self.tabs.addTab(self.export_tab, "Export")
        self.tabs.addTab(self.browse_tab, "Browse BMF")

        self._build_export_tab()
        self._build_browse_tab()
        self.statusBar().showMessage("Ready")

    def _build_export_tab(self) -> None:
        layout = QtWidgets.QVBoxLayout(self.export_tab)

        paths_group = QtWidgets.QGroupBox("Paths")
        paths_form = QtWidgets.QFormLayout(paths_group)
        self.export_input_edit = QtWidgets.QLineEdit()
        self.export_output_edit = QtWidgets.QLineEdit()
        self.export_summary_edit = QtWidgets.QLineEdit()
        paths_form.addRow("Input CSV", self._path_row(self.export_input_edit, self._browse_input_csv))
        paths_form.addRow("Output BMF", self._path_row(self.export_output_edit, self._browse_output_bmf))
        paths_form.addRow("Summary JSON", self._path_row(self.export_summary_edit, self._browse_summary_json))
        layout.addWidget(paths_group)

        options_group = QtWidgets.QGroupBox("Export Options")
        options_form = QtWidgets.QFormLayout(options_group)
        self.backend_combo = QtWidgets.QComboBox()
        self.backend_combo.addItems(["tbms-config-text", "tbms-experimental", "vulcan"])
        self.backend_combo.setCurrentText("tbms-config-text")
        self.x_col_edit = QtWidgets.QLineEdit("x")
        self.y_col_edit = QtWidgets.QLineEdit("y")
        self.z_col_edit = QtWidgets.QLineEdit("z")
        self.cell_size_edit = QtWidgets.QLineEdit()
        self.cell_size_edit.setPlaceholderText("dx,dy,dz")
        self.origin_edit = QtWidgets.QLineEdit()
        self.origin_edit.setPlaceholderText("ox,oy,oz")
        self.value_cols_edit = QtWidgets.QLineEdit()
        self.value_cols_edit.setPlaceholderText("grade,domain")
        self.null_float_edit = QtWidgets.QLineEdit("-99")
        self.index_tolerance_edit = QtWidgets.QLineEdit("1e-3")
        self.dry_run_check = QtWidgets.QCheckBox("Validate only (do not write file)")

        options_form.addRow("Backend", self.backend_combo)
        options_form.addRow("X column", self.x_col_edit)
        options_form.addRow("Y column", self.y_col_edit)
        options_form.addRow("Z column", self.z_col_edit)
        options_form.addRow("Cell size", self.cell_size_edit)
        options_form.addRow("Origin", self.origin_edit)
        options_form.addRow("Value columns", self.value_cols_edit)
        options_form.addRow("Null float", self.null_float_edit)
        options_form.addRow("Index tolerance", self.index_tolerance_edit)
        options_form.addRow("", self.dry_run_check)
        layout.addWidget(options_group)

        action_row = QtWidgets.QHBoxLayout()
        action_row.addStretch(1)
        self.export_button = QtWidgets.QPushButton("Run Export")
        self.export_button.clicked.connect(self._run_export)
        action_row.addWidget(self.export_button)
        layout.addLayout(action_row)

        self.export_result_text = QtWidgets.QPlainTextEdit()
        self.export_result_text.setReadOnly(True)
        layout.addWidget(self.export_result_text, stretch=1)

    def _build_browse_tab(self) -> None:
        layout = QtWidgets.QVBoxLayout(self.browse_tab)

        controls_group = QtWidgets.QGroupBox("BMF Reader")
        controls_form = QtWidgets.QFormLayout(controls_group)
        self.browse_input_edit = QtWidgets.QLineEdit()
        controls_form.addRow("BMF file", self._path_row(self.browse_input_edit, self._browse_input_bmf))

        row_limit_row = QtWidgets.QHBoxLayout()
        self.read_all_rows_check = QtWidgets.QCheckBox("Read all rows")
        self.read_all_rows_check.setChecked(True)
        self.read_all_rows_check.toggled.connect(self._toggle_row_limit_enabled)
        self.row_limit_spin = QtWidgets.QSpinBox()
        self.row_limit_spin.setRange(1, 100000000)
        self.row_limit_spin.setValue(1000)
        self.row_limit_spin.setEnabled(False)
        row_limit_row.addWidget(self.read_all_rows_check)
        row_limit_row.addWidget(QtWidgets.QLabel("Rows"))
        row_limit_row.addWidget(self.row_limit_spin)
        row_limit_row.addStretch(1)
        controls_form.addRow("Row limit", self._layout_widget(row_limit_row))

        button_row = QtWidgets.QHBoxLayout()
        self.load_bmf_button = QtWidgets.QPushButton("Open BMF")
        self.load_bmf_button.clicked.connect(self._load_bmf)
        self.reader_mode_label = QtWidgets.QLabel("Reader mode: n/a")
        button_row.addWidget(self.load_bmf_button)
        button_row.addWidget(self.reader_mode_label)
        button_row.addStretch(1)
        controls_form.addRow("", self._layout_widget(button_row))
        layout.addWidget(controls_group)

        self.browser_tabs = QtWidgets.QTabWidget()
        self.metadata_view = self._make_table_view(self.metadata_model)
        self.fields_view = self._make_table_view(self.fields_model)
        self.rows_view = self._make_table_view(self.rows_model)
        self.inspect_text = QtWidgets.QPlainTextEdit()
        self.inspect_text.setReadOnly(True)
        self.browser_tabs.addTab(self.metadata_view, "Metadata")
        self.browser_tabs.addTab(self.fields_view, "Fields")
        self.browser_tabs.addTab(self.rows_view, "Rows")
        self.browser_tabs.addTab(self.inspect_text, "Inspect")
        layout.addWidget(self.browser_tabs, stretch=1)

    def _make_table_view(self, model: DataFrameTableModel) -> QtWidgets.QTableView:
        view = QtWidgets.QTableView()
        view.setModel(model)
        view.setAlternatingRowColors(True)
        view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        view.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        view.horizontalHeader().setStretchLastSection(True)
        view.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        return view

    def _layout_widget(self, layout: QtWidgets.QLayout) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _path_row(self, line_edit: QtWidgets.QLineEdit, callback) -> QtWidgets.QWidget:
        row = QtWidgets.QHBoxLayout()
        row.addWidget(line_edit)
        button = QtWidgets.QPushButton("Browse...")
        button.clicked.connect(callback)
        row.addWidget(button)
        return self._layout_widget(row)

    def _browse_input_csv(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select grid CSV", "", "CSV Files (*.csv);;All Files (*)")
        if path:
            self.export_input_edit.setText(path)
            if not self.export_output_edit.text().strip():
                self.export_output_edit.setText(str(Path(path).with_suffix(".bmf")))

    def _browse_output_bmf(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Select output BMF", self.export_output_edit.text().strip(), "BMF Files (*.bmf);;All Files (*)")
        if path:
            self.export_output_edit.setText(path)

    def _browse_summary_json(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Select summary JSON", self.export_summary_edit.text().strip(), "JSON Files (*.json);;All Files (*)")
        if path:
            self.export_summary_edit.setText(path)

    def _browse_input_bmf(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open BMF file", self.browse_input_edit.text().strip(), "BMF Files (*.bmf);;All Files (*)")
        if path:
            self.browse_input_edit.setText(path)

    def _toggle_row_limit_enabled(self, checked: bool) -> None:
        self.row_limit_spin.setEnabled(not checked)

    def _parse_triplet(self, text: str) -> list[float] | None:
        stripped = text.strip()
        if not stripped:
            return None
        parts = [part.strip() for part in stripped.split(",") if part.strip()]
        if len(parts) != 3:
            raise ValueError("Expected exactly three comma-separated numeric values.")
        return [float(part) for part in parts]

    def _parse_value_cols(self, text: str) -> list[str] | None:
        parts = [part.strip() for part in text.split(",") if part.strip()]
        return parts or None

    def _run_export(self) -> None:
        try:
            result = bmf_tools.export_bmf(
                input_csv=self.export_input_edit.text().strip(),
                output_bmf=self.export_output_edit.text().strip(),
                backend=self.backend_combo.currentText(),
                x_col=self.x_col_edit.text().strip() or "x",
                y_col=self.y_col_edit.text().strip() or "y",
                z_col=self.z_col_edit.text().strip() or "z",
                cell_size=self._parse_triplet(self.cell_size_edit.text()),
                origin=self._parse_triplet(self.origin_edit.text()),
                value_cols=self._parse_value_cols(self.value_cols_edit.text()),
                null_float=float(self.null_float_edit.text().strip()),
                index_tolerance=float(self.index_tolerance_edit.text().strip()),
                dry_run=self.dry_run_check.isChecked(),
                summary_json=self.export_summary_edit.text().strip() or None,
            )
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Export Failed", str(exc))
            self.statusBar().showMessage("Export failed")
            return

        self.export_result_text.setPlainText(json.dumps(result, indent=2, default=str))
        self.statusBar().showMessage("Export completed")
        if not self.dry_run_check.isChecked() and self.export_output_edit.text().strip():
            self.browse_input_edit.setText(self.export_output_edit.text().strip())

    def _load_bmf(self) -> None:
        path = self.browse_input_edit.text().strip()
        if not path:
            QtWidgets.QMessageBox.warning(self, "Missing File", "Select a BMF file first.")
            return
        if self._bmf_load_thread is not None:
            self.statusBar().showMessage("BMF load already in progress")
            return

        row_limit = None if self.read_all_rows_check.isChecked() else int(self.row_limit_spin.value())

        progress = QtWidgets.QProgressDialog("Opening BMF file...", None, 0, 0, self)
        progress.setWindowTitle("Open BMF")
        progress.setWindowModality(QtCore.Qt.WindowModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.show()

        thread = QtCore.QThread(self)
        worker = BmfLoadWorker(path, row_limit)
        worker.moveToThread(thread)

        self._bmf_load_thread = thread
        self._bmf_load_worker = worker
        self._bmf_load_progress = progress
        self.load_bmf_button.setEnabled(False)
        self.reader_mode_label.setText("Reader mode: loading...")
        self.statusBar().showMessage("Opening BMF file...")

        def finish_cleanup() -> None:
            if self._bmf_load_progress is progress:
                progress.close()
                self._bmf_load_progress = None
            self.load_bmf_button.setEnabled(True)
            if self._bmf_load_thread is thread:
                self._bmf_load_thread = None
            if self._bmf_load_worker is worker:
                self._bmf_load_worker = None

        def handle_finished(result: dict) -> None:
            self._apply_bmf_result(result)
            thread.quit()

        def handle_failed(message: str) -> None:
            QtWidgets.QMessageBox.critical(self, "Load Failed", message)
            self.reader_mode_label.setText("Reader mode: n/a")
            self.statusBar().showMessage("BMF load failed")
            thread.quit()

        worker.finished.connect(handle_finished, QtCore.Qt.QueuedConnection)
        worker.failed.connect(handle_failed, QtCore.Qt.QueuedConnection)
        thread.started.connect(worker.run)
        worker.finished.connect(worker.deleteLater)
        worker.failed.connect(worker.deleteLater)
        thread.finished.connect(finish_cleanup)
        thread.finished.connect(thread.deleteLater)
        thread.start()

    def _apply_bmf_result(self, result: dict[str, Any]) -> None:

        self.metadata_model.set_frame(_metadata_frame(result))
        self.fields_model.set_frame(pd.DataFrame(result.get("fields") or []))
        rows_frame = result.get("dataframe")
        if not isinstance(rows_frame, pd.DataFrame):
            rows_frame = pd.DataFrame()
        self.rows_model.set_frame(rows_frame)
        self.inspect_text.setPlainText(json.dumps(result.get("report") or {}, indent=2, default=str))
        self.reader_mode_label.setText(
            f"Reader mode: {result.get('reader_mode', 'unknown')} | loaded {result.get('rows_loaded', 0)} / {result.get('row_count', 0)} rows"
        )

        if result.get("reader_mode") == "basic-inspect":
            self.statusBar().showMessage("Loaded BMF in inspection-only mode")
            QtWidgets.QMessageBox.information(
                self,
                "Inspection-only Mode",
                "This BMF does not expose the experimental section directory, so only metadata/header inspection is available.",
            )
        elif result.get("reader_mode") == "tbms-config":
            self.statusBar().showMessage("Loaded BMF config metadata")
        elif result.get("reader_mode") == "tbms-binary-pages":
            self.statusBar().showMessage("Recognized TBMS binary-page variant")
        elif result.get("reader_mode") == "tbms-binary-flat-roots":
            self.statusBar().showMessage("Loaded TBMS flat-root binary pages")
        elif result.get("reader_mode") == "tbms-string-pages":
            self.statusBar().showMessage("Loaded TBMS text-assignment pages")
        else:
            self.statusBar().showMessage("BMF loaded")


def main(argv: list[str] | None = None) -> int:
    app = QtWidgets.QApplication(argv or sys.argv)
    window = BmfStandaloneWindow()

    if argv is None:
        args = sys.argv[1:]
    else:
        args = argv[1:]

    if args:
        candidate = Path(args[0])
        if candidate.exists():
            window.browse_input_edit.setText(str(candidate))

    window.show()
    return int(app.exec_())


if __name__ == "__main__":
    raise SystemExit(main())