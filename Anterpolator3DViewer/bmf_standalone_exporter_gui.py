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


def _coerce_numeric_series(series: pd.Series, delimiter: str | None = None) -> pd.Series:
    text_series = series.astype(str).str.strip()
    if delimiter == ";":
        text_series = text_series.str.replace(",", ".", regex=False)
    return pd.to_numeric(text_series, errors="coerce")


def _normalize_bmf_field_type(field_type: object) -> str:
    normalized = str(field_type or "").strip().lower()
    aliases = {
        "": "",
        "auto": "",
        "boolean": "boolean",
        "bool": "boolean",
        "int": "int",
        "integer": "int",
        "double": "double",
        "float": "double",
        "real": "double",
        "string": "string",
        "text": "string",
        "namedshort": "string",
    }
    if normalized not in aliases:
        raise ValueError("Unsupported BMF export field type. Use one of: auto, boolean, int, double, string.")
    return aliases[normalized]


def _infer_bmf_field_types_from_preview(frame: pd.DataFrame, candidate_columns: list[str], delimiter: str | None = None) -> dict[str, str]:
    inferred = {}
    truthy = {"1", "true", "t", "yes", "y"}
    falsy = {"0", "false", "f", "no", "n", ""}
    for column_name in candidate_columns:
        if column_name not in frame.columns:
            continue
        series = frame[column_name]
        non_null = series.dropna()
        if len(non_null) == 0:
            inferred[column_name] = ""
            continue
        normalized_text = non_null.astype(str).str.strip().str.lower()
        normalized_text = normalized_text[normalized_text != ""]
        if len(normalized_text) > 0 and normalized_text.isin(truthy | falsy).all():
            inferred[column_name] = "boolean"
            continue
        numeric_series = _coerce_numeric_series(series, delimiter=delimiter)
        blank_mask = series.isna() | (series.astype(str).str.strip() == "")
        if bool((numeric_series.notna() | blank_mask).all()):
            finite = numeric_series.dropna().to_numpy(dtype=float)
            inferred[column_name] = "int" if finite.size and (abs(finite - finite.round()) < 1e-9).all() else "double"
            continue
        inferred[column_name] = "string"
    return inferred


class TrimmedDisplayDoubleSpinBox(QtWidgets.QDoubleSpinBox):
    def textFromValue(self, value: float) -> str:
        text = f"{float(value):.{max(int(self.decimals()), 0)}f}"
        if "." in text:
            text = text.rstrip("0").rstrip(".")
        return "0" if text in {"-0", "-0.0", ""} else text


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


class BmfExportWorker(QtCore.QObject):
    progress = QtCore.pyqtSignal(int, int, str)
    finished = QtCore.pyqtSignal(dict)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, export_kwargs: dict[str, Any]):
        super().__init__()
        self.export_kwargs = dict(export_kwargs)

    @QtCore.pyqtSlot()
    def run(self) -> None:
        def handle_progress(value: int, maximum: int = 100, message: str = "") -> None:
            message_text = str(message or "")
            print(f"BMF export progress: {int(value)}/{int(maximum)} {message_text}")
            self.progress.emit(int(value), int(maximum), message_text)

        try:
            result = bmf_tools.export_bmf(progress_callback=handle_progress, **self.export_kwargs)
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
        self._loaded_bmf_dataframe: pd.DataFrame | None = None
        self._bmf_load_thread: QtCore.QThread | None = None
        self._bmf_load_worker: BmfLoadWorker | None = None
        self._bmf_load_progress: QtWidgets.QProgressDialog | None = None
        self._bmf_export_thread: QtCore.QThread | None = None
        self._bmf_export_worker: BmfExportWorker | None = None
        self._bmf_export_progress: QtWidgets.QProgressDialog | None = None
        self._pending_export_value_cols: list[str] | None = None
        self._pending_export_column_types: dict[str, str] | None = None
        self._export_value_cols_initialized = False
        self._export_columns_source_path = ""
        self._metadata_cell_size_values: list[float] | None = None

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs)

        self.export_tab = QtWidgets.QWidget()
        self.browse_tab = QtWidgets.QWidget()
        self.tabs.addTab(self.export_tab, "csv2bmf")
        self.tabs.addTab(self.browse_tab, "bmf2csv")

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

        csv2bmf_tabs = QtWidgets.QTabWidget()
        csv_format_tab = QtWidgets.QWidget()
        export_options_tab = QtWidgets.QWidget()
        csv_format_form = QtWidgets.QFormLayout(csv_format_tab)
        export_options_form = QtWidgets.QFormLayout(export_options_tab)
        csv2bmf_tabs.addTab(csv_format_tab, "CSV Format")
        csv2bmf_tabs.addTab(export_options_tab, "Export Options")
        self.delimiter_combo = QtWidgets.QComboBox()
        self.delimiter_combo.addItems([",", ";", "\t", "|"])
        self.header_line_spin = QtWidgets.QSpinBox()
        self.header_line_spin.setRange(1, 1_000_000)
        self.header_line_spin.setValue(1)
        self.backend_combo = QtWidgets.QComboBox()
        self.backend_combo.addItems(["tbms-config-text", "tbms-experimental", "vulcan"])
        self.backend_combo.setCurrentText("tbms-config-text")
        self.cell_x_spin = self._make_cell_size_spin()
        self.cell_y_spin = self._make_cell_size_spin()
        self.cell_z_spin = self._make_cell_size_spin()
        self.regularize_base_blocks_check = QtWidgets.QCheckBox("Regularize to base block size")
        self.regularize_base_blocks_check.setChecked(False)
        self.regularize_warning_label = QtWidgets.QLabel(
            "Info: unchecked preserves the source CSV block rows; checked aggregates rows to base blocks and may create a dense regular grid."
        )
        self.regularize_warning_label.setWordWrap(True)
        self.regularize_warning_label.setStyleSheet("color: #b36b00;")
        self.x_col_combo = QtWidgets.QComboBox()
        self.y_col_combo = QtWidgets.QComboBox()
        self.z_col_combo = QtWidgets.QComboBox()
        self.geometry_mode_combo = QtWidgets.QComboBox()
        self.geometry_mode_combo.addItem("Infer from metadata/centroids", "infer")
        self.geometry_mode_combo.addItem("Use size columns dX/dY/dZ", "size")
        self.geometry_mode_combo.addItem("Use lower/upper extent columns", "extent")
        self.size_x_col_combo = QtWidgets.QComboBox()
        self.size_y_col_combo = QtWidgets.QComboBox()
        self.size_z_col_combo = QtWidgets.QComboBox()
        self.lower_x_col_combo = QtWidgets.QComboBox()
        self.upper_x_col_combo = QtWidgets.QComboBox()
        self.lower_y_col_combo = QtWidgets.QComboBox()
        self.upper_y_col_combo = QtWidgets.QComboBox()
        self.lower_z_col_combo = QtWidgets.QComboBox()
        self.upper_z_col_combo = QtWidgets.QComboBox()
        for combo in [
            self.x_col_combo,
            self.y_col_combo,
            self.z_col_combo,
            self.size_x_col_combo,
            self.size_y_col_combo,
            self.size_z_col_combo,
            self.lower_x_col_combo,
            self.upper_x_col_combo,
            self.lower_y_col_combo,
            self.upper_y_col_combo,
            self.lower_z_col_combo,
            self.upper_z_col_combo,
        ]:
            self._configure_column_combo(combo)
        self.value_cols_table = QtWidgets.QTableWidget(0, 2)
        self.value_cols_table.setHorizontalHeaderLabels(["Field", "Type"])
        self.value_cols_table.verticalHeader().setVisible(False)
        self.value_cols_table.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        self.value_cols_table.setMinimumHeight(160)
        self.value_cols_table.horizontalHeader().setStretchLastSection(False)
        self.value_cols_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.value_cols_table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
        self.select_all_value_cols_button = QtWidgets.QPushButton("Select All")
        self.clear_value_cols_button = QtWidgets.QPushButton("Clear")
        self.value_summary_label = QtWidgets.QLabel("No export value columns loaded.")
        self.value_summary_label.setWordWrap(True)
        self.value_exceptions = {}
        self.exception_table = QtWidgets.QTableWidget(0, 4)
        self.exception_table.setHorizontalHeaderLabels(["Field", "CSV Value", "Replacement", "Use In Regularization"])
        self.exception_table.verticalHeader().setVisible(False)
        self.exception_table.setMinimumHeight(120)
        self.exception_table.horizontalHeader().setStretchLastSection(False)
        self.exception_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        self.exception_table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        self.exception_table.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
        self.exception_table.horizontalHeader().setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeToContents)
        self.add_exception_button = QtWidgets.QPushButton("Add Exception")
        self.remove_exception_button = QtWidgets.QPushButton("Remove Selected")
        self.exception_summary_label = QtWidgets.QLabel("No BMF value exceptions configured.")
        self.exception_summary_label.setWordWrap(True)
        self.origin_edit = QtWidgets.QLineEdit()
        self.origin_edit.setPlaceholderText("ox,oy,oz")
        self.null_float_edit = QtWidgets.QLineEdit("-99")
        self.index_tolerance_edit = QtWidgets.QLineEdit("1e-3")
        self.dry_run_check = QtWidgets.QCheckBox("Validate only (do not write file)")

        csv_format_form.addRow("Input Delimiter", self.delimiter_combo)
        csv_format_form.addRow("Input Header Line", self.header_line_spin)
        export_options_form.addRow("Backend", self.backend_combo)
        cell_layout = QtWidgets.QHBoxLayout()
        cell_layout.addWidget(self.cell_x_spin)
        cell_layout.addWidget(self.cell_y_spin)
        cell_layout.addWidget(self.cell_z_spin)
        csv_format_form.addRow("Cell Size X / Y / Z", self._layout_widget(cell_layout))
        regularize_layout = QtWidgets.QVBoxLayout()
        regularize_layout.addWidget(self.regularize_base_blocks_check)
        regularize_layout.addWidget(self.regularize_warning_label)
        export_options_form.addRow("", self._layout_widget(regularize_layout))
        coords_layout = QtWidgets.QHBoxLayout()
        coords_layout.addWidget(self.x_col_combo)
        coords_layout.addWidget(self.y_col_combo)
        coords_layout.addWidget(self.z_col_combo)
        csv_format_form.addRow("X / Y / Z Columns", self._layout_widget(coords_layout))
        csv_format_form.addRow("Block Size Source", self.geometry_mode_combo)
        size_cols_layout = QtWidgets.QHBoxLayout()
        size_cols_layout.addWidget(self.size_x_col_combo)
        size_cols_layout.addWidget(self.size_y_col_combo)
        size_cols_layout.addWidget(self.size_z_col_combo)
        csv_format_form.addRow("Size dX / dY / dZ Columns", self._layout_widget(size_cols_layout))
        extent_cols_layout = QtWidgets.QGridLayout()
        extent_cols_layout.addWidget(QtWidgets.QLabel("Lower X"), 0, 0)
        extent_cols_layout.addWidget(QtWidgets.QLabel("Upper X"), 0, 1)
        extent_cols_layout.addWidget(QtWidgets.QLabel("Lower Y"), 0, 2)
        extent_cols_layout.addWidget(QtWidgets.QLabel("Upper Y"), 0, 3)
        extent_cols_layout.addWidget(QtWidgets.QLabel("Lower Z"), 0, 4)
        extent_cols_layout.addWidget(QtWidgets.QLabel("Upper Z"), 0, 5)
        extent_cols_layout.addWidget(self.lower_x_col_combo, 1, 0)
        extent_cols_layout.addWidget(self.upper_x_col_combo, 1, 1)
        extent_cols_layout.addWidget(self.lower_y_col_combo, 1, 2)
        extent_cols_layout.addWidget(self.upper_y_col_combo, 1, 3)
        extent_cols_layout.addWidget(self.lower_z_col_combo, 1, 4)
        extent_cols_layout.addWidget(self.upper_z_col_combo, 1, 5)
        csv_format_form.addRow("Lower/Upper Extent Columns", self._layout_widget(extent_cols_layout))
        export_options_form.addRow("Value Columns", self.value_cols_table)
        value_buttons_layout = QtWidgets.QHBoxLayout()
        value_buttons_layout.addWidget(self.select_all_value_cols_button)
        value_buttons_layout.addWidget(self.clear_value_cols_button)
        export_options_form.addRow("", self._layout_widget(value_buttons_layout))
        export_options_form.addRow("", self.value_summary_label)
        export_options_form.addRow("Value Exceptions", self.exception_table)
        exception_buttons_layout = QtWidgets.QHBoxLayout()
        exception_buttons_layout.addWidget(self.add_exception_button)
        exception_buttons_layout.addWidget(self.remove_exception_button)
        export_options_form.addRow("", self._layout_widget(exception_buttons_layout))
        export_options_form.addRow("", self.exception_summary_label)
        export_options_form.addRow("Origin", self.origin_edit)
        export_options_form.addRow("Null float", self.null_float_edit)
        export_options_form.addRow("Index tolerance", self.index_tolerance_edit)
        export_options_form.addRow("", self.dry_run_check)
        layout.addWidget(csv2bmf_tabs)

        self.export_input_edit.textChanged.connect(lambda _: self._refresh_export_columns())
        self.export_input_edit.textChanged.connect(lambda _: self._sync_export_output_path())
        self.export_input_edit.textChanged.connect(lambda _: self._refresh_metadata_cell_size())
        self.delimiter_combo.currentIndexChanged.connect(lambda _: self._refresh_export_columns())
        self.backend_combo.currentIndexChanged.connect(lambda _: self._refresh_metadata_cell_size())
        self.backend_combo.currentIndexChanged.connect(lambda _: self._update_regularize_warning())
        self.header_line_spin.valueChanged.connect(lambda _: self._refresh_export_columns())
        self.x_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.y_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.z_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.geometry_mode_combo.currentIndexChanged.connect(lambda _: self._update_geometry_mode_controls())
        self.geometry_mode_combo.currentIndexChanged.connect(lambda _: self._refresh_export_columns())
        self.size_x_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.size_y_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.size_z_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.lower_x_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.upper_x_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.lower_y_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.upper_y_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.lower_z_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.upper_z_col_combo.currentTextChanged.connect(lambda _: self._refresh_export_columns())
        self.value_cols_table.itemChanged.connect(lambda _item: self._update_value_summary())
        self.select_all_value_cols_button.clicked.connect(lambda: self._set_all_value_columns_checked(True))
        self.clear_value_cols_button.clicked.connect(lambda: self._set_all_value_columns_checked(False))
        self.add_exception_button.clicked.connect(self._add_exception_row)
        self.remove_exception_button.clicked.connect(self._remove_selected_exception_rows)
        self.exception_table.itemChanged.connect(lambda _item: self._update_exception_summary())
        self.regularize_base_blocks_check.toggled.connect(self._update_regularize_warning)
        self.regularize_base_blocks_check.toggled.connect(lambda _: self._refresh_metadata_cell_size())
        self._update_regularize_warning()
        self._update_geometry_mode_controls()

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
        self.browse_output_edit = QtWidgets.QLineEdit()
        controls_form.addRow("BMF file", self._path_row(self.browse_input_edit, self._browse_input_bmf))
        controls_form.addRow("Output CSV", self._path_row(self.browse_output_edit, self._browse_output_csv))

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
        self.load_bmf_button = QtWidgets.QPushButton("Preview BMF")
        self.load_bmf_button.clicked.connect(self._load_bmf)
        self.save_bmf_csv_button = QtWidgets.QPushButton("Save CSV")
        self.save_bmf_csv_button.clicked.connect(self._save_bmf_csv)
        self.save_bmf_csv_button.setEnabled(False)
        self.reader_mode_label = QtWidgets.QLabel("Reader mode: n/a")
        button_row.addWidget(self.load_bmf_button)
        button_row.addWidget(self.save_bmf_csv_button)
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

    def _make_cell_size_spin(self) -> QtWidgets.QDoubleSpinBox:
        spin = TrimmedDisplayDoubleSpinBox()
        spin.setRange(0.0, 1_000_000_000.0)
        spin.setDecimals(6)
        spin.setSingleStep(1.0)
        spin.setValue(0.0)
        spin.setSpecialValueText("Auto")
        return spin

    def _configure_column_combo(self, combo: QtWidgets.QComboBox) -> None:
        combo.setEditable(True)
        combo.setMinimumWidth(170)
        combo.view().setMinimumWidth(320)
        combo.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)

    def _geometry_mode(self) -> str:
        return str(self.geometry_mode_combo.currentData() or "infer")

    def _update_geometry_mode_controls(self) -> None:
        mode = self._geometry_mode()
        use_size_columns = mode == "size"
        use_extent_columns = mode == "extent"
        for combo in (self.size_x_col_combo, self.size_y_col_combo, self.size_z_col_combo):
            combo.setEnabled(use_size_columns)
        for combo in self._extent_column_combos():
            combo.setEnabled(use_extent_columns)

    def _extent_column_combos(self) -> tuple[QtWidgets.QComboBox, ...]:
        return (
            self.lower_x_col_combo,
            self.upper_x_col_combo,
            self.lower_y_col_combo,
            self.upper_y_col_combo,
            self.lower_z_col_combo,
            self.upper_z_col_combo,
        )

    def _sync_export_output_path(self) -> None:
        input_path = self.export_input_edit.text().strip()
        if input_path and not self.export_output_edit.text().strip():
            self.export_output_edit.setText(str(Path(input_path).with_suffix(".bmf")))

    def _get_cell_size(self) -> list[float] | None:
        values = [float(self.cell_x_spin.value()), float(self.cell_y_spin.value()), float(self.cell_z_spin.value())]
        positive_values = [value for value in values if value > 0]
        if not positive_values:
            return None
        if len(positive_values) != 3:
            raise ValueError("BMF cell size must be set for X, Y, and Z, or left as Auto for all three axes.")
        return values

    def _cell_size_can_accept_metadata(self) -> bool:
        current_values = [float(self.cell_x_spin.value()), float(self.cell_y_spin.value()), float(self.cell_z_spin.value())]
        if all(value <= 0 for value in current_values):
            return True
        previous = self._metadata_cell_size_values
        if previous is None:
            return False
        return all(abs(current - old) <= 1e-9 for current, old in zip(current_values, previous))

    def _set_cell_size_from_metadata(self, cell_size: list[float], source: str = "") -> None:
        if len(cell_size) != 3 or any(float(value) <= 0 for value in cell_size):
            return
        if not self._cell_size_can_accept_metadata():
            return

        blockers = [
            QtCore.QSignalBlocker(self.cell_x_spin),
            QtCore.QSignalBlocker(self.cell_y_spin),
            QtCore.QSignalBlocker(self.cell_z_spin),
        ]
        try:
            for spin, value in zip((self.cell_x_spin, self.cell_y_spin, self.cell_z_spin), cell_size):
                spin.setValue(float(value))
        finally:
            del blockers
        self._metadata_cell_size_values = [float(value) for value in cell_size]
        source_text = f" from {source}" if source else ""
        self.statusBar().showMessage(f"Filled Cell Size X/Y/Z from CSV metadata{source_text}.")

    def _refresh_metadata_cell_size(self) -> None:
        path = self.export_input_edit.text().strip()
        if not path or not Path(path).is_file():
            self._metadata_cell_size_values = None
            self._update_regularize_warning()
            return
        try:
            metadata = bmf_tools.parse_leapfrog_block_metadata(path)
            geometry = bmf_tools.infer_bmf_export_geometry_from_metadata(
                metadata,
                regularize_to_base_block=self.regularize_base_blocks_check.isChecked(),
                dense_regular_grid=self.backend_combo.currentText() in {"tbms-config-text", "vulcan"},
            )
        except Exception:
            return
        cell_size = geometry.get("cell_size") if isinstance(geometry, dict) else None
        if isinstance(cell_size, (list, tuple)):
            self._set_cell_size_from_metadata([float(value) for value in cell_size], str(geometry.get("cell_size_source") or "metadata"))
        self._update_regularize_warning()

    def _get_selected_value_columns(self) -> list[str]:
        selected = []
        for row_index in range(self.value_cols_table.rowCount()):
            item = self.value_cols_table.item(row_index, 0)
            if item is not None and item.checkState() == QtCore.Qt.Checked:
                selected.append(item.text())
        return selected

    def _set_all_value_columns_checked(self, checked: bool) -> None:
        blocker = QtCore.QSignalBlocker(self.value_cols_table)
        try:
            for row_index in range(self.value_cols_table.rowCount()):
                item = self.value_cols_table.item(row_index, 0)
                if item is not None:
                    item.setCheckState(QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
        finally:
            del blocker
        self._update_value_summary()

    def _update_value_summary(self) -> None:
        total = self.value_cols_table.rowCount()
        selected = len(self._get_selected_value_columns())
        if total == 0:
            self.value_summary_label.setText("No export value columns loaded.")
        elif selected == 0:
            self.value_summary_label.setText(f"0 of {total} value columns selected. Select one or more columns to export.")
        else:
            self.value_summary_label.setText(f"{selected} of {total} eligible non-coordinate columns selected for export.")

    def _get_column_type_overrides(self, include_unselected: bool = False) -> dict[str, str]:
        overrides = {}
        for row_index in range(self.value_cols_table.rowCount()):
            item = self.value_cols_table.item(row_index, 0)
            combo = self.value_cols_table.cellWidget(row_index, 1)
            if item is None or combo is None:
                continue
            if not include_unselected and item.checkState() != QtCore.Qt.Checked:
                continue
            field_type = _normalize_bmf_field_type(combo.currentData() or combo.currentText())
            inferred_type = _normalize_bmf_field_type(combo.property("bmf_inferred_type") or "")
            is_explicit = bool(combo.property("bmf_type_is_explicit"))
            if field_type and (field_type != inferred_type or is_explicit):
                overrides[item.text()] = field_type
        return overrides

    def _make_exception_include_item(self, include_in_regularization: bool = False) -> QtWidgets.QTableWidgetItem:
        item = QtWidgets.QTableWidgetItem("")
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(QtCore.Qt.Checked if include_in_regularization else QtCore.Qt.Unchecked)
        item.setToolTip("When checked, the replacement is applied before base-block regularization and can contribute to numeric averages.")
        return item

    def _add_exception_row(self, column_name: str = "", bad_value: str = "", replacement: str = "", include_in_regularization: bool = False) -> None:
        if not column_name:
            selected_columns = self._get_selected_value_columns()
            column_name = selected_columns[0] if selected_columns else ""
        row_index = self.exception_table.rowCount()
        self.exception_table.insertRow(row_index)
        self.exception_table.setItem(row_index, 0, QtWidgets.QTableWidgetItem(str(column_name or "")))
        self.exception_table.setItem(row_index, 1, QtWidgets.QTableWidgetItem(str(bad_value or "")))
        self.exception_table.setItem(row_index, 2, QtWidgets.QTableWidgetItem("" if replacement is None else str(replacement)))
        self.exception_table.setItem(row_index, 3, self._make_exception_include_item(include_in_regularization))
        self.exception_table.setCurrentCell(row_index, 1)
        self._update_exception_summary()

    def _remove_selected_exception_rows(self) -> None:
        rows = sorted({index.row() for index in self.exception_table.selectedIndexes()}, reverse=True)
        if not rows and self.exception_table.currentRow() >= 0:
            rows = [self.exception_table.currentRow()]
        for row_index in rows:
            self.exception_table.removeRow(row_index)
        self._update_exception_summary()

    def _get_value_exceptions(self) -> dict[str, dict[str, object]]:
        rules: dict[str, dict[str, object]] = {}
        for row_index in range(self.exception_table.rowCount()):
            column_item = self.exception_table.item(row_index, 0)
            value_item = self.exception_table.item(row_index, 1)
            replacement_item = self.exception_table.item(row_index, 2)
            include_item = self.exception_table.item(row_index, 3)
            column_name = str(column_item.text() if column_item is not None else "").strip()
            bad_value = str(value_item.text() if value_item is not None else "")
            replacement = str(replacement_item.text() if replacement_item is not None else "")
            include_in_regularization = bool(include_item is not None and include_item.checkState() == QtCore.Qt.Checked)
            if not column_name or not bad_value:
                continue
            if include_in_regularization:
                rules.setdefault(column_name, {})[bad_value] = {
                    "replacement": replacement,
                    "include_in_regularization": True,
                }
            else:
                rules.setdefault(column_name, {})[bad_value] = replacement
        return rules

    def _update_exception_summary(self) -> None:
        self.value_exceptions = self._get_value_exceptions()
        count = sum(len(column_rules) for column_rules in self.value_exceptions.values())
        if count == 0:
            self.exception_summary_label.setText("No BMF value exceptions configured.")
        else:
            self.exception_summary_label.setText(
                f"{count} BMF value exception rule(s) configured. Checked rules are applied before regularization; blank replacements are exported as null/default values."
            )

    def _update_regularize_warning(self) -> None:
        path = self.export_input_edit.text().strip()
        has_subblock_metadata = False
        if path and Path(path).is_file():
            try:
                metadata = bmf_tools.parse_leapfrog_block_metadata(path)
                has_subblock_metadata = bool(metadata.get("subblock_factors"))
            except Exception:
                has_subblock_metadata = False
        backend = self.backend_combo.currentText()
        is_dense_backend = backend in {"tbms-config-text", "vulcan"}

        if has_subblock_metadata and is_dense_backend:
            if self.regularize_base_blocks_check.isChecked():
                self.regularize_warning_label.setText(
                    "Warning: regularization collapses source sub-block rows to parent blocks; mixed parent/sub-block hierarchy is not preserved."
                )
            elif backend == "tbms-config-text":
                self.regularize_warning_label.setText(
                    "Info: tbms-config-text will write source sub-block rows as irregular BMF row extents instead of creating a dense minimum-cell grid."
                )
            else:
                self.regularize_warning_label.setText(
                    "Warning: the Vulcan backend path currently creates a regular model; use tbms-config-text to preserve mixed parent/sub-block row extents."
                )
            self.regularize_warning_label.setVisible(True)
            return

        if self.regularize_base_blocks_check.isChecked():
            self.regularize_warning_label.setText(
                "Warning: regularization aggregates source CSV rows to base blocks and may create a dense regular grid."
            )
        elif backend == "tbms-config-text":
            self.regularize_warning_label.setText(
                "Info: tbms-config-text will preserve source CSV rows as BMF rows instead of regularizing to a dense grid."
            )
        else:
            self.regularize_warning_label.setText(
                "Info: unchecked preserves source CSV rows when the selected backend supports row-indexed or irregular output."
            )
        self.regularize_warning_label.setVisible(True)

    def _refresh_export_columns(self) -> None:
        path = self.export_input_edit.text().strip()
        delimiter = self.delimiter_combo.currentText() or ","
        header_line = int(self.header_line_spin.value())
        pending_selection = self._pending_export_value_cols
        pending_column_types = self._pending_export_column_types
        current_selection = set(self._get_selected_value_columns())
        current_column_types = self._get_column_type_overrides(include_unselected=True)
        path_changed = path != self._export_columns_source_path
        if path_changed and pending_selection is None:
            current_selection = set()
            current_column_types = {}
            self._export_value_cols_initialized = False

        current_x = self.x_col_combo.currentText()
        current_y = self.y_col_combo.currentText()
        current_z = self.z_col_combo.currentText()
        current_size_x = self.size_x_col_combo.currentText()
        current_size_y = self.size_y_col_combo.currentText()
        current_size_z = self.size_z_col_combo.currentText()
        current_extent_cols = [combo.currentText() for combo in self._extent_column_combos()]

        blockers = [
            QtCore.QSignalBlocker(self.x_col_combo),
            QtCore.QSignalBlocker(self.y_col_combo),
            QtCore.QSignalBlocker(self.z_col_combo),
            QtCore.QSignalBlocker(self.size_x_col_combo),
            QtCore.QSignalBlocker(self.size_y_col_combo),
            QtCore.QSignalBlocker(self.size_z_col_combo),
            *[QtCore.QSignalBlocker(combo) for combo in self._extent_column_combos()],
            QtCore.QSignalBlocker(self.value_cols_table),
        ]
        try:
            for combo in [self.x_col_combo, self.y_col_combo, self.z_col_combo]:
                combo.clear()
            for combo in [self.size_x_col_combo, self.size_y_col_combo, self.size_z_col_combo]:
                combo.clear()
            for combo in self._extent_column_combos():
                combo.clear()
            self.value_cols_table.setRowCount(0)

            if not path or not Path(path).is_file():
                self._export_columns_source_path = path
                self._update_value_summary()
                return

            try:
                columns = bmf_tools.parse_effective_header_line(path, delimiter, header_line)
            except Exception:
                self._export_columns_source_path = path
                self._update_value_summary()
                return

            for combo in [self.x_col_combo, self.y_col_combo, self.z_col_combo]:
                combo.addItems(columns)
            for combo in [self.size_x_col_combo, self.size_y_col_combo, self.size_z_col_combo]:
                combo.addItem("(None)")
                combo.addItems(columns)
            for combo in self._extent_column_combos():
                combo.addItem("(None)")
                combo.addItems(columns)

            def restore_or_suggest(combo: QtWidgets.QComboBox, current_value: str, keywords: list[str]) -> None:
                if current_value:
                    index = combo.findText(current_value)
                    if index >= 0:
                        combo.setCurrentIndex(index)
                        return
                for keyword in keywords:
                    for index in range(combo.count()):
                        if combo.itemText(index).lower() == keyword:
                            combo.setCurrentIndex(index)
                            return
                if combo.count() > 0 and combo.currentIndex() < 0:
                    combo.setCurrentIndex(0)

            restore_or_suggest(self.x_col_combo, current_x, ["x", "easting"])
            restore_or_suggest(self.y_col_combo, current_y, ["y", "northing"])
            restore_or_suggest(self.z_col_combo, current_z, ["z", "elevation", "rl"])
            if self._geometry_mode() == "size":
                restore_or_suggest(self.size_x_col_combo, current_size_x, ["dx", "dX", "DX", "size_x", "width_x"])
                restore_or_suggest(self.size_y_col_combo, current_size_y, ["dy", "dY", "DY", "size_y", "width_y"])
                restore_or_suggest(self.size_z_col_combo, current_size_z, ["dz", "dZ", "DZ", "size_z", "width_z"])
            else:
                for combo in (self.size_x_col_combo, self.size_y_col_combo, self.size_z_col_combo):
                    if combo.count() > 0:
                        combo.setCurrentIndex(0)
            if self._geometry_mode() == "extent":
                extent_suggestions = [
                    ["__lower_x", "lower_x", "x_lower", "min_x", "x_min"],
                    ["__upper_x", "upper_x", "x_upper", "max_x", "x_max"],
                    ["__lower_y", "lower_y", "y_lower", "min_y", "y_min"],
                    ["__upper_y", "upper_y", "y_upper", "max_y", "y_max"],
                    ["__lower_z", "lower_z", "z_lower", "min_z", "z_min"],
                    ["__upper_z", "upper_z", "z_upper", "max_z", "z_max"],
                ]
                for combo, current_value, suggestions in zip(self._extent_column_combos(), current_extent_cols, extent_suggestions):
                    restore_or_suggest(combo, current_value, suggestions)
            else:
                for combo in self._extent_column_combos():
                    if combo.count() > 0:
                        combo.setCurrentIndex(0)

            excluded = {
                str(self.x_col_combo.currentText() or "").strip(),
                str(self.y_col_combo.currentText() or "").strip(),
                str(self.z_col_combo.currentText() or "").strip(),
            }
            excluded.update(self._current_geometry_excluded_columns())
            eligible_columns = [column_name for column_name in columns if column_name not in excluded]
            desired_selection = set(pending_selection) if pending_selection is not None else current_selection
            desired_column_types = dict(pending_column_types) if pending_column_types is not None else current_column_types
            try:
                effective_header_line = bmf_tools.resolve_effective_csv_header_line(path, header_line)
                preview = pd.read_csv(
                    path,
                    sep=delimiter,
                    header=effective_header_line - 1,
                    comment="#",
                    nrows=1000,
                    low_memory=False,
                )
                inferred_column_types = _infer_bmf_field_types_from_preview(preview, eligible_columns, delimiter=delimiter)
            except Exception:
                inferred_column_types = {}
            desired_selection = {name for name in desired_selection if name in eligible_columns}
            if not desired_selection and not self._export_value_cols_initialized:
                desired_selection = set(eligible_columns)

            type_options = [("Auto", ""), ("Boolean", "boolean"), ("Integer", "int"), ("Double", "double"), ("String", "string")]
            self.value_cols_table.setRowCount(len(eligible_columns))
            for row_index, column_name in enumerate(eligible_columns):
                item = QtWidgets.QTableWidgetItem(column_name)
                item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.Checked if column_name in desired_selection else QtCore.Qt.Unchecked)
                self.value_cols_table.setItem(row_index, 0, item)

                combo = QtWidgets.QComboBox()
                for label, value in type_options:
                    combo.addItem(label, value)
                has_explicit_type = column_name in desired_column_types
                inferred_type = inferred_column_types.get(column_name, "")
                default_type = desired_column_types.get(column_name, inferred_type)
                combo.setCurrentIndex(max(combo.findData(default_type), 0))
                combo.setProperty("bmf_inferred_type", inferred_type if not has_explicit_type else "")
                combo.setProperty("bmf_type_is_explicit", has_explicit_type)
                combo.activated.connect(lambda _index, cb=combo: cb.setProperty("bmf_type_is_explicit", True))
                self.value_cols_table.setCellWidget(row_index, 1, combo)

            self._pending_export_value_cols = None
            self._pending_export_column_types = None
            self._export_value_cols_initialized = True
            self._export_columns_source_path = path
        finally:
            del blockers
        self._update_value_summary()

    def _browse_input_csv(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select grid CSV", "", "CSV Files (*.csv);;All Files (*)")
        if path:
            self.export_input_edit.setText(path)
            detected = bmf_tools.detect_csv_delimiter(path)
            if detected in [self.delimiter_combo.itemText(index) for index in range(self.delimiter_combo.count())]:
                self.delimiter_combo.setCurrentText(detected)
            if not self.export_output_edit.text().strip():
                self.export_output_edit.setText(str(Path(path).with_suffix(".bmf")))
            self._refresh_export_columns()

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
            if not self.browse_output_edit.text().strip():
                self.browse_output_edit.setText(str(Path(path).with_suffix(".csv")))

    def _browse_output_csv(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Select output CSV", self.browse_output_edit.text().strip(), "CSV Files (*.csv);;All Files (*)")
        if path:
            self.browse_output_edit.setText(path)

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

    def _combo_column_or_none(self, combo: QtWidgets.QComboBox) -> str | None:
        text = str(combo.currentText() or "").strip()
        if not text or text == "(None)":
            return None
        return text

    def _get_size_cols(self) -> list[str] | None:
        if self._geometry_mode() != "size":
            return None
        columns = [
            self._combo_column_or_none(self.size_x_col_combo),
            self._combo_column_or_none(self.size_y_col_combo),
            self._combo_column_or_none(self.size_z_col_combo),
        ]
        if all(column is None for column in columns):
            return None
        if any(column is None for column in columns):
            raise ValueError("Size columns must be selected for dX, dY, and dZ, or left blank for all three axes.")
        return [str(column) for column in columns]

    def _parse_extent_cols(self) -> list[str] | None:
        if self._geometry_mode() != "extent":
            return None
        columns = [self._combo_column_or_none(combo) for combo in self._extent_column_combos()]
        if any(column is None for column in columns):
            raise ValueError("Lower/upper extent column mode requires selecting all six extent columns.")
        return [str(column) for column in columns]

    def _parse_extent_cols_for_exclusion(self) -> list[str] | None:
        columns = [self._combo_column_or_none(combo) for combo in self._extent_column_combos()]
        if any(column is None for column in columns):
            return None
        return [str(column) for column in columns]

    def _current_geometry_excluded_columns(self) -> set[str]:
        excluded = set()
        mode = self._geometry_mode()
        if mode == "size":
            for combo in (self.size_x_col_combo, self.size_y_col_combo, self.size_z_col_combo):
                column_name = self._combo_column_or_none(combo)
                if column_name:
                    excluded.add(column_name)
        elif mode == "extent":
            extent_columns = self._parse_extent_cols_for_exclusion()
            if extent_columns:
                excluded.update(extent_columns)
        return excluded

    def _run_export(self) -> None:
        if self._bmf_export_thread is not None:
            self.statusBar().showMessage("BMF export already in progress")
            return
        try:
            value_columns = self._get_selected_value_columns()
            if not value_columns:
                raise ValueError("Please select at least one value column to export.")
            column_types = self._get_column_type_overrides(include_unselected=False)
            unknown_typed_columns = [column_name for column_name in column_types if column_name not in value_columns]
            if unknown_typed_columns:
                raise ValueError(
                    "Column type overrides must refer only to selected value columns. "
                    f"Unknown or unselected columns: {unknown_typed_columns}"
                )
            export_kwargs = {
                "input_csv": self.export_input_edit.text().strip(),
                "output_bmf": self.export_output_edit.text().strip(),
                "backend": self.backend_combo.currentText(),
                "x_col": self.x_col_combo.currentText().strip() or "x",
                "y_col": self.y_col_combo.currentText().strip() or "y",
                "z_col": self.z_col_combo.currentText().strip() or "z",
                "delimiter": self.delimiter_combo.currentText(),
                "header_line": int(self.header_line_spin.value()),
                "cell_size": self._get_cell_size(),
                "origin": self._parse_triplet(self.origin_edit.text()),
                "size_cols": self._get_size_cols(),
                "extent_cols": self._parse_extent_cols(),
                "value_cols": value_columns,
                "column_types": column_types,
                "value_exceptions": self._get_value_exceptions(),
                "null_float": float(self.null_float_edit.text().strip()),
                "index_tolerance": float(self.index_tolerance_edit.text().strip()),
                "regularize_to_base_block": self.regularize_base_blocks_check.isChecked(),
                "dry_run": self.dry_run_check.isChecked(),
                "summary_json": self.export_summary_edit.text().strip() or None,
            }
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Export Failed", str(exc))
            self.statusBar().showMessage("Export failed")
            return

        progress = QtWidgets.QProgressDialog("Preparing BMF export...", None, 0, 100, self)
        progress.setWindowTitle("Export BMF")
        progress.setWindowModality(QtCore.Qt.WindowModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.show()

        thread = QtCore.QThread(self)
        worker = BmfExportWorker(export_kwargs)
        worker.moveToThread(thread)

        self._bmf_export_thread = thread
        self._bmf_export_worker = worker
        self._bmf_export_progress = progress
        self.export_button.setEnabled(False)
        self.export_result_text.setPlainText("Starting BMF export...\n")
        self.statusBar().showMessage("Starting BMF export...")
        print("Starting BMF export...")

        def finish_cleanup() -> None:
            if self._bmf_export_progress is progress:
                progress.close()
                self._bmf_export_progress = None
            self.export_button.setEnabled(True)
            if self._bmf_export_thread is thread:
                self._bmf_export_thread = None
            if self._bmf_export_worker is worker:
                self._bmf_export_worker = None

        def handle_progress(value: int, maximum: int, message: str) -> None:
            progress.setMaximum(max(int(maximum), 1))
            progress.setValue(max(0, min(int(value), int(maximum))))
            if message:
                progress.setLabelText(message)
                self.statusBar().showMessage(message)
                self.export_result_text.appendPlainText(message)

        def handle_finished(result: dict) -> None:
            self.export_result_text.setPlainText(json.dumps(result, indent=2, default=str))
            self.statusBar().showMessage("Export completed")
            print("BMF export completed.")
            if not self.dry_run_check.isChecked() and self.export_output_edit.text().strip():
                self.browse_input_edit.setText(self.export_output_edit.text().strip())
                if not self.browse_output_edit.text().strip():
                    self.browse_output_edit.setText(str(Path(self.export_output_edit.text().strip()).with_suffix(".csv")))
            thread.quit()

        def handle_failed(message: str) -> None:
            QtWidgets.QMessageBox.critical(self, "Export Failed", message)
            self.export_result_text.appendPlainText(f"\nExport failed: {message}")
            self.statusBar().showMessage("Export failed")
            print(f"BMF export failed: {message}")
            thread.quit()

        worker.progress.connect(handle_progress, QtCore.Qt.QueuedConnection)
        worker.finished.connect(handle_finished, QtCore.Qt.QueuedConnection)
        worker.failed.connect(handle_failed, QtCore.Qt.QueuedConnection)
        thread.started.connect(worker.run)
        worker.finished.connect(worker.deleteLater)
        worker.failed.connect(worker.deleteLater)
        thread.finished.connect(finish_cleanup)
        thread.finished.connect(thread.deleteLater)
        thread.start()

    def _save_bmf_csv(self) -> None:
        input_path = self.browse_input_edit.text().strip()
        if not input_path:
            QtWidgets.QMessageBox.warning(self, "Missing File", "Select a BMF file first.")
            return

        output_path = self.browse_output_edit.text().strip()
        if not output_path:
            output_path = str(Path(input_path).with_suffix(".csv")) if input_path else ""
            if output_path:
                self.browse_output_edit.setText(output_path)
        if not output_path:
            QtWidgets.QMessageBox.warning(self, "Missing Output", "Select an output CSV file first.")
            return

        progress = QtWidgets.QProgressDialog("Reading full BMF file for CSV export...", None, 0, 0, self)
        progress.setWindowTitle("Save CSV")
        progress.setWindowModality(QtCore.Qt.WindowModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.show()
        QtWidgets.QApplication.processEvents()

        self.load_bmf_button.setEnabled(False)
        self.save_bmf_csv_button.setEnabled(False)
        self.statusBar().showMessage("Reading full BMF file for CSV export...")

        try:
            result = bmf_tools.load_bmf_table(input_path, row_limit=None)
            rows_frame = result.get("dataframe")
            if not isinstance(rows_frame, pd.DataFrame):
                raise ValueError("This BMF reader mode did not produce tabular row data to save as CSV.")
            progress.setLabelText("Writing CSV file...")
            QtWidgets.QApplication.processEvents()
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            rows_frame.to_csv(output_path, index=False)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Save Failed", str(exc))
            self.statusBar().showMessage("CSV save failed")
            return
        finally:
            progress.close()
            self.load_bmf_button.setEnabled(True)
            self.save_bmf_csv_button.setEnabled(self._loaded_bmf_dataframe is not None)

        self.statusBar().showMessage(f"Saved CSV to {output_path}")
        QtWidgets.QMessageBox.information(self, "CSV Saved", f"Saved {len(rows_frame):,} rows to:\n{output_path}")

    def _load_bmf(self) -> None:
        path = self.browse_input_edit.text().strip()
        if not path:
            QtWidgets.QMessageBox.warning(self, "Missing File", "Select a BMF file first.")
            return
        if self._bmf_load_thread is not None:
            self.statusBar().showMessage("BMF load already in progress")
            return

        row_limit = None if self.read_all_rows_check.isChecked() else int(self.row_limit_spin.value())
        if row_limit is None:
            response = QtWidgets.QMessageBox.warning(
                self,
                "Read All Rows",
                (
                    "Read all rows is selected. Previewing the full BMF can take a long time and may use a lot of memory.\n\n"
                    "Continue loading the full file?"
                ),
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel,
                QtWidgets.QMessageBox.Cancel,
            )
            if response != QtWidgets.QMessageBox.Ok:
                self.statusBar().showMessage("BMF preview canceled")
                return
        self._loaded_bmf_dataframe = None

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
        self.save_bmf_csv_button.setEnabled(False)
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
            self.save_bmf_csv_button.setEnabled(self._loaded_bmf_dataframe is not None)

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
        self._loaded_bmf_dataframe = rows_frame.copy()
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
            window.browse_output_edit.setText(str(candidate.with_suffix(".csv")))

    window.show()
    return int(app.exec_())


if __name__ == "__main__":
    raise SystemExit(main())