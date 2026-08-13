#!/usr/bin/env python3
"""Standalone GUI for drillhole limit interval merging and behavior analysis."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
from PyQt5 import QtCore, QtWidgets

from drillhole_limit_analysis import BETWEEN_BOUNDARIES_LABEL, analyze_drillhole_limit_behavior


_VIEWER_CSV_HELPERS: tuple[Any, Any] | None = None


def _display_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if pd.isna(value):
            return ""
        return f"{value:g}"
    return str(value)


def _coerce_numeric_series(series: pd.Series) -> pd.Series:
    numeric_series = pd.to_numeric(series, errors="coerce")
    missing_mask = numeric_series.isna()
    if not missing_mask.any():
        return numeric_series

    text_series = series.astype(str)
    normalized = text_series.str.replace("\u00a0", "", regex=False).str.replace(" ", "", regex=False)
    alternate_numeric = pd.to_numeric(normalized, errors="coerce")
    return numeric_series.where(~missing_mask, alternate_numeric)


def _guess_column(columns: list[str], tokens: tuple[str, ...]) -> str:
    lowered = {column: column.strip().lower() for column in columns}
    for token in tokens:
        for column, normalized in lowered.items():
            if normalized == token:
                return column
    for token in tokens:
        for column, normalized in lowered.items():
            if token in normalized:
                return column
    return ""


def _detect_numeric_columns(frame: pd.DataFrame) -> list[str]:
    numeric_columns: list[str] = []
    for column in frame.columns:
        series = frame[column]
        numeric_values = _coerce_numeric_series(series)
        valid_count = int(numeric_values.notna().sum())
        if valid_count == 0:
            continue
        if pd.api.types.is_numeric_dtype(series):
            numeric_columns.append(str(column))
            continue
        if valid_count >= max(3, int(len(frame) * 0.6)):
            numeric_columns.append(str(column))
    return numeric_columns


def _get_viewer_csv_helpers() -> tuple[Any, Any]:
    global _VIEWER_CSV_HELPERS
    if _VIEWER_CSV_HELPERS is None:
        import anterpolator3DViewer as viewer_app

        _VIEWER_CSV_HELPERS = (viewer_app.read_autodetect_csv, viewer_app.load_csv_preview_dataframe)
    return _VIEWER_CSV_HELPERS


def _read_drillhole_csv(csv_path: Path, separator: str | None) -> pd.DataFrame:
    read_autodetect_csv, _load_csv_preview_dataframe = _get_viewer_csv_helpers()
    return read_autodetect_csv(csv_path, min_cols=1, forced_delimiter=separator)


class DataFrameTableModel(QtCore.QAbstractTableModel):
    def __init__(self, frame: pd.DataFrame | None = None):
        super().__init__()
        self._frame = frame if frame is not None else pd.DataFrame()

    @property
    def frame(self) -> pd.DataFrame:
        return self._frame

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
        return _display_value(self._frame.iat[index.row(), index.column()])

    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = QtCore.Qt.DisplayRole) -> Any:
        if role != QtCore.Qt.DisplayRole:
            return None
        if orientation == QtCore.Qt.Horizontal:
            if 0 <= section < len(self._frame.columns):
                return str(self._frame.columns[section])
            return ""
        return str(section)


class DrillholeAnalysisWorker(QtCore.QObject):
    progress = QtCore.pyqtSignal(int, int, str)
    finished = QtCore.pyqtSignal(dict)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, analysis_kwargs: dict[str, Any]):
        super().__init__()
        self.analysis_kwargs = dict(analysis_kwargs)

    @QtCore.pyqtSlot()
    def run(self) -> None:
        try:
            result = analyze_drillhole_limit_behavior(
                progress_callback=lambda current, total, message: self.progress.emit(int(current), int(total), str(message)),
                **self.analysis_kwargs,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


class DrillholeLimitAnalysisWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Drillhole Limit Analysis")
        self.resize(1440, 920)

        self.source_frame: pd.DataFrame | None = None
        self.analysis_result: dict[str, pd.DataFrame] | None = None
        self.numeric_candidates: list[str] = []
        self.categorical_candidates: list[str] = []
        self._analysis_thread: QtCore.QThread | None = None
        self._analysis_worker: DrillholeAnalysisWorker | None = None
        self._analysis_progress_dialog: QtWidgets.QProgressDialog | None = None
        self._analysis_started_at: float | None = None

        self.preview_model = DataFrameTableModel()
        self.quality_model = DataFrameTableModel()
        self.merged_model = DataFrameTableModel()
        self.summary_model = DataFrameTableModel()
        self.group_model = DataFrameTableModel()
        self.tests_model = DataFrameTableModel()
        self.boundary_model = DataFrameTableModel()

        self._build_ui()
        self.progress_status_label = QtWidgets.QLabel("Idle")
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFixedWidth(260)
        self.statusBar().addPermanentWidget(self.progress_status_label, 1)
        self.statusBar().addPermanentWidget(self.progress_bar)
        self.statusBar().showMessage("Load a CSV to begin.")

    def _build_ui(self) -> None:
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        help_label = QtWidgets.QLabel(
            "Selected boundary categories stay as standalone intervals. All intervals between those boundaries are merged per hole using interval-length weighted averages for numeric fields and weighted-majority values for selected categorical fields."
        )
        help_label.setWordWrap(True)
        layout.addWidget(help_label)

        source_group = QtWidgets.QGroupBox("Source CSV")
        source_form = QtWidgets.QFormLayout(source_group)
        self.csv_path_edit = QtWidgets.QLineEdit()
        self.delimiter_combo = QtWidgets.QComboBox()
        self.delimiter_combo.addItem("Auto", None)
        self.delimiter_combo.addItem("Comma ,", ",")
        self.delimiter_combo.addItem("Semicolon ;", ";")
        self.delimiter_combo.addItem("Tab", "\t")
        self.delimiter_combo.addItem("Pipe |", "|")
        self.load_csv_button = QtWidgets.QPushButton("Load CSV")
        self.load_csv_button.clicked.connect(self._load_csv)
        source_form.addRow("CSV file", self._path_row(self.csv_path_edit, self._browse_csv))
        source_form.addRow("Delimiter", self.delimiter_combo)
        source_form.addRow("", self.load_csv_button)
        layout.addWidget(source_group)

        mapping_group = QtWidgets.QGroupBox("Analysis Setup")
        mapping_layout = QtWidgets.QVBoxLayout(mapping_group)
        mapping_form = QtWidgets.QFormLayout()
        self.hole_id_combo = QtWidgets.QComboBox()
        self.from_combo = QtWidgets.QComboBox()
        self.to_combo = QtWidgets.QComboBox()
        self.limit_combo = QtWidgets.QComboBox()
        self.limit_combo.currentTextChanged.connect(self._refresh_boundary_categories)
        self.min_group_size_spin = QtWidgets.QSpinBox()
        self.min_group_size_spin.setRange(1, 999)
        self.min_group_size_spin.setValue(3)
        self.skip_holes_without_boundary_check = QtWidgets.QCheckBox("Skip drillholes without selected boundary")
        self.skip_holes_without_boundary_check.setChecked(True)
        mapping_form.addRow("Hole ID", self.hole_id_combo)
        mapping_form.addRow("From", self.from_combo)
        mapping_form.addRow("To", self.to_combo)
        mapping_form.addRow("Limit field", self.limit_combo)
        mapping_form.addRow("Min group size", self.min_group_size_spin)
        mapping_form.addRow("", self.skip_holes_without_boundary_check)
        mapping_layout.addLayout(mapping_form)

        selector_row = QtWidgets.QHBoxLayout()
        self.numeric_list = self._make_multi_select_list()
        self.composite_list = self._make_multi_select_list()
        self.boundary_list = self._make_multi_select_list()
        selector_row.addWidget(self._grouped_list("Numeric fields", self.numeric_list))
        selector_row.addWidget(self._grouped_list("Composite categorical fields", self.composite_list))
        selector_row.addWidget(self._grouped_list("Boundary categories", self.boundary_list))
        mapping_layout.addLayout(selector_row)

        action_row = QtWidgets.QHBoxLayout()
        self.load_config_button = QtWidgets.QPushButton("Load Config")
        self.load_config_button.clicked.connect(self._load_config)
        self.save_config_button = QtWidgets.QPushButton("Save Config")
        self.save_config_button.clicked.connect(self._save_config)
        self.run_button = QtWidgets.QPushButton("Run Analysis")
        self.run_button.clicked.connect(self._run_analysis)
        self.export_button = QtWidgets.QPushButton("Export Current Table")
        self.export_button.clicked.connect(self._export_current_table)
        action_row.addWidget(self.load_config_button)
        action_row.addWidget(self.save_config_button)
        action_row.addWidget(self.run_button)
        action_row.addWidget(self.export_button)
        action_row.addStretch(1)
        mapping_layout.addLayout(action_row)
        layout.addWidget(mapping_group)

        self.results_tabs = QtWidgets.QTabWidget()
        self.preview_view = self._make_table_view(self.preview_model)
        self.quality_view = self._make_table_view(self.quality_model)
        self.merged_view = self._make_table_view(self.merged_model)
        self.summary_view = self._make_table_view(self.summary_model)
        self.group_view = self._make_table_view(self.group_model)
        self.tests_view = self._make_table_view(self.tests_model)
        self.boundary_view = self._make_table_view(self.boundary_model)
        self.results_tabs.addTab(self.preview_view, "Input Preview")
        self.results_tabs.addTab(self.quality_view, "Quality")
        self.results_tabs.addTab(self.merged_view, "Merged Intervals")
        self.results_tabs.addTab(self.summary_view, "Summary")
        self.results_tabs.addTab(self.group_view, "Group Summary")
        self.results_tabs.addTab(self.tests_view, "Statistical Tests")
        self.results_tabs.addTab(self.boundary_view, "Boundary Transitions")
        layout.addWidget(self.results_tabs, stretch=1)

    def _make_multi_select_list(self) -> QtWidgets.QListWidget:
        widget = QtWidgets.QListWidget()
        widget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        return widget

    def _grouped_list(self, title: str, widget: QtWidgets.QListWidget) -> QtWidgets.QGroupBox:
        group = QtWidgets.QGroupBox(title)
        group_layout = QtWidgets.QVBoxLayout(group)
        group_layout.addWidget(widget)
        return group

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

    def _browse_csv(self) -> None:
        path, _selected_filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Open Drillhole CSV",
            str(Path(self.csv_path_edit.text()).parent if self.csv_path_edit.text().strip() else Path.cwd()),
            "CSV files (*.csv *.txt);;All files (*.*)",
        )
        if path:
            self.csv_path_edit.setText(path)

    def _load_csv(self) -> None:
        path_text = self.csv_path_edit.text().strip()
        if not path_text:
            QtWidgets.QMessageBox.warning(self, "Missing Path", "Choose a CSV file first.")
            return

        csv_path = Path(path_text)
        if not csv_path.exists():
            QtWidgets.QMessageBox.warning(self, "Missing File", f"CSV file not found:\n{csv_path}")
            return

        separator = self.delimiter_combo.currentData()
        try:
            frame = _read_drillhole_csv(csv_path, separator)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Load Failed", str(exc))
            return

        self.source_frame = frame
        self.preview_model.set_frame(frame.head(500))
        self._populate_field_controls(frame)
        self.statusBar().showMessage(f"Loaded {len(frame):,} rows and {len(frame.columns):,} columns from {csv_path.name}.")
        self.progress_status_label.setText("CSV loaded")

    def _populate_field_controls(self, frame: pd.DataFrame) -> None:
        columns = [str(column) for column in frame.columns]
        self.numeric_candidates = _detect_numeric_columns(frame)
        numeric_set = set(self.numeric_candidates)
        self.categorical_candidates = [column for column in columns if column not in numeric_set]

        for combo in [self.hole_id_combo, self.from_combo, self.to_combo, self.limit_combo]:
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(columns)
            combo.blockSignals(False)

        self._fill_list(self.numeric_list, self.numeric_candidates)
        composite_candidates = [column for column in self.categorical_candidates if column not in {self.hole_id_combo.currentText(), self.limit_combo.currentText()}]
        self._fill_list(self.composite_list, composite_candidates)

        self._set_combo_value(self.hole_id_combo, _guess_column(columns, ("hole_id", "holeid", "dhid", "bhid", "id")))
        self._set_combo_value(self.from_combo, _guess_column(columns, ("from", "depth_from", "from_m", "desde")))
        self._set_combo_value(self.to_combo, _guess_column(columns, ("to", "depth_to", "to_m", "hasta")))

        limit_guess = _guess_column(
            self.categorical_candidates,
            ("fault_zone", "fault", "domain", "zone", "limit", "structure"),
        )
        self._set_combo_value(self.limit_combo, limit_guess or (self.categorical_candidates[0] if self.categorical_candidates else ""))

        selected_numeric = [column for column in self.numeric_candidates if column not in {self.from_combo.currentText(), self.to_combo.currentText()}]
        self._set_selected_items(self.numeric_list, selected_numeric)
        self._fill_list(
            self.composite_list,
            [column for column in self.categorical_candidates if column not in {self.hole_id_combo.currentText(), self.limit_combo.currentText()}],
        )
        self._refresh_boundary_categories()

    def _fill_list(self, widget: QtWidgets.QListWidget, values: list[str]) -> None:
        current_selection = {item.text() for item in widget.selectedItems()}
        widget.clear()
        for value in values:
            item = QtWidgets.QListWidgetItem(value)
            widget.addItem(item)
            if value in current_selection:
                item.setSelected(True)

    def _set_combo_value(self, combo: QtWidgets.QComboBox, value: str) -> None:
        if not value:
            if combo.count() > 0:
                combo.setCurrentIndex(0)
            return
        index = combo.findText(value)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _set_selected_items(self, widget: QtWidgets.QListWidget, values: list[str]) -> None:
        selected = set(values)
        for row_index in range(widget.count()):
            item = widget.item(row_index)
            item.setSelected(item.text() in selected)

    def _selected_item_texts(self, widget: QtWidgets.QListWidget) -> list[str]:
        return [item.text() for item in widget.selectedItems()]

    def _refresh_boundary_categories(self) -> None:
        if self.source_frame is None:
            self.boundary_list.clear()
            return

        limit_field = self.limit_combo.currentText().strip()
        composite_fields = [column for column in self.categorical_candidates if column not in {self.hole_id_combo.currentText(), limit_field}]
        self._fill_list(self.composite_list, composite_fields)

        if not limit_field or limit_field not in self.source_frame.columns:
            self.boundary_list.clear()
            return

        series = self.source_frame[limit_field].fillna("").astype(str).str.strip()
        unique_values = sorted(value for value in series.unique().tolist() if value)
        self._fill_list(self.boundary_list, unique_values)

        if unique_values:
            suggested_boundaries = [
                value
                for value in unique_values
                if any(token in value.lower() for token in ("fault", "dyke", "shear", "zone", "contact"))
            ]
            self._set_selected_items(self.boundary_list, suggested_boundaries)

    def _run_analysis(self) -> None:
        if self._analysis_thread is not None:
            QtWidgets.QMessageBox.information(self, "Analysis Running", "Wait for the current analysis to finish.")
            return

        if self.source_frame is None:
            QtWidgets.QMessageBox.warning(self, "No Data", "Load a CSV before running the analysis.")
            return

        hole_id_col = self.hole_id_combo.currentText().strip()
        from_col = self.from_combo.currentText().strip()
        to_col = self.to_combo.currentText().strip()
        limit_field = self.limit_combo.currentText().strip()
        numeric_fields = self._selected_item_texts(self.numeric_list)
        composite_fields = self._selected_item_texts(self.composite_list)
        boundary_categories = self._selected_item_texts(self.boundary_list)

        if not all([hole_id_col, from_col, to_col, limit_field]):
            QtWidgets.QMessageBox.warning(self, "Missing Fields", "Select the hole, from, to, and limit fields.")
            return
        if not numeric_fields:
            QtWidgets.QMessageBox.warning(self, "Missing Numeric Fields", "Select at least one numeric field to analyze.")
            return

        analysis_kwargs = {
            "df": self.source_frame,
            "hole_id_col": hole_id_col,
            "from_col": from_col,
            "to_col": to_col,
            "limit_fields": [limit_field],
            "numeric_fields": numeric_fields,
            "min_group_size": int(self.min_group_size_spin.value()),
            "boundary_categories_by_field": {limit_field: boundary_categories} if boundary_categories else None,
            "composite_categorical_fields": composite_fields,
            "skip_holes_without_boundary": bool(self.skip_holes_without_boundary_check.isChecked()),
        }
        self._start_analysis_worker(analysis_kwargs)

    def _start_analysis_worker(self, analysis_kwargs: dict[str, Any]) -> None:
        self._analysis_started_at = time.monotonic()
        self._set_analysis_running_state(True)
        self.progress_status_label.setText("Running analysis")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setValue(0)

        dialog = QtWidgets.QProgressDialog("Preparing analysis...", None, 0, 0, self)
        dialog.setWindowTitle("Running Analysis")
        dialog.setMinimumDuration(0)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)
        dialog.setCancelButton(None)
        dialog.setWindowModality(QtCore.Qt.WindowModal)
        dialog.show()
        self._analysis_progress_dialog = dialog

        self._analysis_thread = QtCore.QThread(self)
        self._analysis_worker = DrillholeAnalysisWorker(analysis_kwargs)
        self._analysis_worker.moveToThread(self._analysis_thread)
        self._analysis_thread.started.connect(self._analysis_worker.run)
        self._analysis_worker.progress.connect(self._update_analysis_progress)
        self._analysis_worker.finished.connect(self._handle_analysis_finished)
        self._analysis_worker.failed.connect(self._handle_analysis_failed)
        self._analysis_worker.finished.connect(self._analysis_thread.quit)
        self._analysis_worker.failed.connect(self._analysis_thread.quit)
        self._analysis_thread.finished.connect(self._cleanup_analysis_thread)
        self._analysis_thread.start()

    def _set_analysis_running_state(self, running: bool) -> None:
        for widget in [
            self.load_csv_button,
            self.load_config_button,
            self.save_config_button,
            self.run_button,
            self.export_button,
            self.csv_path_edit,
            self.delimiter_combo,
            self.hole_id_combo,
            self.from_combo,
            self.to_combo,
            self.limit_combo,
            self.min_group_size_spin,
            self.skip_holes_without_boundary_check,
            self.numeric_list,
            self.composite_list,
            self.boundary_list,
        ]:
            widget.setEnabled(not running)
        self.progress_bar.setVisible(running)

    def _format_remaining_time(self, current: int, total: int) -> str:
        if self._analysis_started_at is None or current <= 0 or total <= 0 or current >= total:
            return "estimating..." if current <= 0 else "less than 1s"
        elapsed_seconds = max(0.001, time.monotonic() - self._analysis_started_at)
        remaining_seconds = max(0.0, elapsed_seconds * (total - current) / current)
        if remaining_seconds < 1.0:
            return "less than 1s"
        if remaining_seconds < 60.0:
            return f"{int(round(remaining_seconds))}s"
        minutes = int(remaining_seconds // 60)
        seconds = int(round(remaining_seconds % 60))
        return f"{minutes}m {seconds:02d}s"

    @QtCore.pyqtSlot(int, int, str)
    def _update_analysis_progress(self, current: int, total: int, message: str) -> None:
        bounded_total = max(1, int(total))
        bounded_current = max(0, min(int(current), bounded_total))
        remaining_text = self._format_remaining_time(bounded_current, bounded_total)
        progress_text = f"{bounded_current:,}/{bounded_total:,} | ETA {remaining_text}"
        self.progress_status_label.setText(progress_text)
        self.progress_bar.setRange(0, bounded_total)
        self.progress_bar.setValue(bounded_current)
        self.progress_bar.setFormat(progress_text)
        self.statusBar().showMessage(message)
        if self._analysis_progress_dialog is not None:
            self._analysis_progress_dialog.setRange(0, bounded_total)
            self._analysis_progress_dialog.setValue(bounded_current)
            self._analysis_progress_dialog.setLabelText(f"{message}\n{progress_text}")

    @QtCore.pyqtSlot(dict)
    def _handle_analysis_finished(self, result: dict[str, pd.DataFrame]) -> None:
        self.analysis_result = result
        self.quality_model.set_frame(result.get("quality_summary"))
        self.merged_model.set_frame(result.get("merged_intervals"))
        self.summary_model.set_frame(result.get("summary"))
        self.group_model.set_frame(result.get("group_summary"))
        self.tests_model.set_frame(result.get("statistical_tests"))
        self.boundary_model.set_frame(result.get("boundary_transitions"))

        quality_counts = dict(zip(result.get("quality_summary", pd.DataFrame()).get("metric", []), result.get("quality_summary", pd.DataFrame()).get("value", [])))
        limit_field = self.limit_combo.currentText().strip()
        skipped_holes = int(quality_counts.get(f"holes_skipped_without_boundary_{limit_field}", 0) or 0)
        merged_frame = result.get("merged_intervals", pd.DataFrame())
        boundary_mode = "boundary-aware merged intervals"
        if not self._selected_item_texts(self.boundary_list):
            boundary_mode = "contiguous category runs"
        self.progress_status_label.setText("Complete")
        self.progress_bar.setFormat("Done")
        self.statusBar().showMessage(
            f"Analysis complete using {boundary_mode}. Produced {len(merged_frame):,} merged intervals and skipped {skipped_holes:,} drillholes without the selected boundary. Category '{BETWEEN_BOUNDARIES_LABEL}' marks composite intervals between boundaries."
        )
        self.results_tabs.setCurrentWidget(self.merged_view)
        if self._analysis_progress_dialog is not None:
            self._analysis_progress_dialog.setValue(self.progress_bar.maximum())
            self._analysis_progress_dialog.close()

    @QtCore.pyqtSlot(str)
    def _handle_analysis_failed(self, message: str) -> None:
        self.progress_status_label.setText("Failed")
        self.statusBar().showMessage("Analysis failed")
        if self._analysis_progress_dialog is not None:
            self._analysis_progress_dialog.close()
        QtWidgets.QMessageBox.critical(self, "Analysis Failed", message)

    @QtCore.pyqtSlot()
    def _cleanup_analysis_thread(self) -> None:
        if self._analysis_worker is not None:
            self._analysis_worker.deleteLater()
        if self._analysis_thread is not None:
            self._analysis_thread.deleteLater()
        self._analysis_worker = None
        self._analysis_thread = None
        self._analysis_progress_dialog = None
        self._analysis_started_at = None
        self._set_analysis_running_state(False)

    def _config_dict(self) -> dict[str, Any]:
        return {
            "config_type": "drillhole_limit_analysis",
            "csv_path": self.csv_path_edit.text().strip(),
            "delimiter": self.delimiter_combo.currentData(),
            "hole_id_col": self.hole_id_combo.currentText().strip(),
            "from_col": self.from_combo.currentText().strip(),
            "to_col": self.to_combo.currentText().strip(),
            "limit_field": self.limit_combo.currentText().strip(),
            "min_group_size": int(self.min_group_size_spin.value()),
            "skip_holes_without_boundary": bool(self.skip_holes_without_boundary_check.isChecked()),
            "numeric_fields": self._selected_item_texts(self.numeric_list),
            "composite_categorical_fields": self._selected_item_texts(self.composite_list),
            "boundary_categories": self._selected_item_texts(self.boundary_list),
        }

    def _apply_config_dict(self, config: dict[str, Any]) -> None:
        self.csv_path_edit.setText(str(config.get("csv_path", "") or ""))
        delimiter_value = config.get("delimiter", None)
        delimiter_index = self.delimiter_combo.findData(delimiter_value)
        if delimiter_index >= 0:
            self.delimiter_combo.setCurrentIndex(delimiter_index)

        csv_path = Path(self.csv_path_edit.text().strip()) if self.csv_path_edit.text().strip() else None
        if csv_path is not None and csv_path.exists():
            self._load_csv()

        if self.source_frame is None:
            self.min_group_size_spin.setValue(int(config.get("min_group_size", 3) or 3))
            self.skip_holes_without_boundary_check.setChecked(bool(config.get("skip_holes_without_boundary", True)))
            return

        self._set_combo_value(self.hole_id_combo, str(config.get("hole_id_col", "") or ""))
        self._set_combo_value(self.from_combo, str(config.get("from_col", "") or ""))
        self._set_combo_value(self.to_combo, str(config.get("to_col", "") or ""))
        self._set_combo_value(self.limit_combo, str(config.get("limit_field", "") or ""))
        self.min_group_size_spin.setValue(int(config.get("min_group_size", 3) or 3))
        self.skip_holes_without_boundary_check.setChecked(bool(config.get("skip_holes_without_boundary", True)))
        self._refresh_boundary_categories()
        self._set_selected_items(self.numeric_list, [str(value) for value in config.get("numeric_fields", [])])
        self._set_selected_items(self.composite_list, [str(value) for value in config.get("composite_categorical_fields", [])])
        self._set_selected_items(self.boundary_list, [str(value) for value in config.get("boundary_categories", [])])

    def _save_config(self) -> None:
        suggested_path = Path(self.csv_path_edit.text().strip()).with_suffix(".limit_analysis.json") if self.csv_path_edit.text().strip() else Path.cwd() / "drillhole_limit_analysis_config.json"
        save_path, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Config",
            str(suggested_path),
            "JSON files (*.json)",
        )
        if not save_path:
            return
        try:
            Path(save_path).write_text(json.dumps(self._config_dict(), indent=2), encoding="utf-8")
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Save Config Failed", str(exc))
            return
        self.statusBar().showMessage(f"Saved config to {save_path}")

    def _load_config(self) -> None:
        start_dir = str(Path(self.csv_path_edit.text()).parent) if self.csv_path_edit.text().strip() else str(Path.cwd())
        config_path, _selected_filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Load Config",
            start_dir,
            "JSON files (*.json)",
        )
        if not config_path:
            return
        try:
            config = json.loads(Path(config_path).read_text(encoding="utf-8"))
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Load Config Failed", str(exc))
            return
        self._apply_config_dict(dict(config or {}))
        self.statusBar().showMessage(f"Loaded config from {config_path}")

    def _current_model(self) -> tuple[str, DataFrameTableModel] | None:
        current_widget = self.results_tabs.currentWidget()
        mapping = {
            self.preview_view: ("input_preview", self.preview_model),
            self.quality_view: ("quality_summary", self.quality_model),
            self.merged_view: ("merged_intervals", self.merged_model),
            self.summary_view: ("summary", self.summary_model),
            self.group_view: ("group_summary", self.group_model),
            self.tests_view: ("statistical_tests", self.tests_model),
            self.boundary_view: ("boundary_transitions", self.boundary_model),
        }
        return mapping.get(current_widget)

    def _export_current_table(self) -> None:
        active = self._current_model()
        if active is None:
            return
        table_name, model = active
        if model.frame.empty:
            QtWidgets.QMessageBox.information(self, "Nothing To Export", "The current table is empty.")
            return

        base_name = Path(self.csv_path_edit.text().strip()).stem or "drillhole_limit_analysis"
        suggested_path = Path.cwd() / f"{base_name}_{table_name}.csv"
        save_path, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export Current Table",
            str(suggested_path),
            "CSV files (*.csv)",
        )
        if not save_path:
            return

        try:
            model.frame.to_csv(save_path, index=False)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Export Failed", str(exc))
            return
        self.statusBar().showMessage(f"Exported {table_name} to {save_path}")


def main(argv: list[str] | None = None) -> int:
    app = QtWidgets.QApplication(argv or sys.argv)
    window = DrillholeLimitAnalysisWindow()

    args = list(sys.argv[1:] if argv is None else argv[1:])
    if args:
        candidate = Path(args[0])
        if candidate.exists():
            window.csv_path_edit.setText(str(candidate))
            window._load_csv()

    window.show()
    return int(app.exec_())


if __name__ == "__main__":
    raise SystemExit(main())