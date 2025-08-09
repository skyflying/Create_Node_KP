#Made by Mingyi Hsu

import sys
import os
import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QLabel, QPushButton,
    QVBoxLayout, QWidget, QLineEdit, QCheckBox, QComboBox,
    QMessageBox, QTextEdit, QHBoxLayout, QSizePolicy,
    QTableWidget, QTableWidgetItem, QHeaderView
)
from PyQt5.QtGui import QFont

from processor import process_line_data


class MapCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.updateGeometry()

    def clear(self, message=None):
        self.ax.clear()
        if message:
            self.ax.text(0.5, 0.5, message, ha="center", va="center",
                         transform=self.ax.transAxes, color="gray")
        self.draw()

    def plot_gdf(self, gdf, color="tab:blue"):
        self.ax.clear()
        try:
            gdf.plot(ax=self.ax, edgecolor=color, linewidth=1.0)
            self.ax.set_aspect("equal", adjustable="datalim")
            #self.ax.grid(True, linestyle="--", alpha=0.3)
            self.ax.set_title("Geometry Preview", fontsize=10)
        except Exception as e:
            self.ax.clear()
            self.ax.text(0.5, 0.5, f"Preview failed:\n{e}", ha="center", va="center",
                         transform=self.ax.transAxes, color="crimson")
        self.draw()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Line Node Processor - Created by Mingyi Hsu © 2025")
        self.setGeometry(200, 200, 820, 840)

        self.gis_path = ""
        self.tiff_path = ""
        self.fields = []
        self.gdf = None

        self.init_ui()

    def init_ui(self):
        root = QWidget()
        layout = QVBoxLayout()

        # GIS file selection
        row_gis = QHBoxLayout()
        row_gis.addWidget(QLabel("Select GIS File (.shp / .geojson / .gpkg):"))
        self.gis_path_edit = QLineEdit()
        self.gis_path_edit.setReadOnly(True)
        row_gis.addWidget(self.gis_path_edit)
        btn_gis = QPushButton("Browse")
        btn_gis.clicked.connect(self.load_gis)
        row_gis.addWidget(btn_gis)
        layout.addLayout(row_gis)

        # Info box
        layout.addWidget(QLabel("File Summary"))
        self.info_display = QTextEdit()
        self.info_display.setReadOnly(True)
        self.info_display.setMinimumHeight(130)
        layout.addWidget(self.info_display)

        # --- Attribute headers + first 5 rows (table-like preview) ---
        layout.addWidget(QLabel("Attribute Preview (headers + first 5 rows)"))
        self.preview_table = QTableWidget()
        self.preview_table.setAlternatingRowColors(True)
        self.preview_table.setShowGrid(True)
        # 讓格線淡一點（透明感）
        self.preview_table.setStyleSheet("QTableWidget { gridline-color: rgba(0,0,0,0.2); }")
        self.preview_table.setMinimumHeight(180)
        self.preview_table.horizontalHeader().setStretchLastSection(True)
        self.preview_table.verticalHeader().setVisible(False)
        layout.addWidget(self.preview_table)
        # Inline geometry preview
        layout.addWidget(QLabel("Geometry Preview"))
        self.canvas = MapCanvas(self, width=6, height=4, dpi=100)
        self.canvas.clear("Load a GIS file to preview geometry")
        layout.addWidget(self.canvas)

        # Group by field
        row_group = QHBoxLayout()
        row_group.addWidget(QLabel("Field to Group By:"))
        self.group_combo = QComboBox()
        row_group.addWidget(self.group_combo)
        layout.addLayout(row_group)

        # Fixed distance
        row_dist = QHBoxLayout()
        row_dist.addWidget(QLabel("Fixed Distance (meters):"))
        self.distance_input = QLineEdit("10")
        row_dist.addWidget(self.distance_input)
        layout.addLayout(row_dist)

        # GeoTIFF selection (optional)
        row_tiff = QHBoxLayout()
        row_tiff.addWidget(QLabel("Optional GeoTIFF (for elevation):"))
        self.tiff_path_edit = QLineEdit()
        self.tiff_path_edit.setReadOnly(True)
        row_tiff.addWidget(self.tiff_path_edit)
        btn_tiff = QPushButton("Browse")
        btn_tiff.clicked.connect(self.load_tiff)
        row_tiff.addWidget(btn_tiff)
        layout.addLayout(row_tiff)

        # Options
        self.preserve_nodes_cb = QCheckBox("Preserve Original Nodes")
        self.preserve_nodes_cb.setChecked(True)
        layout.addWidget(self.preserve_nodes_cb)

        self.preserve_attrs_cb = QCheckBox("Preserve Original Attributes")
        self.preserve_attrs_cb.setChecked(True)
        layout.addWidget(self.preserve_attrs_cb)

        self.output_shp_cb = QCheckBox("Output Shapefiles")
        self.output_shp_cb.setChecked(True)
        layout.addWidget(self.output_shp_cb)

        # Run
        btn_run = QPushButton("Run Processing")
        btn_run.clicked.connect(self.run_processing)
        layout.addWidget(btn_run)

        # Footer
        credit = QLabel("Created by Mingyi Hsu © 2025")
        f = QFont()
        f.setPointSize(9)
        credit.setFont(f)
        layout.addWidget(credit)

        root.setLayout(layout)
        self.setCentralWidget(root)

    def load_gis(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select GIS File", "",
            "GIS Files (*.shp *.geojson *.gpkg);;All Files (*)"
        )
        if not path:
            return

        try:
            gdf = gpd.read_file(path)
        except Exception as e:
            QMessageBox.critical(self, "Read Error", f"Failed to read file:\n{e}")
            self.canvas.clear("Failed to read file")
            return

        # Set default CRS if missing
        if gdf.crs is None:
            try:
                gdf.set_crs(epsg=3826, inplace=True)
                crs_info = "None (set to EPSG:3826)"
            except Exception:
                crs_info = "None (failed to set EPSG:3826)"
        else:
            try:
                crs_info = gdf.crs.to_string()
            except Exception:
                crs_info = str(gdf.crs)

        # Warn if geographic CRS (e.g., EPSG:4326)
        try:
            epsg_code = gdf.crs.to_epsg()
        except Exception:
            epsg_code = None

        is_geographic = False
        try:
            is_geographic = getattr(gdf.crs, "is_geographic", False)
        except Exception:
            if epsg_code == 4326:
                is_geographic = True

        if is_geographic or epsg_code == 4326:
            QMessageBox.warning(
                self,
                "Projected CRS Required",
                ("Your input layer appears to use a geographic CRS (e.g., EPSG:4326 / WGS84).\n\n"
                 "Distance-based sampling requires a *projected* coordinate system (planar), "
                 "such as UTM or your local EPSG (meters). Please reproject your data and try again.")
            )

        # Update state
        self.gis_path = path
        self.gdf = gdf
        self.gis_path_edit.setText(path)

        # Populate group fields
        self.fields = list(gdf.columns)
        self.group_combo.clear()
        self.group_combo.addItems(self.fields)

        # Summary info
        geom_counts = gdf.geometry.geom_type.value_counts(dropna=False).to_dict()
        geom_info = ", ".join([f"{k}:{v}" for k, v in geom_counts.items()])
        total_features = len(gdf)

        summary_lines = [
            f"CRS: {crs_info}",
            f"Total Features: {total_features}",
            f"Geometry Types: {geom_info}",
        ]
        self.info_display.setText("\n".join(summary_lines))

        # Attribute preview (first 5 rows, up to 5 columns)
        try:
            cols = [c for c in gdf.columns if c != gdf.geometry.name]  # 去掉幾何欄
            head_df = gdf[cols].head(5)  # 只取前 5 列，但保留全部欄位
        
            # 設定表格大小與表頭
            self.preview_table.clear()
            self.preview_table.setColumnCount(len(cols))
            self.preview_table.setRowCount(len(head_df))
            self.preview_table.setHorizontalHeaderLabels([str(c) for c in cols])
        
            # 填入資料
            for r in range(len(head_df)):
                for c, col_name in enumerate(cols):
                    val = head_df.iloc[r][col_name]
                    # 簡單轉字串，避免 NaN/None 顯示問題
                    s = "" if pd.isna(val) else str(val)
                    self.preview_table.setItem(r, c, QTableWidgetItem(s))
        
            # 欄寬自動，並保留一點空間
            header = self.preview_table.horizontalHeader()
            header.setSectionResizeMode(QHeaderView.ResizeToContents)
            header.setStretchLastSection(True)
        
        except Exception as e:
            # 若表格失敗，保底清空並顯示訊息
            self.preview_table.clear()
            self.preview_table.setRowCount(0)
            self.preview_table.setColumnCount(1)
            self.preview_table.setHorizontalHeaderLabels(["(preview unavailable)"])
        # Inline preview
        try:
            self.canvas.plot_gdf(gdf)
        except Exception as e:
            self.canvas.clear(f"Preview failed:\n{e}")

        # Warn on non-line geometries
        non_line = set(gdf.geometry.geom_type.unique()) - {"LineString", "MultiLineString"}
        if non_line:
            QMessageBox.warning(
                self, "Geometry Warning",
                "This dataset contains non-line geometry. Only LineString / MultiLineString will be processed."
            )

    def load_tiff(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select GeoTIFF File", "", "GeoTIFF (*.tif *.tiff)"
        )
        if path:
            self.tiff_path = path
            self.tiff_path_edit.setText(path)

    def run_processing(self):
        if not self.gis_path:
            QMessageBox.warning(self, "Missing GIS File", "Please select a GIS file first.")
            return

        try:
            fixed_distance = float(self.distance_input.text())
        except ValueError:
            QMessageBox.warning(self, "Invalid Distance", "Please enter a valid number for fixed distance.")
            return

        group_field = self.group_combo.currentText().strip() if self.group_combo.count() else None
        if not group_field:
            QMessageBox.warning(self, "Group Field", "Please select a field to group by.")
            return

        confirm = QMessageBox.question(
            self, "Confirm",
            f"Input: {self.gis_path}\n"
            f"GeoTIFF: {getattr(self, 'tiff_path', None) or '(none)'}\n"
            f"Distance: {fixed_distance}\n"
            f"Group by: {group_field}\n"
            f"Export SHP: {'Yes' if self.output_shp_cb.isChecked() else 'No'}\n\n"
            f"Start processing?",
            QMessageBox.Ok | QMessageBox.Cancel
        )
        if confirm != QMessageBox.Ok:
            return

        ok = process_line_data(
            gis_path=self.gis_path,
            group_field=group_field,
            fixed_distance=fixed_distance,
            output_shp=self.output_shp_cb.isChecked(),
            tiff_path=getattr(self, "tiff_path", None) or None,
            preserve_nodes=self.preserve_nodes_cb.isChecked(),
            preserve_attributes=self.preserve_attrs_cb.isChecked()
        )

        if ok:
            out_dir = os.path.join(os.path.dirname(self.gis_path) or os.getcwd(), "output")
            QMessageBox.information(self, "Done", f"Processing completed.\nOutputs saved in:\n{out_dir}")
        else:
            QMessageBox.critical(self, "Failed", "Processing failed. Check console logs for details.")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())
