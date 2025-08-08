import sys
import os
import geopandas as gpd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QFileDialog, QLabel, QPushButton,
                             QVBoxLayout, QWidget, QLineEdit, QCheckBox, QComboBox, QMessageBox, QTextEdit)
from PyQt5.QtGui import QFont
from processor import process_line_data

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Create Node from Line Processor")
        self.setGeometry(200, 200, 600, 500)

        self.gis_path = ""
        self.tiff_path = ""
        self.fields = []

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        self.gis_label = QLabel("Select GIS File:")
        layout.addWidget(self.gis_label)
        self.gis_button = QPushButton("Browse")
        self.gis_button.clicked.connect(self.load_gis)
        layout.addWidget(self.gis_button)

        self.info_display = QTextEdit()
        self.info_display.setReadOnly(True)
        layout.addWidget(self.info_display)

        self.group_label = QLabel("Field to Group By:")
        layout.addWidget(self.group_label)
        self.group_combo = QComboBox()
        layout.addWidget(self.group_combo)

        self.distance_label = QLabel("Fixed Distance (meters):")
        layout.addWidget(self.distance_label)
        self.distance_input = QLineEdit("10")
        layout.addWidget(self.distance_input)

        self.tiff_label = QLabel("Optional GeoTIFF for Elevation:")
        layout.addWidget(self.tiff_label)
        self.tiff_button = QPushButton("Browse")
        self.tiff_button.clicked.connect(self.load_tiff)
        layout.addWidget(self.tiff_button)

        self.preserve_nodes_cb = QCheckBox("Preserve Original Nodes")
        self.preserve_nodes_cb.setChecked(True)
        layout.addWidget(self.preserve_nodes_cb)

        self.preserve_attrs_cb = QCheckBox("Preserve Original Attributes")
        self.preserve_attrs_cb.setChecked(True)
        layout.addWidget(self.preserve_attrs_cb)

        self.output_shp_cb = QCheckBox("Output Shapefiles")
        self.output_shp_cb.setChecked(True)
        layout.addWidget(self.output_shp_cb)

        self.run_button = QPushButton("Run Processing")
        self.run_button.clicked.connect(self.run_processing)
        layout.addWidget(self.run_button)

        self.credit = QLabel("Created by Mingyi Hsu © 2025")
        font = QFont()
        font.setPointSize(9)
        self.credit.setFont(font)
        layout.addWidget(self.credit)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_gis(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select GIS File", "", "Shapefile (*.shp);;GeoJSON (*.geojson);;All Files (*)")
        if file_path:
            self.gis_path = file_path
            gdf = gpd.read_file(file_path)
            self.fields = list(gdf.columns)
            self.group_combo.clear()
            self.group_combo.addItems(self.fields)
            info = f"CRS: {gdf.crs}\nTotal Features: {len(gdf)}\nPreview Fields: {self.fields[:5]}"
            self.info_display.setText(info)

    def load_tiff(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select GeoTIFF File", "", "TIFF Files (*.tif *.tiff)")
        if file_path:
            self.tiff_path = file_path

    def run_processing(self):
        if not self.gis_path:
            QMessageBox.warning(self, "Missing GIS File", "Please select a GIS file to process.")
            return

        try:
            fixed_distance = float(self.distance_input.text())
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number for fixed distance.")
            return

        group_field = self.group_combo.currentText()
        success = process_line_data(
            gis_path=self.gis_path,
            group_field=group_field,
            fixed_distance=fixed_distance,
            output_shp=self.output_shp_cb.isChecked(),
            tiff_path=self.tiff_path if self.tiff_path else None,
            preserve_nodes=self.preserve_nodes_cb.isChecked(),
            preserve_attributes=self.preserve_attrs_cb.isChecked()
        )

        if success:
            QMessageBox.information(self, "Done", "Processing completed. Outputs saved in 'output/' folder.")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
