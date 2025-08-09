# Create_Node_KP

A user-friendly  GUI tool (built with Python) that allows you to process vector line data (e.g., roads, pipelines) and optionally extract elevation from GeoTIFF files. It supports creating interpolated nodes or create the original vertex along line segments and exporting the results as CSV and Shapefile formats.

---

## ✨ Features

- Load vector line files (`.shp`, `.geojson`, `.gpkg`)
- Optional elevation sampling from a GeoTIFF
- Customize node spacing along lines
- Export processed data to CSV and Shapefile
- Automatically computes:
  - Azimuths
  - KP (distance along line)
  - Elevation (if GeoTIFF provided)
  - 2D/3D distances between nodes
- Preview lines with matplotlib before processing

---
## ✨ Function

- Line Listing
- RPL(Route Position Line)
- Line Segement
- Cable length
