#Made by Mingyi Hsu

# gui.py
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import geopandas as gpd
import matplotlib.pyplot as plt
import os
from processor import process_line_data

gdf = None
segment_fields = []

def is_valid_line_geometry(gdf_local):
    allowed = ["LineString", "MultiLineString"]
    return all(g in allowed for g in gdf_local.geometry.geom_type.unique())

def preview_map(gdf_local, field_name=None):
    plt.figure(figsize=(8, 8))
    if field_name and field_name in gdf_local.columns:
        gdf_local.plot(column=field_name, legend=True)
    else:
        gdf_local.plot(edgecolor="blue")
    plt.title("Line Segment Preview")
    plt.axis("equal")
    plt.grid(True)
    plt.show()

def summarize_gdf(gdf_local):
    info = []
    geom_types = gdf_local.geometry.geom_type.value_counts().to_dict()
    total_length = gdf_local.length.sum()
    total_features = len(gdf_local)

    info.append(f"✅ Total line segments: {total_features}")
    info.append("✅ Geometry types:")
    for k, v in geom_types.items():
        info.append(f"  - {k}: {v}")

    if gdf_local.crs is None:
        gdf_local.set_crs("EPSG:3826", inplace=True)
        info.append("⚠️ No CRS found, defaulting to EPSG:3826")
    else:
        info.append(f"✅ CRS: {gdf_local.crs.to_string()}")

    info.append(f"📏 Total Length (2D): {round(total_length, 2)} meters")
    return "\n".join(info)

def open_gis_file():
    global gdf, segment_fields

    filetypes = [("GIS files", "*.shp *.geojson *.gpkg")]
    filepath = filedialog.askopenfilename(title="Select GIS file", filetypes=filetypes)

    if not filepath:
        return

    try:
        gdf = gpd.read_file(filepath)
    except Exception as e:
        messagebox.showerror("Read Error", f"Unable to read file:\n{e}")
        return

    if not is_valid_line_geometry(gdf):
        messagebox.showerror("Geometry Error", "❌ Only LineString or MultiLineString are supported")
        return

    segment_fields = list(gdf.columns)
    field_combo['values'] = segment_fields
    field_combo.set(segment_fields[0])

    summary = summarize_gdf(gdf)
    summary_text.config(state="normal")
    summary_text.delete("1.0", tk.END)
    summary_text.insert(tk.END, summary)
    summary_text.config(state="disabled")

    preview_map(gdf)

    file_path_entry.delete(0, tk.END)
    file_path_entry.insert(0, filepath)

def open_geotiff_file():
    filetypes = [("GeoTIFF files", "*.tif *.tiff")]
    filepath = filedialog.askopenfilename(title="Select GeoTIFF file", filetypes=filetypes)
    geotiff_entry.delete(0, tk.END)
    geotiff_entry.insert(0, filepath)

def run_processing():
    if gdf is None:
        messagebox.showerror("Error", "Please import a GIS file first")
        return

    file_path = file_path_entry.get()
    geotiff_path = geotiff_entry.get() or None
    try:
        distance = float(distance_entry.get())
    except ValueError:
        messagebox.showerror("Error", "Distance must be a number")
        return
    include_original = include_original_var.get()
    retain_attr = retain_attr_var.get()
    export_shp = export_shp_var.get()
    field_name = field_combo.get()

    try:
        process_line_data(
            file_path=file_path,
            geotiff_path=geotiff_path,
            fixed_distance=distance,
            include_original_nodes=include_original,
            retain_attributes=retain_attr,
            export_shapefile=export_shp,
            field_name=field_name
        )
        messagebox.showinfo("Done", "✅ Processing completed! Output saved to 'output' folder.")
    except Exception as e:
        messagebox.showerror("Processing Error", f"Processing failed:\n{e}")

def launch_gui():
    global summary_text, distance_entry, include_original_var, retain_attr_var, export_shp_var, field_combo
    global file_path_entry, geotiff_entry

    root = tk.Tk()
    root.title("GIS Line Processor")
    root.geometry("700x700")

    tk.Label(root, text="📂 GIS File:").pack()
    file_path_entry = tk.Entry(root, width=80)
    file_path_entry.pack()
    tk.Button(root, text="Select GIS File", command=open_gis_file).pack(pady=5)

    tk.Label(root, text="🌄 GeoTIFF File (optional):").pack()
    geotiff_entry = tk.Entry(root, width=80)
    geotiff_entry.pack()
    tk.Button(root, text="Select GeoTIFF File", command=open_geotiff_file).pack(pady=5)

    tk.Label(root, text="📐 Node Spacing (meters, 0 = original nodes only)").pack()
    distance_entry = tk.Entry(root)
    distance_entry.insert(0, "10")
    distance_entry.pack()

    include_original_var = tk.BooleanVar(value=True)
    retain_attr_var = tk.BooleanVar(value=True)
    export_shp_var = tk.BooleanVar(value=True)

    tk.Checkbutton(root, text="✅ Include Original Nodes", variable=include_original_var).pack()
    tk.Checkbutton(root, text="✅ Retain Attributes", variable=retain_attr_var).pack()
    tk.Checkbutton(root, text="✅ Export Shapefile", variable=export_shp_var).pack()

    tk.Label(root, text="📊 Classification Field (segment name field)").pack()
    field_combo = ttk.Combobox(root, state="readonly")
    field_combo.pack(pady=5)

    tk.Button(root, text="▶️ Run Processing", font=("Arial", 14), command=run_processing).pack(pady=10)

    tk.Label(root, text="📋 File Summary", font=("Arial", 12, "bold")).pack()
    summary_text = tk.Text(root, height=15, width=80, wrap="word", state="disabled")
    summary_text.pack(padx=10, pady=10)

    root.mainloop()

if __name__ == "__main__":
    launch_gui()
