#Made by Mingyi Hsu

# processor.py
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point
from pyproj import Transformer
from math import sqrt
import rasterio
import os

def sample_elevation(point, src_array, transform, nodata_value):
    if src_array is None or transform is None:
        return None
    px, py = ~transform * (point.x, point.y)
    px, py = int(px), int(py)
    if 0 <= px < src_array.shape[1] and 0 <= py < src_array.shape[0]:
        value = src_array[py, px]
        if value != nodata_value and value != -9999 and not np.isnan(value):
            return value
    return None

def process_line_data(
    file_path: str,
    geotiff_path: str = None,
    fixed_distance: float = 0,
    include_original_nodes: bool = True,
    retain_attributes: bool = True,
    export_shapefile: bool = True,
    field_name: str = None
):
    gdf = gpd.read_file(file_path)
    if gdf.crs is None:
        gdf.set_crs("EPSG:3826", inplace=True)
    if field_name not in gdf.columns:
        gdf["segment_name"] = gdf.index.astype(str)
        field_name = "segment_name"

    use_geotiff = bool(geotiff_path)
    if use_geotiff:
        with rasterio.open(geotiff_path) as src:
            geotiff_array = src.read(1)
            transform = src.transform
            nodata_value = src.nodata
    else:
        geotiff_array, transform, nodata_value = None, None, None

    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    transformer = Transformer.from_crs(gdf.crs, "EPSG:4326", always_xy=True)
    report_data = []

    for segment_id in gdf[field_name].unique():
        segment_gdf = gdf[gdf[field_name] == segment_id]
        nodes_data = []
        geometry_data = []
        total_2d = 0
        total_3d = 0

        for _, row in segment_gdf.iterrows():
            line = row.geometry
            if line is None:
                continue
            original_points = []
            if line.geom_type == 'LineString':
                original_points = list(line.coords)
            elif line.geom_type == 'MultiLineString':
                for ls in line.geoms:
                    original_points.extend(ls.coords)
            points = set()
            if fixed_distance == 0:
                for pt in original_points:
                    p = Point(pt)
                    d = line.project(p)
                    points.add((p, d))
            else:
                d = 0
                while d <= line.length:
                    p = line.interpolate(d)
                    points.add((p, d))
                    d += fixed_distance
                if include_original_nodes:
                    for pt in original_points:
                        p = Point(pt)
                        d = line.project(p)
                        points.add((p, d))
                if d - fixed_distance < line.length:
                    p = line.interpolate(line.length)
                    points.add((p, line.length))
            sorted_pts = sorted(points, key=lambda x: x[1])
            prev = None
            azimuths = []

            for i, (pt, dist) in enumerate(sorted_pts):
                lon, lat = transformer.transform(pt.x, pt.y)
                elev = sample_elevation(pt, geotiff_array, transform, nodata_value) if use_geotiff else None
                azimuth = None
                dist_2d = None
                dist_3d = None
                if prev:
                    dx = pt.x - prev[0].x
                    dy = pt.y - prev[0].y
                    if dx != 0 or dy != 0:
                        angle = np.arctan2(dx, dy)
                        azimuth = (np.degrees(angle) + 360) % 360
                    dist_2d = round(pt.distance(prev[0]), 3)
                    total_2d += dist_2d
                    if use_geotiff and elev is not None and prev[2] is not None:
                        dz = elev - prev[2]
                        dist_3d = sqrt(dist_2d ** 2 + dz ** 2)
                        total_3d += dist_3d
                azimuths.append(azimuth)
                prev = (pt, lon, lat, elev)

            for i, (pt, dist) in enumerate(sorted_pts):
                lon, lat = transformer.transform(pt.x, pt.y)
                elev = sample_elevation(pt, geotiff_array, transform, nodata_value) if use_geotiff else None
                az = azimuths[i] if i < len(azimuths) else None
                dist_2d = None
                dist_3d = None
                if i > 0:
                    dist_2d = round(pt.distance(sorted_pts[i-1][0]), 3)
                    if use_geotiff and elev is not None:
                        prev_elev = sorted_pts[i-1][0].z if hasattr(sorted_pts[i-1][0], 'z') else None
                        dz = elev - (prev_elev or 0)
                        dist_3d = sqrt(dist_2d ** 2 + dz ** 2)
                node = {
                    "Longitude": lon,
                    "Latitude": lat,
                    "Easting": pt.x,
                    "Northing": pt.y,
                    "Elevation": elev,
                    "Distance": dist_2d,
                    "Length_3D": dist_3d,
                    "Azimuth": round(az, 3) if az is not None else None,
                    "KP": round(dist, 3),
                    "Total_3D_Length": total_3d if use_geotiff else None,
                    field_name: segment_id
                }
                if retain_attributes:
                    for col in gdf.columns:
                        if col not in node:
                            node[col] = row[col]
                nodes_data.append(node)
                geometry_data.append(Point(pt.x, pt.y))

        if nodes_data:
            df = pd.DataFrame(nodes_data)
            csv_path = os.path.join(output_dir, f"{segment_id}_KP.csv")
            df.to_csv(csv_path, index=False)
            shp_path = None
            if export_shapefile:
                geo_df = gpd.GeoDataFrame(df, geometry=geometry_data, crs=gdf.crs)
                shp_path = os.path.join(output_dir, f"{segment_id}_KP.shp")
                geo_df.to_file(shp_path, driver="ESRI Shapefile")
            report_data.append({
                "Segment": segment_id,
                "Total_2D_Distance": total_2d,
                "Total_3D_Distance": total_3d if use_geotiff else None,
                "CSV": csv_path,
                "Shapefile": shp_path
            })

    pd.DataFrame(report_data).to_csv(os.path.join(output_dir, "processing_report.csv"), index=False)
