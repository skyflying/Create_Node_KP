import os
import math
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString
from pyproj import Transformer

# --- Rounding policy: by fixed_distance ---
def round_by_distance(value: float, step: float):
    if value is None:
        return None
    try:
        s = float(step)
    except Exception:
        s = 0.0
    if s < 1:
        return round(value, 2)   # e.g., 0.5 -> 2 decimals
    elif s < 10:
        return round(value, 1)   # e.g., 1 -> 1 decimal
    else:
        return round(value)      # integer

def _azimuth(p_prev: Point, p_cur: Point):
    if p_prev is None:
        return None
    dx = p_cur.x - p_prev.x
    dy = p_cur.y - p_prev.y
    ang = math.degrees(math.atan2(dx, dy))
    return (ang + 360.0) % 360.0

def _sanitize_name(v):
    s = str(v).strip()
    for ch in r'\/:*?"<>|':
        s = s.replace(ch, "_")
    return s or "group"

def _identity_xy(x, y):
    return x, y

def _build_to_raster_transformer(gdf_crs, ds):
    try:
        raster_crs = getattr(ds, "crs", None)
    except Exception:
        raster_crs = None

    if raster_crs is None:
        # No CRS on raster → assume same as GIS
        return _identity_xy

    if str(raster_crs) == str(gdf_crs):
        return _identity_xy

    tr = Transformer.from_crs(gdf_crs, raster_crs, always_xy=True)
    return tr.transform

def _within_bounds(ds, x, y) -> bool:
    b = ds.bounds  # left, bottom, right, top
    return (b.left <= x <= b.right) and (b.bottom <= y <= b.top)

def _sample_elev_lazy(ds, x, y, to_raster_xy):
    """
    x, y in GIS CRS → transform to raster CRS (if needed) → sample.
    Returns float or None.
    """
    try:
        if ds is None:
            return None

        rx, ry = to_raster_xy(x, y)
        if not _within_bounds(ds, rx, ry):
            return None

        row, col = ds.index(rx, ry)
        if row < 0 or col < 0 or row >= ds.height or col >= ds.width:
            return None

        val = ds.read(1)[row, col]
        # Common NoData handling
        if val is None:
            return None
        if isinstance(val, float):
            if math.isnan(val) or val >= 3.4e38:
                return None
        return float(val)
    except Exception:
        return None

# --- Safe step and sampling along a line ---
def _safe_step(step: float, L: float) -> float:
    """Ensure step is positive; if <=0 use whole length (only endpoints)."""
    try:
        s = float(step)
    except Exception:
        s = 0.0
    if s <= 0:
        s = L
    return s

def _sample_points_along_line(line: LineString, step: float, preserve_nodes: bool):
    """
    Return sorted & deduplicated list of (Point, kp) along line.
    Uses fixed spacing + (optional) original vertices. Handles zero-length lines.
    """
    L = float(line.length or 0.0)
    if L == 0:
        p = line.interpolate(0.0)
        return [(p, 0.0)]

    s = _safe_step(step, L)

    pts = []
    d = 0.0
    while d < L - 1e-9:
        p = line.interpolate(d)
        pts.append((p, d))
        d += s
    # ensure end point
    pts.append((line.interpolate(L), L))

    if preserve_nodes and line.coords:
        for xy in line.coords:
            p = Point(xy)
            kp = line.project(p)
            pts.append((p, kp))

    # sort & dedup
    pts.sort(key=lambda t: t[1])
    dedup, seen = [], set()
    for p, kp in pts:
        key = (round(kp, 6), round(p.x, 6), round(p.y, 6))
        if key in seen:
            continue
        seen.add(key)
        dedup.append((p, kp))
    return dedup

def process_line_data(
    gis_path,
    group_field,
    fixed_distance,
    output_shp=True,
    tiff_path=None,
    preserve_nodes=True,
    preserve_attributes=True,
):
    # Read GIS
    gdf = gpd.read_file(gis_path)
    if gdf.crs is None:
        gdf.set_crs(epsg=3826, inplace=True)

    out_dir = os.path.join(os.path.dirname(gis_path) or os.getcwd(), "output")
    os.makedirs(out_dir, exist_ok=True)

    # Lazy import rasterio and open dataset
    ds = None
    use_geotiff = False
    to_raster_xy = _identity_xy
    if tiff_path:
        try:
            import rasterio
            ds = rasterio.open(tiff_path)
            use_geotiff = True
            to_raster_xy = _build_to_raster_transformer(gdf.crs, ds)
        except Exception as e:
            print(f"[WARN] Elevation disabled (raster not available or invalid): {e}")
            ds = None
            use_geotiff = False
            to_raster_xy = _identity_xy

    # lon/lat transformer for output
    to_wgs84 = Transformer.from_crs(gdf.crs, "EPSG:4326", always_xy=True)

    supported = {"LineString", "MultiLineString"}
    if group_field not in gdf.columns:
        gdf = gdf.assign(_group_fallback=gdf.index.astype(str))
        group_field = "_group_fallback"

    for group_val, sub in gdf.groupby(group_field):
        records = []
        shp_points = []

        for _, feat in sub.iterrows():
            geom = feat.geometry
            if geom is None or geom.is_empty:
                continue
            if geom.geom_type not in supported:
                continue

            lines = [geom] if isinstance(geom, LineString) else list(geom.geoms)
            for line in lines:
                samples = _sample_points_along_line(line, float(fixed_distance), preserve_nodes)
                prev_pt = None
                prev_elev = None
                total_3d = 0.0

                for pt, kp in samples:
                    lon, lat = to_wgs84.transform(pt.x, pt.y)
                    elev = _sample_elev_lazy(ds, pt.x, pt.y, to_raster_xy) if use_geotiff else None

                    if prev_pt is not None:
                        d2 = prev_pt.distance(pt)
                        if (prev_elev is not None) and (elev is not None):
                            dz = elev - prev_elev
                            seg3d = math.sqrt(d2 ** 2 + dz ** 2)
                            total_3d += seg3d
                        else:
                            seg3d = None
                        az = _azimuth(prev_pt, pt)
                    else:
                        d2 = 0.0
                        seg3d = None
                        az = None

                    rec = {
                        "Longitude": round(lon, 8),
                        "Latitude": round(lat, 8),
                        "Easting": pt.x,
                        "Northing": pt.y,
                        "Elevation": elev,
                        "Distance": round_by_distance(d2, fixed_distance),
                        "Length_3D": round_by_distance(seg3d, fixed_distance) if seg3d is not None else None,
                        "Azimuth": round(az, 3) if az is not None else None,
                        "KP": round_by_distance(kp, fixed_distance),
                        "Total_3D_Length": round_by_distance(total_3d, fixed_distance) if total_3d > 0 else None,
                        group_field: group_val,
                    }

                    if preserve_attributes:
                        for col in sub.columns:
                            if col != sub.geometry.name and col not in rec:
                                rec[col] = feat[col]

                    records.append(rec)
                    shp_points.append(pt)

                    prev_pt = pt
                    prev_elev = elev

        # filename suffix
        suffix = f"_{int(fixed_distance) if float(fixed_distance).is_integer() else fixed_distance}_node"
        safe_group = _sanitize_name(group_val)

        # CSV — always written
        csv_path = os.path.join(out_dir, f"{safe_group}{suffix}.csv")
        if records:
            pd.DataFrame(records).to_csv(csv_path, index=False)
        else:
            pd.DataFrame([], columns=[
                "Longitude","Latitude","Easting","Northing","Elevation","Distance","Length_3D",
                "Azimuth","KP","Total_3D_Length", group_field
            ]).to_csv(csv_path, index=False)

        # SHP — optional
        if output_shp and records:
            try:
                out_gdf = gpd.GeoDataFrame(records, geometry=shp_points, crs=gdf.crs)
                shp_path = os.path.join(out_dir, f"{safe_group}{suffix}.shp")
                out_gdf.to_file(shp_path)
            except Exception as e:
                print(f"[WARN] Failed to write SHP for group {group_val}: {e}")

    # close raster
    try:
        if ds is not None:
            ds.close()
    except Exception:
        pass

    return True
