import os
import math
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString
from pyproj import Transformer

# NOTE: rasterio is imported lazily only if a GeoTIFF path is provided.

def round_by_distance(value: float, step: float) -> float:
    if value is None:
        return None
    try:
        s = float(step)
    except Exception:
        s = 0.0
    if s < 1:
        return round(value, 2)
    elif s < 10:
        return round(value, 1)
    else:
        return round(value)

def _azimuth(p_prev: Point, p_cur: Point):
    if p_prev is None:
        return None
    dx = p_cur.x - p_prev.x
    dy = p_cur.y - p_prev.y
    ang = math.degrees(math.atan2(dx, dy))
    return (ang + 360.0) % 360.0

def _within_bounds(ds, x, y) -> bool:
    b = ds.bounds  # left, bottom, right, top
    return (b.left <= x <= b.right) and (b.bottom <= y <= b.top)

def _sample_elev_lazy(ds, x, y):
    try:
        if ds is None:
            return None
        if not _within_bounds(ds, x, y):
            return None
        row, col = ds.index(x, y)
        if row < 0 or col < 0 or row >= ds.height or col >= ds.width:
            return None
        val = ds.read(1)[row, col]
        # 處理常見 NoData
        if val is None:
            return None
        if isinstance(val, float) and (math.isnan(val) or val >= 3.4e38):
            return None
        return float(val)
    except Exception:
        return None

def _safe_step(step: float, L: float) -> float:
    """保證步長安全：非數值/負數 → 使用 L（只取端點）；L=0 → 回傳 0 表示只有一個點可取。"""
    try:
        s = float(step)
    except Exception:
        s = 0.0
    if s <= 0:
        s = L  # 只取端點
    return s

def _sample_points_along_line(line: LineString, step: float, preserve_nodes: bool):
    """
    回傳 [(Point, kp)] 排序於線上，避免 0 長度/非法步長造成錯誤。
    規則：固定距離 +（可選）原始節點 + 去重。
    """
    L = float(line.length or 0.0)
    # 處理退化線 (L == 0)：只回傳單點
    if L == 0:
        p = line.interpolate(0.0)
        pts = [(p, 0.0)]
        if preserve_nodes:
            # 若原本就只有一個座標，這裡也只會有 1 點
            pass
        return pts

    s = _safe_step(step, L)

    pts = []
    # 產生等距 KP；使用 while 避免 np.arange 的浮點坑
    d = 0.0
    while d < L - 1e-9:
        p = line.interpolate(d)
        pts.append((p, d))
        d += s
    # 確保加入末端點
    pts.append((line.interpolate(L), L))

    # 加入原始節點（可選）
    if preserve_nodes and line.coords:
        for xy in line.coords:
            p = Point(xy)
            kp = line.project(p)
            pts.append((p, kp))

    # 排序 + 去重
    pts.sort(key=lambda t: t[1])
    dedup = []
    seen = set()
    for p, kp in pts:
        key = (round(kp, 6), round(p.x, 6), round(p.y, 6))
        if key in seen:
            continue
        seen.add(key)
        dedup.append((p, kp))
    return dedup

def _sanitize_name(v):
    s = str(v).strip()
    # Windows 檔名不允許字元清掉
    for ch in r'\/:*?"<>|':
        s = s.replace(ch, "_")
    return s or "group"

def process_line_data(
    gis_path,
    group_field,
    fixed_distance,
    output_shp=True,
    tiff_path=None,
    preserve_nodes=True,
    preserve_attributes=True,
):
    # 讀檔
    gdf = gpd.read_file(gis_path)
    if gdf.crs is None:
        gdf.set_crs(epsg=3826, inplace=True)

    out_dir = os.path.join(os.path.dirname(gis_path) or os.getcwd(), "output")
    os.makedirs(out_dir, exist_ok=True)

    # 懶載入 rasterio（沒 GeoTIFF 也能跑）
    ds = None
    use_geotiff = False
    if tiff_path:
        try:
            import rasterio  # lazy import
            ds = rasterio.open(tiff_path)
            use_geotiff = True
        except Exception as e:
            print(f"[WARN] Elevation disabled: {e}")
            ds = None
            use_geotiff = False

    transformer = Transformer.from_crs(gdf.crs, "EPSG:4326", always_xy=True)

    supported = {"LineString", "MultiLineString"}
    if group_field not in gdf.columns:
        gdf = gdf.assign(_group_fallback=gdf.index.astype(str))
        group_field = "_group_fallback"

    # 逐群組處理
    for group_val, sub in gdf.groupby(group_field):
        records = []
        shp_points = []

        for _, feat in sub.iterrows():
            geom = feat.geometry
            if geom is None or geom.is_empty:
                continue
            if geom.geom_type not in supported:
                # 跳過非線
                continue

            lines = [geom] if isinstance(geom, LineString) else list(geom.geoms)
            for line in lines:
                samples = _sample_points_along_line(line, float(fixed_distance), preserve_nodes)
                prev_point = None
                prev_elev = None
                total_3d = 0.0

                for i, (pt, kp) in enumerate(samples):
                    lon, lat = transformer.transform(pt.x, pt.y)
                    elev = _sample_elev_lazy(ds, pt.x, pt.y) if use_geotiff else None

                    # 距離（2D/3D）
                    if prev_point is not None:
                        d2 = prev_point.distance(pt)
                        if (prev_elev is not None) and (elev is not None):
                            dz = elev - prev_elev
                            seg3d = math.sqrt(d2 ** 2 + dz ** 2)
                            total_3d += seg3d
                        else:
                            seg3d = None
                        az = _azimuth(prev_point, pt)
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

                    prev_point = pt
                    prev_elev = elev

        # 檔名：{group}_{distance}_node
        suffix = f"_{int(fixed_distance) if float(fixed_distance).is_integer() else fixed_distance}_node"
        safe_group = _sanitize_name(group_val)

        # 一定輸出 CSV（即使沒有點，也給出 header 方便排錯）
        csv_path = os.path.join(out_dir, f"{safe_group}{suffix}.csv")
        if records:
            pd.DataFrame(records).to_csv(csv_path, index=False)
        else:
            pd.DataFrame([], columns=[
                "Longitude","Latitude","Easting","Northing","Elevation","Distance","Length_3D",
                "Azimuth","KP","Total_3D_Length", group_field
            ]).to_csv(csv_path, index=False)

        # Shapefile（可選）
        if output_shp and records:
            try:
                out_gdf = gpd.GeoDataFrame(records, geometry=shp_points, crs=gdf.crs)
                shp_path = os.path.join(out_dir, f"{safe_group}{suffix}.shp")
                out_gdf.to_file(shp_path)
            except Exception as e:
                print(f"[WARN] Failed to write SHP for group {group_val}: {e}")

    # 關閉 raster
    try:
        if ds is not None:
            ds.close()
    except Exception:
        pass

    return True
