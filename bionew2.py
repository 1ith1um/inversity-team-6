
import geopandas as gpd
import h3
import networkx as nx
from shapely.geometry import Polygon, Point, LineString
import matplotlib.pyplot as plt

# Provided coordinates for Offshore Wind Farms and Coastal Landing Sites
shore_points = gpd.GeoDataFrame(geometry=[
    Point(50.99204, -2.13736),  
    Point(50.92476, -1.74743),
    Point(50.19311, -4.24092),
    Point(50.00774, -1.61027),
    Point(50.99334, -2.77223)
])

wind_farm_points = gpd.GeoDataFrame(geometry=[
    Point(51.34001, -1.94266),  
    Point(51.73717, -3.53176),
    Point(51.05769, -0.57098),
    Point(51.97663, -2.69583),
    Point(51.20804, -3.12798)
])

# Example biodiverse areas
biodiverse_points = gpd.GeoDataFrame(geometry=[
    Point(51.5, -2.5),
    Point(52.0, -3.0)
])

# Generate a Hexagonal Grid
def create_hex_grid(bbox, resolution):
    hexes = h3.polyfill({
        'type': 'Polygon',
        'coordinates': [[
            [bbox[0], bbox[1]],
            [bbox[2], bbox[1]],
            [bbox[2], bbox[3]],
            [bbox[0], bbox[3]],
            [bbox[0], bbox[1]]
        ]]
    }, resolution)
    
    hexagons = []
    for hex in hexes:
        hex_boundary = h3.h3_to_geo_boundary(hex, geo_json=True)
        hexagons.append(Polygon(hex_boundary))
        
    return gpd.GeoDataFrame(geometry=hexagons)

bbox = [-5, 50, 0, 54]  # Adjusted bounding box to cover the points provided
hex_grid = create_hex_grid(bbox, resolution=5)

def snap_points_to_hex_centers(points, resolution):
    hex_centers = []
    for point in points.geometry:
        hex_id = h3.geo_to_h3(point.y, point.x, resolution)
        hex_center = h3.h3_to_geo(hex_id)
        hex_centers.append(Point(hex_center[::-1]))
    return gpd.GeoDataFrame(geometry=hex_centers)

snapped_shore_points = snap_points_to_hex_centers(shore_points, resolution=5)
snapped_wind_farm_points = snap_points_to_hex_centers(wind_farm_points, resolution=5)
snapped_biodiverse_points = snap_points_to_hex_centers(biodiverse_points, resolution=5)

def create_graph_from_hex_grid(hex_grid, biodiverse_points, buffer_radius=2):
    G = nx.Graph()
    for idx, hex in hex_grid.iterrows():
        G.add_node(idx, geometry=hex.geometry)
    
    biodiverse_hexes = set()
    for point in biodiverse_points.geometry:
        hex_id = h3.geo_to_h3(point.y, point.x, resolution=5)
        biodiverse_hexes.add(hex_id)
        neighbors = h3.k_ring(hex_id, buffer_radius)
        biodiverse_hexes.update(neighbors)

    for idx, hex in hex_grid.iterrows():
        hex_id = h3.geo_to_h3(hex.geometry.centroid.y, hex.geometry.centroid.x, resolution=5)
        if hex_id in biodiverse_hexes:
            continue
        neighbors = h3.k_ring(hex_id, 1)
        for neighbor in neighbors:
            if neighbor in biodiverse_hexes:
                continue
            neighbor_center = Point(h3.h3_to_geo(neighbor)[::-1])
            neighbor_hex = hex_grid[hex_grid.contains(neighbor_center)]
            if not neighbor_hex.empty:
                neighbor_idx = neighbor_hex.index[0]
                if G.has_node(neighbor_idx):
                    G.add_edge(idx, neighbor_idx)
    
    return G

G = create_graph_from_hex_grid(hex_grid, snapped_biodiverse_points)

def find_closest_hex(point, hex_grid):
    min_distance = float('inf')
    closest_hex = None
    for idx, hex in hex_grid.iterrows():
        distance = point.distance(hex.geometry.centroid)
        if distance < min_distance:
            min_distance = distance
            closest_hex = idx
    return closest_hex

def plot_hex_grid_with_points_and_paths(hex_grid, offshore_wind_farms, coastal_landing_sites, biodiverse_points, paths):
    fig, ax = plt.subplots()
    hex_grid.plot(ax=ax, color='white', edgecolor='black')
    offshore_wind_farms.plot(ax=ax, color='blue', markersize=5)
    coastal_landing_sites.plot(ax=ax, color='orange', markersize=5)
    biodiverse_points.plot(ax=ax, color='red', markersize=10)

    for path in paths:
        if len(path) > 1:
            path_coords = [hex_grid.loc[idx].geometry.centroid for idx in path]
            path_line = LineString(path_coords)
            gpd.GeoSeries([path_line]).plot(ax=ax, color='green')
        else:
            print("Path is too short to display as a line.")
    
    plt.show()

paths = []
for wind_farm_point in snapped_wind_farm_points.geometry:
    best_path = None
    best_length = float('inf')
    best_target = None

    wind_farm_hex = find_closest_hex(wind_farm_point, hex_grid)
    if wind_farm_hex is None:
        print(f"Skipping wind farm {wind_farm_point} - no corresponding hex found.")
        continue

    for coastal_site in snapped_shore_points.geometry:
        coastal_site_hex = find_closest_hex(coastal_site, hex_grid)
        if coastal_site_hex is None:
            print(f"Skipping coastal site {coastal_site} - no corresponding hex found.")
            continue

        try:
            path = nx.shortest_path(G, source=wind_farm_hex, target=coastal_site_hex)
            if len(path) < best_length:
                best_path = path
                best_length = len(path)
                best_target = coastal_site
        except nx.NetworkXNoPath:
            continue

    if best_path:
        paths.append(best_path)
        print(f"Path from wind farm {wind_farm_point} to coastal site {best_target} established.")
    else:
        print(f"No path found from wind farm {wind_farm_point}. Attempting alternatives.")

plot_hex_grid_with_points_and_paths(hex_grid, snapped_wind_farm_points, snapped_shore_points, snapped_biodiverse_points, paths)
