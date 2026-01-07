import os
from obspy import read_inventory

stationxml_dir = "20250425/stations"
output_file = os.path.join(stationxml_dir, "STATION0")

def decimal_to_deg_min(value, is_lon=False):
    """Convert decimal degrees to (deg, min, hemisphere)."""
    if is_lon:
        hemi = "E" if value >= 0 else "W"
    else:
        hemi = "N" if value >= 0 else "S"
    abs_val = abs(value)
    deg = int(abs_val)
    minutes = (abs_val - deg) * 60
    return deg, minutes, hemi

def format_station_line(station, lat, lon, elev):
    """Return a properly formatted SEISAN station line."""
    lat_deg, lat_min, lat_hemi = decimal_to_deg_min(lat, is_lon=False)
    lon_deg, lon_min, lon_hemi = decimal_to_deg_min(lon, is_lon=True)
    elev = int(round(elev))

    # Format latitude/longitude with fixed widths
    lat_str = f"{lat_deg:02d}{lat_min:05.2f}{lat_hemi}"
    lon_str = f"{lon_deg:02d}{lon_min:05.2f}{lon_hemi}"

    # Adjust spacing based on station code length (strict SEISAN columns)
    if len(station) == 5:
        # 1 leading space for 5-character station name
        line = f" {station}{lat_str} {lon_str}{elev:4d} 0.00   0.00\n"
    elif len(station) == 4:
        # 2 leading spaces for 4-character station name
        line = f"  {station}{lat_str} {lon_str}{elev:4d} 0.00   0.00\n"
    elif len(station) == 3:
        # 3 leading spaces for 3-character station name
        line = f"   {station}{lat_str} {lon_str}{elev:4d} 0.00   0.00\n"
    else:
        # Default fallback (rare)
        line = f"  {station[:5]:5s}{lat_str} {lon_str}{elev:4d} 0.00   0.00\n"

    return line

# Collect XML files
xml_files = sorted(f for f in os.listdir(stationxml_dir) if f.endswith(".xml"))
lines = []

for xml_file in xml_files:
    inv = read_inventory(os.path.join(stationxml_dir, xml_file))
    sta = inv[0][0].code
    lat = inv[0][0].latitude
    lon = inv[0][0].longitude
    elev = inv[0][0].elevation
    lines.append(format_station_line(sta, lat, lon, elev))

# Write STATION0
with open(output_file, "w") as f:
    f.writelines(lines)

print(f"✅ STATION0 created at {output_file}")
for l in lines:
    print(l.strip())

