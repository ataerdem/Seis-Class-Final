import os
import re
import pandas as pd
from obspy import read, UTCDateTime, read_inventory
from obspy.core.event import Catalog, Event, Pick, Origin
from obspy.core.event.base import WaveformStreamID

# Parameters -----------------------------
stationxml_dir = "20250425/stations"
# eventdir = "."
marker_file = "20250425/markers.txt"
hypout_file = "20250425/hyp.in"
# ---------------------------------------

def get_station_location(network, station, stationxml_dir):
    stationxml_file = os.path.join(stationxml_dir, f"{network}.{station}.xml")
    if not os.path.exists(stationxml_file):
        return None, None
    inv = read_inventory(stationxml_file)
    latitude = inv[0][0].latitude
    longitude = inv[0][0].longitude
    return latitude, longitude


def parse_markers(marker_file):
    """
    Parse the markers file and return a DataFrame
    with UTCDateTime objects and phase, network, station info.
    """
    pattern = re.compile(
        r"phase:\s+([\d\-:\. ]+)\s+\d+\s+([A-Z0-9]+)\.([A-Z0-9]+)\.\.[A-Z0-9]+\s+None\s+None\s+None\s+([PS])"
    )

    pick_times, networks, stations, phases = [], [], [], []

    with open(marker_file, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                pick_time_str, network, station, phase = match.groups()
                try:
                    pick_time = UTCDateTime(pick_time_str.strip())
                    pick_times.append(pick_time)
                    networks.append(network.strip())
                    stations.append(station.strip())
                    phases.append(phase.strip())
                except Exception as e:
                    print(f"Could not parse line: {line.strip()}")
                    continue

    df = pd.DataFrame({
        "pick_time": pick_times,
        "network": networks,
        "station": stations,
        "phase": phases
    })
    return df


def buildevent(marker_file, stationxml_dir):
    df = parse_markers(marker_file)
    event = Event()
    earliest_P_time = None
    lat_orig, lon_orig = None, None

    for _, row in df.iterrows():
        pick_time = row["pick_time"]
        network = row["network"]
        station = row["station"]
        phase = row["phase"]

        WaveID = WaveformStreamID(network_code=network, station_code=station)
        pick = Pick(time=pick_time, phase_hint=phase, waveform_id=WaveID, evaluation_mode="manual")

        # Weight based on phase type
        if phase == "P":
            pick.extra = {'nordic_pick_weight': {'value': 0, 'namespace': 'https://seis.geus.net/software/seisan/node239.html'}}
            if earliest_P_time is None or pick_time < earliest_P_time:
                earliest_P_time = pick_time
                lat_orig, lon_orig = get_station_location(network, station, stationxml_dir)
        else:
            pick.extra = {'nordic_pick_weight': {'value': 2, 'namespace': 'https://seis.geus.net/software/seisan/node239.html'}}

        event.picks.append(pick)

    if earliest_P_time and lat_orig is not None:
        origin = Origin()
        origin.time = earliest_P_time
        origin.latitude = lat_orig
        origin.longitude = lon_orig
        origin.depth = 0
        origin.method_id = "Earliest P phase arrival"
        event.origins.append(origin)
        print(f"Earliest P arrival: {earliest_P_time} at station near ({lat_orig}, {lon_orig})")
    else:
        print("No P phase found or station location missing")

    return event


# --- Run the event builder ---
catalog = Catalog()
event = buildevent(marker_file, stationxml_dir)
catalog.append(event)

# Write the output in NORDIC format
catalog.write(hypout_file, format="NORDIC")
print(f"Written to {hypout_file}")

