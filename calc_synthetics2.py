import os
from pyrocko import gf
from pyrocko import trace
from pyrocko import io
import obspy
from pathlib import Path

engine = gf.LocalEngine(store_superdirs=['.'])
store_id = 'marmara_gfs'

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def read_station_metadata(stations_dir):
    """
    Read station metadata from XML files in the stations directory.
    
    Parameters:
    -----------
    stations_dir : str
        Path to the directory containing station XML files
        
    Returns:
    --------
    dict : Dictionary with station codes as keys and station information as values
           {
               'SLVT': {'lat': float, 'lon': float, 'network': str, 'station': str},
               'BOTS': {'lat': float, 'lon': float, 'network': str, 'station': str},
               ...
           }
    """
    stations_dict = {}
    stations_path = Path(stations_dir)
    
    # Find all XML files in the directory
    for xml_file in stations_path.glob('*.xml'):
        try:
            # Use obspy to read the station XML file
            inv = obspy.read_inventory(str(xml_file))
            
            # Extract station information
            for network in inv.networks:
                for station in network.stations:
                    station_code = station.code
                    stations_dict[station_code] = {
                        'lat': station.latitude,
                        'lon': station.longitude,
                        'network': network.code,
                        'station': station_code,
                        'elevation': station.elevation
                    }
                    print(f"Loaded station: {network.code}.{station_code} - Lat: {station.latitude:.4f}, Lon: {station.longitude:.4f}")
        except Exception as e:
            print(f"Warning: Could not read {xml_file}: {e}")
    
    return stations_dict


def create_targets_for_station(station_code, lat, lon, elevation, store_id, channels=['Z', 'N', 'E']):
    """
    Create GF targets for a single station across specified channels.
    
    Parameters:
    -----------
    station_code : str
        Station code (e.g., 'SLVT')
    lat : float
        Station latitude
    lon : float
        Station longitude
    store_id : str
        Green's function store ID
    channels : list
        List of channel codes (default: ['Z', 'N', 'E'])
        
    Returns:
    --------
    list : List of gf.Target objects
    """
    targets = []
    for channel in channels:
        target = gf.Target(
            lat=lat,
            lon=lon,
            elevation=elevation,
            store_id=store_id,
            codes=('', station_code, '', channel)
        )
        targets.append(target)
    return targets


def read_observed_waveforms(waveforms_dir, station_code, channels=None):
    """
    Read observed waveform data for a station from mseed files.
    Automatically detects available channel types (HH or HN).
    
    Parameters:
    -----------
    waveforms_dir : str
        Path to directory containing waveform files
    station_code : str
        Station code (e.g., 'SLVT')
    channels : list, optional
        List of channel codes to read. If None, auto-detects available channels.
        
    Returns:
    --------
    obspy.Stream : Stream containing the waveform traces
    """
    stream = obspy.Stream()
    waveforms_path = Path(waveforms_dir)
    
    # Auto-detect available channel types if not specified
    if channels is None:
        # Check which channel type is available for this station
        hhz_files = list(waveforms_path.glob(f"KO.{station_code}..HHZ*.mseed"))
        hnz_files = list(waveforms_path.glob(f"KO.{station_code}..HNZ*.mseed"))
        
        if hhz_files:
            channels = ['HHE', 'HHN', 'HHZ']
        elif hnz_files:
            channels = ['HNE', 'HNN', 'HNZ']
        else:
            print(f"Warning: Could not auto-detect channel type for {station_code}")
            return stream
    
    for channel in channels:
        # Try to find and read the mseed file
        try:
            # Use glob to find the mseed file with wildcard date/time
            mseed_files = list(waveforms_path.glob(f"KO.{station_code}..{channel}*.mseed"))
            if mseed_files:
                filename = str(mseed_files[0])
                trace_obj = obspy.read(filename)
                stream += trace_obj
            else:
                print(f"Warning: No mseed file found for {station_code} {channel}")
        except Exception as e:
            print(f"Warning: Could not read {station_code} {channel}: {e}")
    
    return stream


def process_waveforms(stream, detrend_methods=['demean', 'linear'], 
                      taper_percentage=0.05, lowpass_freq=0.2):
    """
    Process waveform data (detrending, tapering, filtering).
    
    Parameters:
    -----------
    stream : obspy.Stream
        Input waveform stream
    detrend_methods : list
        Detrending methods to apply
    taper_percentage : float
        Taper percentage (default: 0.05)
    lowpass_freq : float
        Lowpass filter frequency in Hz
        
    Returns:
    --------
    obspy.Stream : Processed waveform stream
    """
    processed = stream.copy()
    
    for method in detrend_methods:
        processed = processed.detrend(method)
    
    processed = processed.taper(max_percentage=taper_percentage, type='cosine')
    processed = processed.filter('lowpass', freq=lowpass_freq)
    
    return processed


# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Load station metadata from XML files
stations_dict = read_station_metadata('20250425/stations')
print(f"\nLoaded {len(stations_dict)} stations from metadata.\n")

# Source parameters
stf = gf.HalfSinusoidSTF(duration=5.0)

event_lat = 40.8558
event_lon = 28.4311
event_depth = 10000.0  # meters (adjusted to 10 km - max depth in GF store)
event_time = '2025-04-25 17:33:15'

# Convert event time to timestamp
from pyrocko.util import str_to_time
event_time_stamp = str_to_time(event_time)

source = gf.DCSource(
    lat=event_lat,
    lon=event_lon,
    depth=event_depth,
    strike=168.0,
    dip=47.0,
    rake=-37.0,
    magnitude=4.2,
    time=event_time_stamp,  # Use actual event time
    stf=stf
)

# Create targets for the stations we want
targets = []

# Automatically use all stations from the loaded metadata
# station_codes_to_use = list(stations_dict.keys())
station_codes_to_use = ['SLVT', 'BOTS', 'BGKT', 'KAVV']
print(f"Processing stations: {', '.join(station_codes_to_use)}\n")

for station_code in station_codes_to_use:
    if station_code in stations_dict:
        station_info = stations_dict[station_code]
        station_targets = create_targets_for_station(
            station_code=station_code,
            lat=station_info['lat'],
            lon=station_info['lon'],
            elevation=station_info['elevation'],
            store_id=store_id
        )
        targets.extend(station_targets)
        print(f"Created targets for station {station_code}")



print("Calculating synthetics...")
response = engine.process(source, targets)
synthetic_traces = response.pyrocko_traces()

print(f"Generated {len(synthetic_traces)} traces.")
for tr in synthetic_traces:
    print(f"Trace: {tr.nslc_id}")
    print(f"  Min Amp: {tr.ydata.min()}")
    print(f"  Max Amp: {tr.ydata.max()}")
    print(f"  Start time (tmin): {tr.tmin} seconds")
    print(f"  End time (tmax): {tr.tmax} seconds")
    print(f"  Duration: {tr.tmax - tr.tmin} seconds")
    print(f"  Sample rate: {1.0/tr.deltat} Hz")

output_filename = 'PYROCKO.syn.mseed'
io.save(synthetic_traces, output_filename)
print(f"Saved traces to {output_filename}\n")

synth_traces = obspy.read("PYROCKO.syn.mseed")

# Read and process observed waveforms for each station
print("Reading observed waveforms...")
observed_traces = {}
for station_code in station_codes_to_use:
    print(f"\nProcessing station: {station_code}")
    stream = read_observed_waveforms('20250425/waveforms', station_code)
    if len(stream) > 0:
        stream = process_waveforms(stream)
        observed_traces[station_code] = stream
        print(f"  Loaded and processed {len(stream)} traces")
    else:
        print(f"  No waveforms found for station {station_code}")



# Check timing information for observed traces
print("\n--- Observed Traces Timing ---")
for station_code, stream in observed_traces.items():
    for tr in stream:
        print(f"{station_code} {tr.stats.component}: {tr.stats.starttime} to {tr.stats.endtime}")
        print(f"  Duration: {tr.stats.endtime - tr.stats.starttime} seconds")

print("\n--- Synthetic Traces Timing ---")
for tr in synth_traces:
    print(f"Synthetic {tr.stats.station}.{tr.stats.component}: starts at {tr.stats.starttime}")
    print(f"  End time: {tr.stats.endtime}")


# Plot comparison for each station
import matplotlib.pyplot as plt

for station_code in station_codes_to_use:
    if station_code not in observed_traces:
        continue
    
    obs_stream = observed_traces[station_code]
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    components = ['Z', 'N', 'E']
    for i, comp in enumerate(components):
        # Get synthetic trace
        synth_trace = synth_traces.select(station=station_code, component=comp)
        if len(synth_trace) > 0:
            comp_data = synth_trace[0].data
            comp_data /= max(abs(comp_data))  # Normalize synthetic data
            comp_data *= 1
            comp_data = comp_data[175:]
            comp_times = synth_trace[0].times()
            comp_times = comp_times[:-175]
            axs[i].plot(comp_times, comp_data, label='Synthetic', color='blue')
        
        # Get observed trace
        obs_trace = obs_stream.select(component=comp)
        if len(obs_trace) > 0:
            obs_data = obs_trace[0].data
            obs_data /= max(abs(obs_data))  # Normalize observed data
            obs_times = obs_trace[0].times()
            obs_data = obs_data[:12000]
            obs_times = obs_times[:12000]
            axs[i].plot(obs_times, obs_data, label='Observed', color='red', alpha=0.7)
        
        axs[i].set_title(f'{station_code} Station - {comp} Component')
        axs[i].legend()
        axs[i].set_ylabel('Amplitude')
        axs[i].tick_params(labelbottom=False, labelleft=False)  # Remove numbers from axes

    plt.xlabel('Time (s)')
    plt.tight_layout()
    plt.show()
