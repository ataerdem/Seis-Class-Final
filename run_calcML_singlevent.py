import os
import sys
import numpy as np
from obspy import read, read_events, UTCDateTime
from obspy.core.event import Magnitude
from obspy.core.event.base import WaveformStreamID
from obspy import Catalog

# ---------------- PARAMETERS ----------------
inputcatalogfile = "20250425/hyp.out"
outputcatalogfile = 'hypM.out'
hyp_in_file = "20250425/hyp.in"  # File with pick times

event_dir = '20250425/waveforms'
xmldir = '20250425/stations'
# path_to_calcML = '/Users/aliozgunkonca/Documents/Teaching/seismicdata_tutorial/magnitude/'
SNR_threshold = 1.3

# List of stations to exclude from magnitude calculation
# Add station names here (e.g., ['BGKT', 'SILV', 'MRMT'])
excluded_stations = []
# --------------------------------------------

# sys.path.append(path_to_calcML)
import calcMLsingle as cML  # assumes this has a `calc_mL()` function

def parse_hyp_in_picks(filename):
    """Parse hyp.in file to extract P-wave pick times for each station"""
    station_picks = {}
    with open(filename, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith(' 2025') or 'Action:' in line or 'STAT SP' in line:
                continue
            if not line.strip():
                continue
            
            parts = line.split()
            if len(parts) >= 4 and parts[0].isalpha():
                station = parts[0]
                phase = parts[1]
                
                # Only process P-wave picks
                if phase == 'P':
                    try:
                        # Parse time: HHMM SS.SSS format (e.g., 173321.041)
                        time_str = parts[3]
                        hhmm = int(time_str[:4])
                        seconds = float(time_str[4:])
                        hour = hhmm // 100
                        minute = hhmm % 100
                        
                        # Store pick time info
                        station_picks[station] = {
                            'hour': hour,
                            'minute': minute,
                            'second': seconds,
                            'time_string': time_str
                        }
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing pick for station {station}: {e}")
                        continue
    
    return station_picks

inputcatalogpath= os.path.join(inputcatalogfile)
outputcatalogpath= os.path.join(event_dir,outputcatalogfile)
cat=read_events(inputcatalogpath)
event=cat[0]
origin_time = event.origins[0].time

# Load pick times from hyp.in
station_picks = parse_hyp_in_picks(hyp_in_file)
print(f"\nLoaded {len(station_picks)} P-wave picks from {hyp_in_file}")
for station, pick_info in station_picks.items():
    print(f"  {station}: {pick_info['time_string']}")


event_id = os.path.basename(event_dir)
eventmagfile = os.path.join(event_dir, 'magnitude.txt')
hypM_event = os.path.join(event_dir,outputcatalogfile)

weightsum, magsum = 0, 0
magnitudes, distances, stnames, weights = [], [], [], []

with open(eventmagfile, 'w') as f:
    f.write("eventID stname channel distance mag weight peak2peak tdif\n")
    for filename in sorted(os.listdir(event_dir)):
        if filename.endswith(("HH1.SAC", "HHE.SAC", "BHE.SAC", "EHE.SAC", 
                              "HH2.SAC", "HHN.SAC", "BHN.SAC", "EHN.SAC" )):
            sac_file_path = os.path.join(event_dir, filename)
            
            # Extract station name from filename
            # Assumes format like: KO.BGKT..HHE.SAC or similar
            filename_base = os.path.basename(sac_file_path)
            parts = filename_base.split('.')
            if len(parts) >= 2:
                station_name = parts[1]
            else:
                print(f"Cannot parse station name from {filename}")
                continue
            
            # Only process stations that have P-wave picks in hyp.in
            if station_name not in station_picks:
                print(f"Skipping {filename} - no P-wave pick found in hyp.in for {station_name}")
                continue
            
            # Skip excluded stations
            if station_name in excluded_stations:
                print(f"Skipping {filename} - {station_name} is in excluded stations list")
                continue
            
            print(f"working on file  {sac_file_path} (station {station_name} with pick)")
            try:
                st = read(sac_file_path)
                for tr in st:
                    if tr.stats.starttime > tr.stats.endtime:
                        print(f"Bad time in {sac_file_path}. Skipping.")
                        continue
            except Exception as e:
                print(f"Error reading {sac_file_path}: {e}")
                continue

            # Create pick time as UTCDateTime for this station
            pick_info = station_picks[station_name]
            pick_time = UTCDateTime(
                year=origin_time.year,
                month=origin_time.month,
                day=origin_time.day,
                hour=pick_info['hour'],
                minute=pick_info['minute'],
                second=int(pick_info['second']),
                microsecond=int((pick_info['second'] % 1) * 1e6)
            )
            
            # Write the pick time and other required info into the SAC header
            try:
                from obspy import read_inventory
                from obspy.geodetics import gps2dist_azimuth
                import glob as glob_module
                
                st_copy = st.copy()
                for tr in st_copy:
                    # Calculate time difference from trace start
                    time_diff = pick_time - tr.stats.starttime
                    
                    # Initialize SAC headers if not present
                    if not hasattr(tr.stats, 'sac') or tr.stats.sac is None:
                        tr.stats.sac = {}
                    
                    # Set P arrival time
                    tr.stats.sac['a'] = time_diff  # P arrival time relative to reference
                    tr.stats.sac['ka'] = 'P'  # Label for 'a' header
                    tr.stats.sac['t0'] = time_diff  # Also set t0 as some codes use this
                    tr.stats.sac['kt0'] = 'P'
                    
                    # Set event coordinates
                    tr.stats.sac['evla'] = event.origins[0].latitude
                    tr.stats.sac['evlo'] = event.origins[0].longitude
                    tr.stats.sac['evdp'] = event.origins[0].depth / 1000.0  # convert to km
                    
                    # Get station coordinates from XML
                    station_xml_pattern = os.path.join(xmldir, f"*.{station_name}.xml")
                    station_xml_files = glob_module.glob(station_xml_pattern)
                    
                    if station_xml_files:
                        inv = read_inventory(station_xml_files[0])
                        station_obj = inv[0][0]  # First network, first station
                        
                        tr.stats.sac['stla'] = station_obj.latitude
                        tr.stats.sac['stlo'] = station_obj.longitude
                        tr.stats.sac['stel'] = station_obj.elevation
                        
                        # Calculate epicentral distance
                        dist_m, az, baz = gps2dist_azimuth(
                            event.origins[0].latitude, event.origins[0].longitude,
                            station_obj.latitude, station_obj.longitude
                        )
                        tr.stats.sac['dist'] = dist_m / 1000.0  # convert to km
                        tr.stats.sac['az'] = az
                        tr.stats.sac['baz'] = baz
                    else:
                        print(f"  Warning: No station XML found for {station_name}")
                    
                # Write back to SAC file with updated headers
                st_copy.write(sac_file_path, format='SAC')
                print(f"  Added P-wave pick and metadata to SAC header")
            except Exception as e:
                print(f"Warning: Could not write headers to SAC: {e}")
                import traceback
                traceback.print_exc()
            
            try:
                # Try passing pick_time to calc_mL if it accepts it
                stname, channel, distance, mag, weight, peak2peak, tdif = cML.calc_mL(
                    sac_file_path, xmldir, SNR_threshold, origin_time, pick_time=pick_time
                )
            except TypeError:
                # If calc_mL doesn't accept pick_time parameter, try without it
                try:
                    stname, channel, distance, mag, weight, peak2peak, tdif = cML.calc_mL(sac_file_path, xmldir, SNR_threshold,origin_time)
                except Exception as e:
                    print(f"calc_mL failed on {sac_file_path}: {e}")
                    continue
            except Exception as e:
                print(f"calc_mL failed on {sac_file_path}: {e}")
                continue
            except Exception as e:
                print(f"calc_mL failed on {sac_file_path}: {e}")
                continue

            if weight > 1 and not np.isnan(mag):
                magnitudes.append(mag)
                distances.append(distance)
                stnames.append(stname)
                weights.append(weight)
                weightsum += weight
                magsum += weight * mag

            f.write(f"{event_id} {stname} {channel} {distance} {mag} {weight} {peak2peak} {tdif}\n")

    # Compute ML
    if weightsum > 1:
        mL_mag = magsum / weightsum
    else:
        mL_mag = -2

    mL_std = np.std(magnitudes) if len(magnitudes) > 2 else 0.25
    f.write(f"{mL_mag} {mL_std}\n")

    if len(magnitudes) > 3:
        threshold = mL_std
        filtered = [(mag, dist, stn, w) for mag, dist, stn, w in zip(magnitudes, distances, stnames, weights)
                    if abs(mag - mL_mag) <= threshold]
        if filtered:
            f_mags, f_weights = zip(*[(m, w) for m, _, _, w in filtered])
            total_weight = np.sum(f_weights)
            if total_weight > 0:
                norm_weights = [w / total_weight for w in f_weights]
                filtered_mag = np.sum(np.array(f_mags) * np.array(norm_weights))
                mL_mag = filtered_mag
                f.write(f"filteredML= {filtered_mag}\n")

    if mL_mag > 3:
        threshold_dist = 40 if mL_mag > 3.5 else 30
        new_mags = [m for d, m in zip(distances, magnitudes) if d >= threshold_dist and 3.2 < m < 7]
        if len(new_mags) > 1:
            avg = np.mean(new_mags)
            if avg > mL_mag:
                mL_mag = avg
                f.write(f"magnitude after second iteration: {mL_mag}\n")

        # Save magnitude to event if valid
if mL_mag > 0:
    magnitude = Magnitude()
    magnitude.mag = mL_mag
    magnitude.magnitude_type = 'ML'
    magnitude.mag_errors = mL_std
    event.magnitudes.append(magnitude)
    event.magnitudes.append(magnitude)
    event.magnitudes.append(magnitude)
else:
    print(f"ML not valid for event {event_id} (mL_mag = {mL_mag})")
eventcat = Catalog() 
eventcat.append(event)
eventcat.write(hypM_event, format='NORDIC')


print("\nDone.")
