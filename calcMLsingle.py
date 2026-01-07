import math
from obspy.signal.invsim import estimate_magnitude
import obspy
from obspy import read, read_inventory
import numpy as np
import os
import sys

def calculate_snr(sac_file_path, datawin):
    # Read SAC file
    print(f"Working on calculate_snr {sac_file_path}")
    try:
        st = read(sac_file_path)
    except Exception as e:
        print(f"Error reading SAC file {sac_file_path}: {str(e)}")
        return 0
    
    st.taper(max_percentage=0.05, type='cosine')
    #st.filter("bandpass", freqmin=0.2,freqmax=5,corners=6, zerophase=True)
    st.filter("bandpass", freqmin=1,freqmax=8,corners=6, zerophase=True)
    tr = st[0]

    # Get P and/or S arrival times
    if 't0' in tr.stats.sac:  # Assuming 't0' corresponds to P arrival time
        tarr = tr.stats.sac.t0
    elif 't2' in tr.stats.sac:  # Assuming 't2' corresponds to S arrival time
        dist = tr.stats.sac.dist
        tarr = tr.stats.sac.t2 - dist / 8  # Rough calculation for P arrival time from S arrival time
    else:
        print(f'Warning: P or S arrival time not found for {sac_file_path}. Exiting...')
        return 0

    print(f"In calculate_snr tarr={tarr}")
    noise_starttime = tr.stats.starttime
    noise_endtime = noise_starttime + datawin
    data_starttime = tr.stats.starttime + tarr
    data_endtime = tr.stats.starttime + tarr + datawin

    # Check if the data window goes beyond trace's endtime
    if data_endtime > tr.stats.endtime:
        print(f"Warning: your seismogram's end time is smaller than data window.")
        return 0  # Returning zero for SNR when the data window is out of bounds

    # Get the data within the time window
    noise_tr = tr.copy()
    data_tr = tr.copy()
    noise_window = noise_tr.trim(starttime=noise_starttime, endtime=noise_endtime)
    data_window = data_tr.trim(starttime=data_starttime, endtime=data_endtime)

    # Calculate SNR
    noise_amp = np.mean(np.abs(noise_window.data - np.mean(noise_window.data)))
    signal_amp = np.mean(np.abs(data_window.data - np.mean(data_window.data)))
    
    if noise_amp == 0:  # Handle division by zero if noise_amp is zero
        return 0

    snr = signal_amp / noise_amp
    return snr

def read_event_data(sac_file_path,origin_time):
    try:
        stream = read(sac_file_path)
        if len(stream) == 0:
            print(f"No data found in the SAC file {sac_file_path}.")
            return None, None, None, None, None, None, None

        stream.detrend('demean')
        stream.detrend('linear')
        stream.taper(max_percentage=0.05, type='cosine')
        stream.filter("bandpass", freqmin=0.25, freqmax=5)
        trace = stream[0]
        tr_endtime = trace.stats.endtime
        distance = trace.stats.sac.dist
        depth = trace.stats.sac.evdp
        data_starttime = origin_time + distance / 6
        calc_endtime = data_starttime + distance / 8 + 15
        data_endtime = min(tr_endtime, calc_endtime)

        trace_dat = trace.trim(starttime=data_starttime, endtime=data_endtime)
        time_values = trace_dat.times()
        amplitude_values = trace_dat.data

        channel = trace.stats.channel
        stname = trace.stats.station
        network = trace.stats.network

        return time_values, amplitude_values, network, stname, channel, distance, depth
    except Exception as e:
        print(f"An error occurred while reading the SAC file: {str(e)}")
        return None, None, None, None, None, None, None

def find_max_min_amplitude(time_values, amplitude_values):
    if time_values is None or amplitude_values is None:
        return None, None, None
    max_amplitude_index = np.argmax(amplitude_values)
    min_amplitude_index = np.argmin(amplitude_values)

    max_amplitude = amplitude_values[max_amplitude_index]
    min_amplitude = amplitude_values[min_amplitude_index]

    max_amplitude_time = time_values[max_amplitude_index]
    min_amplitude_time = time_values[min_amplitude_index]

    return max_amplitude, min_amplitude, abs(max_amplitude_time - min_amplitude_time)

def extract_response(xmldir, netwk, station_code, channel_code):
    try:
        # Read the StationXML inventory file
        xmlname = f"{netwk}.{station_code}.xml"
        xmlfile = os.path.join(xmldir, xmlname)
        inventory = read_inventory(xmlfile)

        for network in inventory:
            for station in network:
                if station.code == station_code:
                    for channel in station:
                        if channel.code == channel_code:
                            response = channel.response
                            return response
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

def calc_mL(sac_file_path, xmldir, SNR_threshold,origin_time):
    datawin = 5  # Window length for noise and data time windows for SNR
    # Calculate signal-to-noise ratio
    SNR = calculate_snr(sac_file_path, datawin)

    if SNR < SNR_threshold:
        weight = 0
        magnitude = 0
        time_values, amplitude_values, network, stname, channel, distance, depth = read_event_data(sac_file_path,origin_time)
        peaktopeak = 0
        time_difference = 1
    else:
        # Read SAC file data
        time_values, amplitude_values, network, stname, channel, distance, depth = read_event_data(sac_file_path,origin_time)
        if time_values is None or amplitude_values is None:
            print("No valid event data found. Returning default values.")
            return None, None, None, None, None, None, None

        # Extract station response
        response = extract_response(xmldir, network, stname, channel)
        if response is None:
            print(f"Error: No response data found for {network}.{stname}.{channel}. Exiting...")
            return None, None, None, None, None, None, None
        
        max_amplitude, min_amplitude, time_difference = find_max_min_amplitude(time_values, amplitude_values)
        
        if abs(max_amplitude) < 0.1 or abs(min_amplitude) < 0.1:
            magnitude = 0
            peaktopeak = 0
            time_difference = 1
        else:
            peaktopeak = abs(max_amplitude - min_amplitude)
            h_dist = math.sqrt(distance ** 2 + depth ** 2)
            print(f"{stname}, {channel}, P2P {peaktopeak}, tdif={time_difference}, h_dist={h_dist}")
            magnitude = estimate_magnitude(response, peaktopeak, time_difference, h_dist)

        # Set the weight based on SNR
        weight = 2 if SNR > 2 else SNR

        # If the channel is a specific type, set weight to 0
        if channel in ['HHZ', 'HNZ', 'BHZ', 'SHZ', 'EHZ', 'HH3']:
            weight = 0

    return stname, channel, distance, magnitude, weight, peaktopeak, time_difference
