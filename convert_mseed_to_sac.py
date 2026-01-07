#!/usr/bin/env python
"""
Convert MiniSEED waveforms to SAC format for magnitude calculation
"""
import os
import glob
from obspy import read
from obspy.io.sac import SACTrace

# Parameters
waveform_dir = "20250425/waveforms"
station_xml_dir = "20250425/stations"

def convert_mseed_to_sac(mseed_file, output_dir=None):
    """Convert a MiniSEED file to SAC format with proper headers"""
    try:
        # Read the MiniSEED file
        st = read(mseed_file)
        
        for tr in st:
            # Create output filename
            if output_dir is None:
                output_dir = os.path.dirname(mseed_file)
            
            # Generate SAC filename: NET.STA.LOC.CHA.SAC
            sac_filename = f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel}.SAC"
            sac_path = os.path.join(output_dir, sac_filename)
            
            # Write to SAC format
            tr.write(sac_path, format='SAC')
            print(f"Converted: {os.path.basename(mseed_file)} -> {sac_filename}")
            
    except Exception as e:
        print(f"Error converting {mseed_file}: {e}")

def main():
    # Find all MiniSEED files
    mseed_files = glob.glob(os.path.join(waveform_dir, "*.mseed"))
    
    if not mseed_files:
        print(f"No .mseed files found in {waveform_dir}")
        return
    
    print(f"Found {len(mseed_files)} MiniSEED files to convert")
    print("-" * 60)
    
    # Convert each file
    for mseed_file in sorted(mseed_files):
        convert_mseed_to_sac(mseed_file)
    
    print("-" * 60)
    print(f"Conversion complete! {len(mseed_files)} files converted to SAC format")

if __name__ == "__main__":
    main()
