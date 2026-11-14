import h5py

def patch_header_time(fn):
    with h5py.File(fn,'r+') as f:
        hdr = f['Header'].attrs
        print("Old Time:", hdr['Time'])
        hdr['Time'] = 0.0
        hdr['Redshift'] = 0.0
        hdr['HubbleParam'] = 1.0
        print("Patched Time to:", hdr['Time'])



#patch_header_time('snapshot_038_cutout_split.hdf5')



with h5py.File("snapshot_038_cutout_split.hdf5", "r") as f:
    header = f["Header"].attrs
    print("Header attributes in")
    for key, val in header.items():
        print(f"{key}: {val}")




import numpy as np

def inspect_snapshot(fn):
    with h5py.File(fn,'r') as f:
        hdr = dict(f['Header'].attrs)
        print("Header['Time'] =", hdr.get('Time'))
        print("Header['Redshift'] =", hdr.get('Redshift'))
        print("Header['HubbleParam'] =", hdr.get('HubbleParam'))
        print("NumPart_ThisFile =", hdr.get('NumPart_ThisFile'))
        for pt in ['PartType0','PartType1','PartType3']:
            if pt in f:
                n = f[pt]['Coordinates'].shape[0]
                print(f"{pt} has {n} particles")
        if 'PartType0' in f:
            print("Sample vels:", f['PartType0']['Velocities'][:5])
            print("Sample coords:", f['PartType0']['Coordinates'][:5])



#inspect_snapshot('snapshot_038_cutout_split.hdf5')

def years_to_code_time(years, h):
    seconds = years * 365.25 * 24 * 3600
    unit_time = 3.085677580962325e16 / h   # your code time unit (s/h)
    return seconds / unit_time

fn = "snapshot_038_cutout_split.hdf5"

# Choose starting time
t_start_yrs = 5000
h = 0.6774
t_code = years_to_code_time(t_start_yrs, h)

with h5py.File(fn,'r+') as f:
    hdr = f['Header'].attrs
    print("Old Time:", hdr['Time'])
    hdr['Time'] = 1
    hdr['Redshift'] = 0.0
    hdr['HubbleParam'] = 1.0
    print("New Time:", hdr['Time'])
