import h5py, numpy as np

fn = "snapshot_038_cutout_split_physical.hdf5"   # set this to your file
with h5py.File(fn,'r') as f:
    hdr = f['Header'].attrs
    print("Header NumPart_ThisFile:", hdr.get('NumPart_ThisFile'))
    print("Header NumPart_Total:", hdr.get('NumPart_Total'))
    print("Header MassTable:", hdr.get('MassTable'))
    print("Header NumFilesPerSnapshot:", hdr.get('NumFilesPerSnapshot'))
    print("Header Time, Redshift, HubbleParam:", hdr.get('Time'), hdr.get('Redshift'), hdr.get('HubbleParam'))
    for pt in range(6):
        gname = f'PartType{pt}'
        if gname in f:
            grp = f[gname]
            print("\nGroup", gname, "has datasets:", list(grp.keys()))
            lengths = {}
            for d in grp.keys():
                arr = grp[d]
                try:
                    n = arr.shape[0]
                except Exception as e:
                    # might be scalar or attribute-like
                    n = None
                lengths[d] = n
            print(" lengths:", lengths)
        else:
            print("\nGroup", gname, "MISSING")


with h5py.File(fn,'r') as f:
    for pt in range(6):
        g = f.get(f'PartType{pt}')
        if g is None: continue
        if 'Masses' in g:
            print(pt, "Masses dtype:", g['Masses'].dtype, "shape:", g['Masses'].shape)
        else:
            print(pt, "no Masses; header mass table entry:", f['Header'].attrs['MassTable'][pt])


with h5py.File("snapshot_038_cutout_split_physical.hdf5","r") as f:
    ids = []
    for pt in f.keys():
        if pt.startswith("PartType") and "ParticleIDs" in f[pt]:
            ids.append(f[pt]["ParticleIDs"][:])
    ids = np.concatenate(ids)
    print(len(ids), len(np.unique(ids)))

    print("gas", f['PartType0']['Masses'].shape)
    print("dm ", f['PartType1']['Masses'].shape)
    print("bh ", f['PartType3']['Masses'].shape)

f =h5py.File("snapshot_038_cutout_split_physical.hdf5","r")
for pt in ["PartType0","PartType1","PartType3"]:
    if pt in f:
        print(pt, {d: f[pt][d].shape for d in f[pt].keys()})
print("Header NumPart_ThisFile:", f['Header'].attrs['NumPart_ThisFile'])
f.close()
