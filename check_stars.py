import yt
import sys
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=67.74, Om0=0.315)
import numpy as np
import h5py

ds = yt.load("snapshot_020_cutout_split_physical.hdf5")
ad = ds.all_data()

z = ds['Redshift']
a = ds['Time']
h = ds['HubbleParam']

print(z)
print("Length unit:", ds.length_unit)
print("Mass unit:", ds.mass_unit)
print("Time unit:", ds.time_unit)
print("Velocity unit:", ds.velocity_unit)
print("Density unit:", ds.quan(1.0, "code_density").in_units("g/cm**3"))
print("Energy unit:", ds.quan(1.0, "code_specific_energy").in_units("km**2/s**2"))
#print("Energy unit:", ds.quan(1.0, "code_energy").in_units("km/s"))

# Check if 'PartType3', 'Masses' field exists
#if ('PartType4', 'Masses') not in ds.field_list:
 #   print("PartType4 particles not found in this snapshot. Exiting.")
  #  sys.exit(1)



for field in ds.field_list:
#    if field[0] == "PartType4":
        print(field)
        
#for key in ds.parameters:
 #   print(f"{key}: {ds.parameters[key]}")

print(ad[('PartType0','Density')])
print(np.max(ad[('PartType0','Density')]))
print(ad[('PartType0','Coordinates')])
print(ad[('PartType0','InternalEnergy')])
print(ad[('PartType0','InternalEnergy')].to('km**2/s**2'))
print(np.max(ad[('PartType0','InternalEnergy')]))
print(ad[('PartType0','Masses')])
print(ad[('PartType0','Masses')].to('Msun'))

#print(ad[('PartType4','ParticleIDs')])
#print(ad[('PartType4','Stellar Age')])
#print(ad[('PartType4','Star Type')])
#print(ad[('PartType4','StellarFormationTime')])
#print(1/ad[('PartType4','StellarFormationTime')]-1)

#print(ad[('PartType4','Masses')].to('Msun/h'))
# 586.05725331


#print(ad[('PartType0','Masses')].to('Msun'))
#print(ad[('PartType3','Masses')].to('Msun'))

#print(ad[('PartType0','Velocities')]/np.sqrt((1+z)))
print(ad[('PartType0','Velocities')].to('km/s'))
