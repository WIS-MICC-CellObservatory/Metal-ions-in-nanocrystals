# Metal Ions Dynamic Nuclear Polarization in Mn(II) Doped CdS Nanocrystals: Atomic Scale Investigation of the Dopant and its Host
Ran E. Abutbul†, Daniel Jardon-Alvarez†, Stephen Fox‡, Lothar Houben‡, Ofra Golani#, Ehud Sivan#, Raanan Carmieli‡, Ilia Kaminker§, Michal Leskes†*
## Scripts Overview
For the processing and analysis of the STEM-EDS data we wrote a Fiji script that does the following:
1. It identifies the borders of the nanocrystal and the metal atoms within it
2. It calculates the distance of the atoms to the nanocrystal edge

The script can be found at the [Fiji folder](../../tree/main/Fiji). 

The output of the Fiji script is then used by a MATLAB script to do the following:
1. Simulate the distribution of dopant atoms within the core of the nanocrystal and on its surface
2. The distribution of number of dopants per particle for each particle size is extracted. 

The scripts can be found at the [MATLAB folder](../../tree/main/Matlab). 
