# Metal-ions-in-nanocrystals
## Overview
For the processing and analysis of the STEM-EDS data we wrote a Fiji script that does the following:
1. It identifies the borders of the nanocrystal and the metal atoms within it
2. It calculates the distance of the atoms to the nanocrystal edge
The script can be found at the [Fiji folder](../../tree/main/Fiji). 
The output of the Fiji script is then used by a Matlab script to do teh following:
1. Simulate the distribution of dopant atoms within the core of the nanocrystal and on its surface
2. The distribution of number of dopants per particle for each particle size is extracted.
