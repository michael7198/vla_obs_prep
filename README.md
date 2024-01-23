# Code to help with preparing VLA observations

This repo contains code I've used to assist in prepping observations with the VLA. 


find_vla_calibrators.py can be used to point to a target list (ideally everything in a planned schedule block) and find potential complex gain calibrators for the targets of interest. By default it will also output an 'optimised' list of the minimum amount of calibrators that can be used which may be helpful when there are a large number of targets in the same area of the sky that can use the same calibrator. Run as:

    >python find_vla_calibrators.py targets

where *targets* is the filename of the list of targets including columns for object name, ra and dec. There are also a number of optional arguments that can be called.


<br/>

## Code Dependencies

This code was developed using Python 3.10 and the following packages are required to run (version used in development):
* argparse (1.1)
* astropy (5.0.2)
* distutils (3.10.12) 
* matplotlib (3.5.1)
* numpy (1.22.3)

