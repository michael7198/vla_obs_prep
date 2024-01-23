# Code to help with preparing VLA observations

This repo contains code I've used to assist in prepping observations with the VLA. In the examples below *targets* is the filename of the list of targets including columns for object name, ra and dec, and *date* is string of the desired observing date in the form yyyy-mm-dd. There are also a number of optional arguments that can be called for the python scripts in this repo.


_find_vla_calibrators.py_ can be used to point to a target list (ideally everything in a planned schedule block) and find potential complex gain calibrators for the targets of interest. By default it will also output an 'optimised' list of the minimum amount of calibrators that can be used which may be helpful when there are a large number of targets in the same area of the sky that can use the same calibrator. Run as:

    >python find_vla_calibrators.py targets


_VLA_calibrator_list.txt_ is a text file containing the list of complex gain calibrators for the VLA, in the same format as at https://science.nrao.edu/facilities/vla/observing/callist


_target_elevation.py_ is a script to output elevation vs LST plots for targets on a given date. Run as

    >python target_elevation.py targets date
    



<br/>

## Code Dependencies

This code was developed using Python 3.10 and the following packages are required to run (version used in development):
* argparse (1.1)
* astropy (5.0.2)
* distutils (3.10.12) 
* matplotlib (3.6.2)
* numpy (1.22.3)

