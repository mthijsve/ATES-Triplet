# ATES-Triplet is a repository for an ATES-Triplet model that is part of the PhD research of M.S.(Thijs) van Esch.
# Based heavily on the PySeawATES code developed by Bloemendal, Beernink, Olsthoorn and Jaxa-Rozen, that can be found [here](https://github.com/martinbloemendal/PySeawATES/)
The repository consists of the following files:
setup files:
-environment.yml
	is used to install all needed python packages
-swt_v4x64.exe
	is the executable that runs seawat, called in the code
input files:
-sub-surface.csv
	enter sub-surface parameters in this csv, layer by layer
-wells.csv
	enter well data in this csv
-demand.csv
	optional, enter energy demand for the entire run
-monitoring.csv
	optional, enter all points where heads and temperatures should be saved
function files
-PySeawaTriplet.py
	runs the main portion of the code, including the loop where the temperature transport and flow is calculated
-functions_triplet.py
	contains all agent and grid functions from PySeawATES adapted for the Triplet
-flow_function_Triplet.py
	contains the function that calculates the energy profile from the yearly demand.
-Executing_file.py
	contains the grid and timestep parameters as well as the last setup variables needed. This file is used to run the model
