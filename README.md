This Python code was first developed in 2021 to process outputs from the UCL SGNL Orbit Prediction Software (OPS). Its purpose was to provide graphs comparing altitudes at various steps along real and simulated satellite orbits, for inclusion in a master's research project dissertation.

OPS numerically integrates along each step in a satellite's orbit and records the track in an ephemeris output file. The drag force, derived from atmospheric mass density and based on spatial and temporal coordinates, is calculated at each step. The magnitude of this force is the dominant factor determining the rate of decay of the satellite's orbit.

The purpose of the project was to compare the accuracy of predictions based on differing values of mass density from two atmosphere models, using differing configurations of various solar and geomagnetic input parameters.

This code is an extract only, published here as an example of work completed as part of the master's project. The file will not run without the OPS output files.
