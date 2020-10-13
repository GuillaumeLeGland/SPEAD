# SPEAD
Code and data for the phytoplankton community model SPEAD (Simulating Plankton Evolution with Adaptive Dynamics)

The temporary version currently available (SPEAD1.0_temp) is the version used to write the article submitted to Geoscientific Model Development by Guillaume Le Gland, Sergio M. Vallina, S. Lan Smith and Pedro Cermeno in early September 2020. It is built in MATLAB R2010b.
Its main module is jamstecrest_gausssecomodel1D.m. The model parameters can be modified in jamstecrest_gaussecomodel1D_keys and jamstecrest_gaussecomodel1D_parameters. 

A more user-friendly version (SPEAD_1.0) able to compute all the quantities and plot all the figures mentioned in the article was uploaded on the 13th of October 2020. For this definitive version, compatibility with GNU Octave has been checked. The execution has been tested on Windows with a 2.5 GHz Intel i5-3210M processor and on Linux Ubuntu with a 2.4 GHz Intel Xeon E5645 processor. SPEAD 1_0 requires the auxiliary folders present in the SPEAD1.0_temp directory.

Environmental variables (temperature, vertical diffusivity, mixed layer depth and Photosynthetically Available Radiation) and observational data to validate the model are located in the INPUTS Directory. The INPUTS directory must be at the same level as the SPEAD1.0_temp or SPEAD1.0 directory.

Please contact me (legland@icm.csic.es) if you cannot download, run or understand the model.
