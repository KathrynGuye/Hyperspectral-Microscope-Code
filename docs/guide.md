**Basic functions for working with the h5 files output by the IMA Hyperspectral Microscope**  
Here you can import the h5 file into a workable array, extract relevant experimental information from the h5 file, replicate simple images and plots, and extract spectral features. 

***

**Dependencies**

	import numpy as np  
	import matplotlib.pyplot as plt  
	import h5py  
	import pyUSID as usid  

***
**User Input** 
	file = r'some_.h5_file'  
	
**Data Import**  
`import_h5(file,mode = 'r')` | imports h5 file; default read-only  
`get_library(file)` | prints the h5 tree  

**Extract Acquisition Information**  
`image_cube(file, transpose = True)` | gets 3D array of monochromatic images  
`wavelengths(file)` | gets 1D array of wavelengths, indexes correlated with z-dimension in image cube  
`data_cube(file,transpose=True)` | outputs (wavelengths,image_cube) as tuple  
`get_exposure_time(file)` | gets the exposure time (default in seconds)  
`get_datetime(file)` | image acquisition time as YYYY/MM/DD HH:MM:SS  
`bounds(file)` | makes array of {min wavelength, max wavelength, stepsize}  
`objective(file)` | finds objective magnification (int)  
`processings(file)` | prints list/makes array of processing steps (i.e. Subtraction, Rectification, Registration, etc)  
