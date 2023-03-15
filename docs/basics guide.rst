**HSFunctions : Basic functions for working with the h5 files output by the IMA Hyperspectral Microscope**   

Here you can import the h5 file into a workable array, extract relevant experimental information from the h5 file, and replicate simple images and plots. 

_________________________  

**Dependant Packages**  

.. code-block:: python  

	import numpy as np  
	import matplotlib.pyplot as plt  
	import h5py  
	import pyUSID as usid  
	from math import prod  

_________________________  

**User Input**  

.. code-block:: python  

	file = r'some_.h5_file'
	calibration_file = r'some_csv/txt_file'
	
_________________________  

**Data Import**  

.. code-block:: python  

	import_h5(file,mode = 'r')  #imports h5 file; default read-only  

	get_library(file) #prints the h5 tree  
_________________________  

**Extract Acquisition Information**  

.. code-block:: python  

	image_cube(file, transpose = True) #gets 3D array of monochromatic images  

	wavelengths(file) #gets 1D array of wavelengths, indexes correlated with z-dimension in image cube  

	data_cube(file,transpose=True) #outputs (wavelengths,image_cube) as tuple  

	get_exposure_time(file) #gets the exposure time (default in seconds)  

	get_datetime(file) #image acquisition time as YYYY/MM/DD HH:MM:SS  

	bounds(file) #makes array of {min wavelength, max wavelength, stepsize}  

	objective(file) #finds objective magnification (int)  

	processings(file) #prints list/makes array of processing steps (i.e. Subtraction, Rectification, Registration, etc)  

	pixel_size(file) #makes array of pixel {size, unit}  
_________________________  

**Extract Data Information**  

.. code-block:: python  

	get_cube_dimensions(data_cube) #outputs array {width,height,num wavelengths}  
_________________________  

**Data Manipulation**  

.. code-block:: python  

	wavelength_correction(data_cube,calibration_file,**kwargs) #corrects for detector sensitivity  

	crop_image_cube(image_cube, min_width=None, max_width=None, min_height=None, max_height=None) #spatial crop of the image_cube  

	crop_data_cube(data_cube, min_width=None, max_width=None, min_height=None, max_height=None, min_wavelength = None, max_wavelength = None) #spatial and wavelength crop of data_cube  
_________________________  

**Replicate PhySpec Plots**  

.. code-block:: python  

	plot_pixel(data_cube,position,fig_num=None,xlabel='None',ylabel=None,**kwargs) #plots spectrum at pixel (x,y); kwargs relevant to matplotlib.pyplot  

	average_spectrum(data_cube,fig_num = None,xlabel = None,ylabel = None,**kwargs) #plots spectrum of whole image; kwargs relevant to matplotlib.pyplot; input crop_data_cube for select region  
