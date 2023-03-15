# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 15:29:30 2023

@author: Kathryn Guye

Version: 1.0.0
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import pyUSID as usid
from math import prod
from scipy.interpolate import InterpolatedUnivariateSpline

#import h5 file, default mode is r
def import_h5(file,mode = 'r'):
    '''
    import h5 file, default mode is r: read only
    '''
    file = h5py.File(file,mode)
    return file

#prints out a library of data cube locations
def get_library(file):
    '''
    prints out a library of data cube locations
    '''
    list(usid.hdf_utils.print_tree(file,rel_paths = True))
    
#outputs the 3D image cube, transposed is default ****note many further functions have dependency on transposed images****
def image_cube(file, transpose = True):
    '''
    outputs the 3D image cube, transposed is default
    '''
    image_cube = np.array(file['Cube/Images'])
    if transpose == True:
        image_cube = image_cube.T       
    return image_cube

#array of wavelengths correlated with z-dimension in image cube
def wavelengths(file):
    '''
    array of wavelengths correlated with z-dimension in image cube
    '''
    wavelength_array = np.array(file['Cube/Wavelength'])
    return wavelength_array

#outputs both wavelength and datacube; i.e. combines image_cube and wavelengths
def data_cube(file, transpose = True):
    '''
    outputs both wavelength and datacube; i.e. combines image_cube and wavelengths
    '''
    image_cube = np.array(file['Cube/Images'])
    if transpose == True:
        image_cube = image_cube.T
    wavelength_array = np.array(file['Cube/Wavelength'])
    return wavelength_array,image_cube
                
    
#finds the exposure time (default in seconds)
def get_exposure_time(file):
    '''
    finds the exposure time (default in seconds)
    '''
    return file['Cube/Info/Cube'].attrs['FixedTimeExposure'][0]

#gets the image acquision time as YYYY/MM/DD HH:MM:SS
def get_datetime(file):
    '''
    image acquision time as YYY/MM/DD HH:MM:SS
    '''
    return file['Cube/Info/Cube'].attrs['CreationDate'][0]

#makes array of [min wavelength,max wavelength,stepsize]
def bounds(file):
    '''
    makes array of [min wavelength,max wavelength,stepsize]
    '''
    return np.array(file['Cube/Info/Cube'].attrs['LowerWavelength'][0],file['Cube/Info/Cube'].attrs['UpperWavelength'][0],file['Cube/Info/Cube'].attrs['WavelengthStep'][0])

#gets objective magnification as int
def objective(file):
    '''
    objective magnification as int
    '''
    return int(file['Cube/Info/Optics'].attrs['Objective'][0][:2])

#gets pixel size as length,unit
def pixel_size(file):
    '''
    gets pixel size as length,unit
    '''
    try:
        return file['Cube/Info/Camera/YAxis/1'].attrs['Coefs'][1],file['Cube/Info/Camera/YAxis/1'].attrs['Unit'][0]
    except:
        return [file['Cube/Info/Processings/000'].attrs['Name']]

#prints a list of the processing steps performed on cube
def processings(file):
    '''
    prints a list of the processing steps performed on cube
    '''
    steps = np.array(file['Cube/Info/Processings'])
    processings = []
    # processings = processings
    for i in range(len(steps)):
        # print(file['Cube/Info/Processings/'+steps[i]].attrs['Name'][0])
        processings.append(file['Cube/Info/Processings/'+steps[i]].attrs['Name'][0].astype(str))
    print("Processing Steps: " + str(processings))
    return processings

'''please note, further definitions depend on a transposed datacube'''

def get_cube_dimensions(data_cube):
    w,h,z = data_cube[1].shape
    return w,h,z

def wavelength_correction(data_cube,calibration_file,**kwargs):
    '''
      Parameters
    ----------
    data_cube : tuple (wavelength_array,image_cube)
    calibration_file : csv or txt
        first column is wavelength; second column is calibrated intensity readings (see SOP)
    **kwargs : keywords
        any keyword entries for np.genfromtxt; especially skip_header, delimiter

    Returns
    -------
    wl : 1D array
        wavelength_array
    data_corrected : 3D array
        corrected data_cube

    '''
    calibration = np.genfromtxt(calibration_file,**kwargs)
    cal_x = calibration[:,0]
    cal_y = calibration[:,1]
    wl = data_cube[0]
    stepsize = wl[1]-wl[0]
    xs = np.linspace(cal_x[0],cal_x[-1],int((cal_x[-1]-cal_x[0])/stepsize)+1)
    spl = InterpolatedUnivariateSpline(cal_x, cal_y)
    correction = [xs,spl(xs)]
    start_index = np.where(xs == wl[0])[0][0]
    end_index = np.where(xs == wl[-1])[0][0]
    coeff_corr = correction[1][start_index:end_index+1]
    data_corrected = data_cube[1]
    for i in range(len(coeff_corr)):
        data_corrected[:,:,i] = data_corrected[:,:,i]/coeff_corr[i]
    return wl,data_corrected

def crop_image_cube(image_cube, min_width=None, 
                    max_width=None, 
                    min_height=None, 
                    max_height=None):
    '''
    crops image, often to capture homogenous illumination region (1024x1024)
    '''
    cropped_image_cube = image_cube[min_width:max_width,min_height:max_height,:]
    return cropped_image_cube

#crops data_cube output spatially and/or by wavelength (data_cube is tuple of wavelength_array and image_cube)
def crop_data_cube(data_cube,
                   min_width=None, 
                   max_width=None, 
                   min_height=None, 
                   max_height=None,
                   min_wavelength=None,
                   max_wavelength=None):
    '''
    crops data_cube output spatially and/or by wavelength (data_cube is tuple of wavelength_array and image_cube) 
    Note: "data_cube" is a tuple of (wavelength_array,image_cube)
    '''
    wl_min = np.where(data_cube[0]==min_wavelength)[0][0]
    wl_max = np.where(data_cube[0]==max_wavelength)[0][0]
    cropped_wl = data_cube[0][wl_min:wl_max]
    cropped_cube = data_cube[1][min_width:max_width,min_height:max_height,wl_min:wl_max]
    return cropped_wl,cropped_cube

#color map of data at a single wavelength slice in pixels
def image_slice_pixels(data_cube,single_wavelength,fig_num = None,xlabel = 'x (pixels)',ylabel = 'y (pixels)',**kwargs):
    '''
    Color map of data at a single wavelength slice in pixels 
    Note: "data_cube" is a tuple of (wavelength_array,image_cube)
    Enter plot params in function
    '''
    wl = np.where(data_cube[0] == single_wavelength)[0][0]
    plt.figure(fig_num)
    plt.imshow(data_cube[1][:,:,wl],**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

#plots spectrum at a single pixel

def plot_pixel(data_cube,position,fig_num=None,xlabel='None',ylabel=None,**kwargs):
    '''
    plots spectrum at a single pixel 
    Note: "data_cube" is a tuple of (wavelength_array,image_cube)
    position: (x,y)
    Enter plot params in function
    '''
    plt.figure(fig_num)
    plt.plot(data_cube[0],data_cube[1][position[0],position[1],:],**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return

def average_spectrum(data_cube,fig_num = None,xlabel = None,ylabel = None,**kwargs):
    '''
    plots average spectrum across all pixels/num_pixels in inputted data cube
    For smaller region, input cropped_data_cube
    Note: "data_cube" is a tuble of (wavelength_array,image_cube)
    Enter plot params in function
    '''
    num_pixels = prod(data_cube[1].shape[:2])
    avg_spectrum = np.sum(np.sum(data_cube[1],axis=0),axis=0)/num_pixels
    
    plt.figure(fig_num)
    plt.plot(data_cube[0],avg_spectrum,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return
































