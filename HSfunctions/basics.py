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
from mpl_toolkits.axes_grid1.anchored_artists import (AnchoredSizeBar)

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
    usid.hdf_utils.print_tree(file,rel_paths = True)
    
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
    return [file['Cube/Info/Cube'].attrs['LowerWavelength'][0],file['Cube/Info/Cube'].attrs['UpperWavelength'][0],file['Cube/Info/Cube'].attrs['WavelengthStep'][0]]

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
        return 'Objective Not Found'

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
    '''
    
    Parameters
    ----------
    data_cube : tuple
        output from data_cube()

    Returns
    -------
    y : integer
        y dimension of cube (formerly the width).
    x : integer
        x dimension of cube (formerly the height).
    z : integer
        wavelength dimension.

    '''
    y,x,z = data_cube[1].shape
    return y,x,z

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
    

    Parameters
    ----------
    image_cube : 3D array
        output from image_cube.
    min_width : integer, optional
        min x boundary. The default is None.
    max_width : integer, optional
        max x boundary. The default is None.
    min_height : integer, optional
        min y boundary. The default is None.
    max_height : integer, optional
        max y boundary. The default is None.

    Returns
    -------
    cropped_image_cube : 3D array
        crops image, often to capture homogenous illumination region (1024x1024).

    '''  

    cropped_image_cube = image_cube[min_height:max_height,min_width:max_width,:]
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
    Parameters
    ----------
    image_cube : 3D array
        output from image_cube.
    min_width : integer, optional
        min x boundary. The default is None.
    max_width : integer, optional
        max x boundary. The default is None.
    min_height : integer, optional
        min y boundary. The default is None.
    max_height : integer, optional
        max y boundary. The default is None.
    min_wavelength : integer, optional
        min z boundary. The default is None.
    max_wavelength : integer, optional
        max z boundary. The default is None.

    Returns
    -------
    cropped_data_cube : tuple 
        crops data_cube output spatially and/or by wavelength (data_cube is tuple of wavelength_array and image_cube) 
        Note: "data_cube" is a tuple of (wavelength_array,image_cube)
    
    '''
    try:
        wl_min = np.where(data_cube[0]==min_wavelength)[0][0]
    except: 
        wl_min = None
    try:
        wl_max = np.where(data_cube[0]==max_wavelength)[0][0]
    except:
        wl_max = None
    cropped_wl = data_cube[0][wl_min:wl_max]
    cropped_cube = data_cube[1][min_height:max_height,min_width:max_width,wl_min:wl_max]
    
    return cropped_wl,cropped_cube

    
#color map of data at a single wavelength slice in pixels
def image_slice(data_cube,single_wavelength,fig_num = None,xlabel = 'x (pixels)',ylabel = 'y (pixels)',**kwargs):
    '''
     Parameters
     ----------
     data_cube : tuple
         (wavelength_array,image_cube)
     single_wavelength : integer
         wavelength (nm).
     fig_num : integer, optional
         assign figure number. The default is None.
     xlabel : string, optional
         x-axis title. The default is 'x (pixels)'.
     ylabel : string, optional
         y-axis title. The default is 'y (pixels)'.
     **kwargs : var
         keyword arguments for matplotlib.pyplot.imshow.
    
     Returns
     -------
     2D Array of image at specified wavelength
     Plots color map of data at a single wavelength slice in pixels 
     Note: "data_cube" is a tuple of (wavelength_array,image_cube)
     Enter plot params in function
    
     '''
    wl = np.where(data_cube[0] == single_wavelength)[0][0]
    image = data_cube[1][:,:,wl]
    plt.figure(fig_num)
    plt.imshow(image,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return image

#plots spectrum at a single pixel

def plot_pixel(data_cube,position,fig_num=None,**kwargs):
    '''

    Parameters
    ----------
    data_cube : tuple
        (wavelength_array,image_cube).
    position : tuple
        (x,y)
    fig_num : integer, optional
        assigned to figure number. The default is None.
    # xlabel : string, optional
    #     x-axis title. The default is 'None'.
    # ylabel : string, optional
    #     y-axis title. The default is None.
    **kwargs : var
        keyword arguments for matplotlib.pyplot.plot.

    Returns
    -------
    2D array (wavelength_array,pixel spectrum)
    plots spectrum at a single pixel 

    '''
    plt.figure(fig_num)
    plt.plot(data_cube[0],data_cube[1][position[0],position[1],:],**kwargs)
    return data_cube[0],data_cube[1][position[0],position[1],:]

def average_spectrum(data_cube,fig_num = None,**kwargs):
    '''

    Parameters
    ----------
    data_cube : tuple
        (wavelength_array,image_cube).
    fig_num : integer, optional
        assigned to figure number. The default is None.
    # xlabel : string, optional
    #     x-axis title. The default is 'None'.
    # ylabel : string, optional
    #     y-axis title. The default is None.
    **kwargs : var
        keyword arguments for matplotlib.pyplot.plot.

    Returns
    -------
    2D array of (wavelength_array,average spectrum)
    plots average spectrum across all pixels/num_pixels in inputted data cube
    For smaller region, input cropped_data_cube

    '''
    num_pixels = prod(data_cube[1].shape[:2])
    avg_spectrum = np.sum(np.sum(data_cube[1],axis=0),axis=0)/num_pixels
    plt.figure(fig_num)
    plt.plot(data_cube[0],avg_spectrum,**kwargs)
    return data_cube[0],avg_spectrum

def x_line_profile(data_cube,wavelength,y_pos,trace_fig_num=None,map_fig_num=None,linecolor='dimgrey',**kwargs):
    '''

    Parameters
    ----------
    data_cube : tuple
        (wavelength_array,image_cube).
    wavelength : integer
        desired wavelength in nm.
    y_pos : integer
        pixel number for line.
    trace_fig_num : integer, optional
        figure number for 2D intensity trace. The default is None.
    map_fig_num : integer, optional
        figure number for sample colormap. The default is None.
    linecolor : matplotlib color, optional
        horizontal line color. The default is 'dimgrey'.
    **kwargs : var
        keyword arguments for matplotlib.pyplot.plot.

    Returns
    -------
    2D array(x position,intensity trace)
    fig1: intensity 2DLine plot
    fig2: sample colormap with horizontal line correlating with intensity trace

    '''
    wl_pos = np.where(data_cube[0]==wavelength)[0][0]
    trace = data_cube[1][y_pos,:,wl_pos]
    plt.figure(trace_fig_num)
    plt.plot(np.linspace(0,len(trace),len(trace)),trace,**kwargs)
    
    plt.figure(map_fig_num)
    plt.imshow(data_cube[1][:,:,wl_pos])
    plt.axhline(y_pos,color = linecolor)
    return np.linspace(0,len(trace),len(trace)),trace

def y_line_profile(data_cube,wavelength,x_pos,trace_fig_num=None,map_fig_num=None,linecolor='dimgrey',**kwargs):
    '''
    
    Parameters
    ----------
    data_cube : tuple
        (wavelength_array,image_cube).
    wavelength : integer
        desired wavelength in nm.
    x_pos : integer
        pixel number for line.
    trace_fig_num : integer, optional
        figure number for 2D intensity trace. The default is None.
    map_fig_num : integer, optional
        figure number for sample colormap. The default is None.
    linecolor : matplotlib color, optional
        DESCRIPTION. The default is 'dimgrey'.
    **kwargs : var
        keyword arguments for matplotlib.pyplot.plot.

    Returns
    -------
    2D array (y position, intensity trace)
    fig1: intensity 2DLine plot
    fig2: sample colormap with vertical line correlating with intensity trace

    '''
    wl_pos = np.where(data_cube[0]==wavelength)[0][0]
    trace = data_cube[1][:,x_pos,wl_pos]
    plt.figure(trace_fig_num)
    plt.plot(np.linspace(0,len(trace),len(trace)),trace,**kwargs)
    
    plt.figure(map_fig_num)
    plt.imshow(data_cube[1][:,:,wl_pos])
    plt.axvline(x_pos,color=linecolor)
    return np.linspace(0,len(trace),len(trace)),trace

def intensity_dist(image,facecolor='grey',fig_num = None,**kwargs):
    '''
    Parameters
    ----------
    image : 2D array
        Any 2D image, e.g. above outputs.
    facecolor : matplotlib color, optional
        Color of histogram. The default is 'grey'.
    **kwargs : var
        Keyword arguments forwarded to matplotlib.pyplot.hist.

    Returns
    -------
    2D Array (counts, bins, bars) from histogram

    '''
    hist_array = np.ravel(image)
    binsize = int(0.01*max(hist_array))
    bins = np.arange(hist_array.min(),hist_array.max()+binsize,binsize)
    plt.figure(fig_num)
    counts, bins, bars = plt.hist(hist_array,bins=bins,facecolor=facecolor,**kwargs)
    return counts,bins,bars

def axis2scalebar(pixel_size,um_length,fig_num=None,loc='lower right',pad = 0.1, borderpad = 0.5, sep=5, frameon=False,size_vertical=20,**kwargs):
    '''
    
    Parameters
    ----------
    pixel_size : float
        pixel size in um (can obtain from pixel_size(file)).
    fig_num : integer, optional
        figure number to assign scalebar to. The default is None.
    um_length : integer or float, optional
        physical length of scalebar in um.
    loc : string, optional
        upper left', 'upper center', 'upper right', 'center left', 'center', 'center right', 'lower left', 'lower center', 'lower right'. For backward compatibility, numeric values are accepted as well. The default is 'lower right'.
    pad : float, optional
        Padding around the label and size bar, in fraction of the font size. The default is 0.1.
    borderpad : float, optional
        Border padding, in fraction of the font size. The default is 0.5.
    sep : float, optional
        Separation between the label and the size bar, in points. The default is 5.
    frameon : boolean, optional
        If True, draw a box around the horizontal bar and label. The default is False.
    size_vertical : integer, optional
        Vertical length of the size bar, given in coordinates of transform. The default is 20.
    **kwargs : var
        Keyword arguments forwarded to AnchoredOffsetbox.

    Returns
    -------
    The input colormap without xy-axes and with a physical scalebar

    '''
    if fig_num != None:
        plt.figure(fig_num)
        plt.axis('off')
    else:
        plt.gca().axis('off')
    plt.fig(fig_num)
    bar_length = um_length/pixel_size
    bar = AnchoredSizeBar(plt.gca().transData,bar_length,str(um_length)+r' $\mu$m',loc=loc,pad=pad,borderpad=borderpad,sep=sep,frameon=frameon,size_vertical=size_vertical,**kwargs)
    plt.gca().add_artist(bar)
    return