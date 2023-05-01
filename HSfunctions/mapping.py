# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 20:16:37 2023

@author: Kathryn Guye
Version: 1.0.0
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sgf


def normalize(y):
    normal_data = (y-min(y))/(max(y)-min(y))
    return normal_data
##############################################################################################################
def peak_map(data_cube,fig_num = None,xlabel = 'x (pixels)',ylabel = 'y (pixels)',cmap = 'inferno',colorbar = True,label = 'Wavelength (nm)',**kwargs):
    wavelength,images = data_cube
    def peak_max(y):
        _y = sgf(y, 5, 1)
        return np.argmax(_y, axis=-1)
    pl_peak = np.zeros_like(images[:,:,0])
    pl_peak[:,:] = [wavelength[peak_max(x)] for x in images[:,:,:]]
    plt.figure(fig_num)
    plt.imshow(pl_peak,cmap = cmap,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colorbar == True:
        plt.colorbar(label = label)
    average = np.average(pl_peak[np.isfinite(pl_peak)])
    std = np.std(pl_peak[np.isfinite(pl_peak)])
    return pl_peak[:,:], average, std


def peak_intensity_map(data_cube,fig_num = None,xlabel = 'x (pixels)',ylabel = 'y (pixels)',cmap = 'inferno',colorbar = True,label = 'Peak Intensity (counts)',**kwargs):
    wavelength,images = data_cube
    x,y,z = images.shape
    peak_int = np.ones((x,y))
    def max_int(y):
        _y = sgf(y,5,1)
        peak_int = np.max(_y,axis = -1)
        return peak_int
    peak_int = max_int(images)
    plt.figure(fig_num)
    plt.imshow(peak_int,cmap = cmap,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colorbar == True:
        plt.colorbar(label = label)
    average = np.average(peak_int[np.isfinite(peak_int)])
    std = np.std(peak_int[np.isfinite(peak_int)])
    return peak_int[:,:], average, std

def integrated_intensity_map(data_cube,fig_num = None,xlabel = 'x (pixels)',ylabel = 'y (pixels)',cmap = 'inferno',colorbar = True,label = 'Integrated Intensity (counts)',**kwargs):
    wavelength_array = data_cube[0]
    intensities = np.sum(data_cube[1],axis=2)
    plt.figure(fig_num)
    plt.imshow(intensities,cmap = cmap,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colorbar == True:
        plt.colorbar(label = label)
    average = np.average(intensities[np.isfinite(intensities)])
    std = np.std(intensities[np.isfinite(intensities)])
    return intensities[:,:], average, std

def linewidth_map(data_cube,fig_num = None,xlabel = "x (pixels)",ylabel = "y (pixels)",cmap = 'inferno',colorbar = True,label = 'Linewidth (nm)',**kwargs):
    def smooth(y,boxcar):
        smoothed = sgf(y,boxcar,1)
        return smoothed
    def cubenormalize(y):
        min_cube = np.repeat(np.min(y,axis = 2)[:,:,np.newaxis],y.shape[2],axis=2)
        max_cube = np.repeat(np.max(y,axis = 2)[:,:,np.newaxis],y.shape[2],axis=2)
        norm = (y-min_cube)/(max_cube-min_cube)
        return norm
    wl,images = data_cube
    stepsize = wl[1]-wl[0]
    x,y,z = images.shape
    fwhm_image = np.ones((x,y))
    spx = smooth(images,11)
    npx = cubenormalize(spx)
    for i in range(x):
        for j in range(y):
            try:
                a = npx[i,j,:]
                peak_pos = np.argmax(a)
                l = a[:peak_pos]
                r = a[peak_pos:]
                crossings_left = np.argwhere(np.diff(np.sign(0.5-l))).flatten()
                crossings_right = np.argwhere(np.diff(np.sign(0.5-r))).flatten()
                idx1 = crossings_left[0]
                idx2 = crossings_left[0]+1
                idx3 = crossings_right[0]+peak_pos
                idx4 = crossings_right[0]+peak_pos+1
                left = np.interp(0.5,[a[idx1],a[idx2]],[idx1,idx2])
                right = np.interp(0.5,[a[idx3],a[idx4]],[idx3,idx4])
                fwhm = stepsize*right - stepsize*left
                fwhm = np.round(fwhm,decimals = 1)
                fwhm_image[i,j] = fwhm
            except:
                fwhm_image[i,j] = np.nan
    plt.figure(fig_num)
    plt.imshow(fwhm_image,cmap=cmap,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colorbar == True:
        cbar = plt.colorbar(label = label)
    average = np.average(fwhm_image[np.isfinite(fwhm_image)])
    std = np.std(fwhm_image[np.isfinite(fwhm_image)])
    return fwhm_image[:,:], average, std

def fast_linewidth(data_cube,fig_num = None,xlabel = "x (pixels)",ylabel = "y (pixels)",cmap = 'inferno',colorbar = True,label = 'Linewidth (nm)',**kwargs):
    def smooth(y,boxcar):
        smoothed = sgf(y,boxcar,1)
        return smoothed
    def cubenormalize(y):
        min_cube = np.repeat(np.min(y,axis = 2)[:,:,np.newaxis],y.shape[2],axis=2)
        max_cube = np.repeat(np.max(y,axis = 2)[:,:,np.newaxis],y.shape[2],axis=2)
        norm = (y-min_cube)/(max_cube-min_cube)
        return norm
    wavelength,images = data_cube #splits up data cube
    stepsize = wavelength[1]-wavelength[0] #finds wavelength stepsize
    x,y,z = images.shape
    spx = smooth(images,11) #super-smoothes plot to prevent noise interference
    npx = cubenormalize(spx) #normalizes to find true 50%
    mins = np.abs(0.5-npx) #makes half max local minima
    ordered = np.argsort(mins,axis=-1)[:,:,:2] #finds the two lowest indecies - common place of error
    ordered = np.sort(ordered,axis = -1) #makes the left-most index come first
    right = stepsize*ordered[:,:,1] #2D array of right intersection
    left = stepsize*ordered[:,:,0] #2D array of left intersection
    fwhm_image = right-left
    
    plt.figure(fig_num)
    plt.imshow(fwhm_image,cmap=cmap,**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if colorbar == True:
        cbar = plt.colorbar(label = label)
    print("User Warning: fast_linewidth for observation only; data not reliable for publication")
    average = np.average(fwhm_image[np.isfinite(fwhm_image)])
    std = np.std(fwhm_image[np.isfinite(fwhm_image)])
    return fwhm_image[:,:], average, std








