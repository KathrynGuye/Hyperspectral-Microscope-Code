Basic functions for working with the h5 files output by the IMA Hyperspectral Microscope  
Here you can import the h5 file into a workable array, extract relevant experimental information from the h5 file, replicate simple images and plots, and extract spectral features. 

***

**Dependencies**

    import numpy as np  
    import matplotlib.pyplot as plt  
    import h5py  
    import pyUSID as usid  

***

    file = some_.h5_file  
**Data Import**  
`import_h5(file,mode = 'r')` | imports h5 file; default read-only  
`get_library(file)` | prints the h5 tree

