#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astroquery.vizier import Vizier
from astropy.io import fits
from matplotlib.colors import LogNorm
from astropy import wcs

# Your base directory here
base_dir = '/Users/esola/Documents/Postdoc/Lectures/DIS_MPhil_2025/notebook/'


# # <a  anchor='#chapter1'> -*-*-*-*- I. Catalogues -*-*-*-*- </a>

# ## <a  anchor='section1_1'> I.1) Direct download of online tables </a>
# 

# In[2]:


# Either you found the table somewhere already...

# Or you can use services such as Vizier to access published catalogues

catalogue_name = 'J/MNRAS/506/5494' # Poulain+2021 MATLAS dwarfs structure and morphology
# https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/MNRAS/506/5494&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
# Download with 'Unlimited' rows the catalogue as VOTable file (and rename it as poulain21.vot)


# In[3]:


## Read the downloaded table
catalogue_path = base_dir + 'poulain21.vot'
catalogue      = Table.read(catalogue_path, format = 'votable', table_id = 0)
catalogue


# ## <a  anchor='section1_2'> I.2) Access with the astroquery module </a>
# 

# In[4]:


catalogue_name = 'J/MNRAS/506/5494' # Poulain+2021 MATLAS dwarfs structure and morphology
catalogue      = Vizier.get_catalogs(catalogue_name) 
print(catalogue)


# In[5]:


## By default, Vizier only queries the 50 first rows. To retrieve all the rows in a query, set ROW_LIMIT to -1
Vizier.ROW_LIMIT = -1
catalogue = Vizier.get_catalogs(catalogue_name) 
print(catalogue)


# In[6]:


## Access the 1st table
table1 = catalogue[0]
table1


# In[7]:


## List the columns names
print(table1.colnames) # Similar to print(table1.keys())


# In[8]:


## Quick statistics of the dwarf morphologies
# dE: dwarf elliptical
# dI: dwarf irregular
# dEN: dwarf elliptical nucleated
# dIN: dwarf irregular nucleated
dE  = table1[np.where(table1['Morph']=='dE')]
dI  = table1[np.where(table1['Morph']=='dI')]
dEN = table1[np.where(table1['Morph']=='dEN')]
dIN = table1[np.where(table1['Morph']=='dIN')]

print('Total number of dwarfs', len(table1))
print('Total number of dE', len(dE), 'of dI', len(dI), 'of dEN', len(dEN), 'of dIN', len(dIN))
if (len(dE) + len(dI) + len(dEN) + len(dIN)) != len(table1):
    print(len(table1) -(len(dE) + len(dI) + len(dEN) + len(dIN)), 'dwarfs do not have a dE,dI,dEN or dIN morphology ')


# In[9]:


## Identify these dwarfs that do not have dE,dI,dEN or dIN morphology
nomorph = table1[np.where((table1['Morph']!='dE') & (table1['Morph']!='dI') & (table1['Morph']!='dEN') & (table1['Morph']!='dIN'))]
print(nomorph)

# Other way of quickly checking the different values of a column
print('Possible morphologies are',set(table1['Morph']))

# > Look at the paper to understand what this means
#dE∗ corresponds to the dwarf whose nucleus matches with a star in the Gaia DR2 catalogue.
#dE∗∗ corresponds to the dwarf for which one of the central bright GCs is possibly an NSC (Forbes et al. 2019)


# ## <a  anchor='section1_3'> I.3) Examples of plots </a>
# 

# In[10]:


# Plot the positions of the dwarfs in RA,DEC
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.scatter(table1['RAJ2000'], table1['DEJ2000'], color='k', marker='.')
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
# Convention: RA increases to the left
plt.gca().invert_xaxis()

plt.subplot(1,2,2)
all_dE = table1[np.where((table1['Morph']=='dE') | (table1['Morph']=='dEN'))]
all_dI = table1[np.where((table1['Morph']=='dI') | (table1['Morph']=='dIN'))]

plt.hist(all_dE['gmag'], bins=17, color='b', edgecolor='b', label='dE/dEN', alpha=0.5)
plt.hist(all_dI['gmag'], bins=17, color='r', edgecolor='r',label='dI/dIN', alpha=0.5)
plt.xlabel('M$_{g}$ [mag]')
plt.ylabel('Number')
plt.legend()


# ## <a  anchor='section1_4'> I.4) Manipulate tables with Topcat </a>
# 

# #### Open the catalogue, double-click on it on the left panel to open the table

# <img src="notebook/images/topcat_table.png" width="1600">
# 

# ### Select the subset of galaxies that are dE/dEN and dI/dIN
# Go to Subset > click on the green + and add the expression (beware of the syntax for str)

# <img src="notebook/images/topcat_subset.png" width="800">
# 

# Create the histogram of the gmag distribution: 
# - Click on the histogram button 
# - In the bottom panel select 'Subset' and select dE/dEN and dI/dIN
# - Adjust the axes title and the label

# <img src="notebook/images/topcat_histogram.png" width="1100">
# 

# Plot the RA,Dec position of the sources:
# - Click on the 'Sky plotting window' icon (the celestial grid)
# - Select RA,Dec columns
# - Do it for all sources, then separating dE/dEN from dI/dIN

# <img src="notebook/images/topcat_scatter.png" width="800">
# 

# # <a  anchor='chapter2'> -*-*-*-*- II. Manipulate FITS files -*-*-*-*- </a>
# 

# ## <a  anchor='section2_1'> II.1) FITS and astropy </a>
# 

# In[11]:


fits_path = base_dir + 'matlas_ngc3230.fits' # 'UGC01245_r.fits'
# Opens the fits file: Header-Data Unit (HDU) contains the header and the image
with fits.open(fits_path) as hdu:
    header = hdu[0].header
    data   = hdu[0].data
    
## Other method to open fits (be careful, less memory efficient)
#header = fits.getheader(fits_path)
#data   = fits.getdata(fits_path)

#header


# In[12]:


# Size of the image
print(data.shape)
# Show the image
plt.imshow(data, cmap='gray')
# What is wrong with it ?


# In[13]:


## Adjust the scaling
plt.imshow(data, origin='lower', norm=LogNorm(), cmap='gray')
plt.xlabel('X [pix]')
plt.ylabel('Y [pix]')


# ## <a  anchor='section2_2'> II.2) Overplot objects' positions on an image </a>
# 

# You have a catalogue of coordinates (RA,DEC) of sources and a fits image. You will plot the position of the sources on that image

# In[14]:


# Opens the FITS image
fits_path = base_dir + 'matlas_ngc3230.fits'
with fits.open(fits_path) as hdu:    
    header = hdu[0].header
    data   = hdu[0].data
    
plt.figure(figsize=(5,5))
plt.imshow(data, origin='lower', cmap='gray', norm=LogNorm())
plt.xlabel('X [pix]')
plt.ylabel('Y [pix]')


# Opens the catalogue
catalogue = Table.read(base_dir + 'matlas_dwarfs_on_ngc3230.txt', format='ascii')
print('Number of sources in the catalogue',len(catalogue))
catalogue
# Initialises a WCS (World Coordinate System) instance with the astronomical calibration of the header


# In[17]:


### To go from the RA,DEC coordinates of the objects to the X,Y pixels in the image above, you need to use the WCS (World Coordinate System)

# Initialise a WCS instance
my_wcs = wcs.WCS(header)
# Examine what the WCS looks like 
print(my_wcs)

# Now that we have the WCS transform element, we can go from (RA,DEC) to (X,Y)
ra  = list(catalogue['RAJ2000'])
dec = list(catalogue['DEJ2000'])
x,y = my_wcs.world_to_pixel_values(ra,dec)

# Iterates through the list of (X,Y) coordinates and overplot them on the image
plt.figure(figsize=(7,7))
plt.imshow(data, origin='lower', cmap='gray', norm=LogNorm())
for i in range(len(x)):
    plt.scatter(x,y, marker='.', color='r')
plt.xlabel('X [pix]')
plt.ylabel('Y [pix]')


# #### Quick check that the positions are correct with Aladin 

# <img src="notebook/images/matlas_dwarfs_on_ngc3230.png" width="600">
# 
# 

# ## <a  anchor='section2_3'> II.3) FITS and DS9 </a>
# 

# - Click on File > Open and select the file
# - Click on File > Header to display the header
# - Change the scaling of the image with scale > histogram (try some others)

# <img src="notebook/images/fits_header.png" width="1200">
# 

# # <a  anchor='chapter3'> -*-*-*-*- III. Aladin -*-*-*-*- </a>
# 

# ### Goal: 
# 1) Load a set of publicly available images in Aladin. You will use MATLAS color images
# 2) Load a catalogue of galaxies. Here you wil
# 3) Find which galaxies are located in the images you have. Export the resulting table and open it with Python

# ## <a  anchor='section3_1'> III.1) Load survey images </a>
# 

# - In the bottom left panel, write MATLAS in the selection box.
# - Then select just above from Collections >  Images > Optical > MATLAS > MATLAS_color
# 

# <img src="notebook/images/matlas_load.png" width="800">

# ## <a  anchor='section3_2'> III.2) Load a catalogue </a>
# 

# - In the bottom left panel, write ATLAS3D
# - Select just above from Collections > Catalog > Vizier > Journal Table > MNRAS > the ATLAS3D project I (Cappellari+2011)
# - Check that the catalogue contains 871 sources (hover your mouse on the catalogue name and look for the number of rows)

# <img src="notebook/images/atlas3d_load.png" width="800">

# You should have two layers loaded on the right-side of the window: the image layer and the catalogue (red triangles)

# <img src="notebook/images/matlas_atlas3d.png" width="800">
# 
# And if you zoom out and display the coordinate grid
# 
# <img src="notebook/images/matlas_atlas3d_zoom_out.png" width="800">
# 
# 

# ## <a  anchor='section3_3'> III.3) Create the MOC (Multi-Order Coverage map) of the images </a>
# 

# - Make sure you selected the MATLAS image layer on the bottom right panel
# - Go to Coverage > Load MOC for current survey
# - The newly created MOC appears as another image layer
# 
# <img src="notebook/images/moc.png" width="800">
# 

# ## <a  anchor='section3_4'> III.4) Intersection between MOC and catalogue </a>
# 

# - Go to Coverage > Filter a table by MOC
# - Correctly select the MOC and catalogue layers
# - A new catalogue resulting from the intersection between the MOC and the original catalogue is created (it appears on the bottom right panel)
# - Export this new catalogue by File > Export layers (FITS, VOTable) and select the layer to export + its output format
# 

# ### Now, load this catalogue with python and inspect its content

# In[18]:


intersection_cat = Table.read(base_dir + 'atlas3d_on_moc_matlas.xml', format='votable')
intersection_cat


# # <a  anchor='chapter4'> -*-*-*-*- IV. Exploration of Euclid data -*-*-*-*- </a>
# 

# Goals with Aladin:
# - Load the Euclid ERO images
# - Search for the Perseus cluster
# - Create the MOC of the ERO images 
# - Overlay the Simbad catalogue in the current MOC
# - Export the catalogue

# <img src="notebook/images/ero_perseus.png" width="1700">

# Visually compare Euclid ERO in VIS band to other surveys in the same view (make sure to click on the small button 'Unif.' at the very bottom of the view)
# 
# <img src="notebook/images/euclid_ero_vs_other.png" width="1700">
# 

# #### Note: Aladin Lite is a lightweight version of Aladin with few options, but you can still do quite a lot quickly !
# 
# <img src="notebook/images/aladin_lite.png" width="1600">

# # <a  anchor='chapter5'> -*-*-*-*- V. Exploration of DESI Legacy Imaging Surveys data -*-*-*-*- </a>
# 

# ## Accessible at: https://www.legacysurvey.org/viewer 

# Goals:
# - Query some objects you want to see
# - Navigate in the images
# - Change the constrat/brightness of the images
# - Try out the different overlay layers

# <img src="notebook/images/desi_ls.png" width="1600">
