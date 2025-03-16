#!/usr/bin/env python
# coding: utf-8

# # DIS GA Minor Module - supervision practice exercise

# -------
# 
# **The goal of this exercise is not to write the best code ever, but to explore the data.**  
# 
# *We will go through how to approach the exercise in the second-to-last lecture on March 14 (the usual room). Please try the exercise beforehand, and see how far you get.*
# 
# -------
# 
# -------
# 
# -------
# 
# **The tasks:**
# 
# **1. Make an account on the Gaia archive (https://gea.esac.esa.int/archive/) and familiarise yourself with ADQL queries**  
#     - To get to the ADQL query page, go to the "Search" tab, then click "Advanced (ADQL)"  
#     - There are many query examples available on the Gaia archive page (top right of the page: "Query examples")  
#     - you can do queries on the Gaia website, but it may be useful to also try out the Python *astroquery* package to be able to do it in Python (optional)  
# 
# 
# **2. Query the Gaia DR3 archive using ADQL and download the Gaia data around the globular cluster NGC 3201**  
#     - hint: you can use Simbad (https://simbad.cds.unistra.fr/simbad/sim-fid) to find out more information about this object (and any other astronomical object you may be interested in), e.g. the coordinates, proper motions, parallax, radial velocity (if available), spatial extent, references to papers  
#     - hint: what is the spatial extent of a globular cluster? So how large a radius (e.g. in degrees) would be useful to query around the coordinates of the cluster? And do you need all the available Gaia columns or only a subset?   
# 
# 
# **3. Create a colour-magnitude diagram (CMD) for this cluster**  
#     - Remember that *lower* magnitudes means *brighter* stars  
#     - Try removing non-members using the proper motions (and maybe the parallax too). All stars in a globular cluster are expected to have the same proper motions and parallax (but of course the measurements come with uncertainties, which are larger for fainter stars).   
#     - Identify the different evolutionary stages of stars on the CMD  
#     - Remember that it is usually necessary to deredden the magnitudes. See below for some hints. What are the differences between the reddened and dereddened CMDs?  
# 
# 
# **4. *(optional)* Derive and plot the Galactic orbit of this cluster**  
#     - You can use the *galpy* package for the orbit integration, which has very good documentation (https://docs.galpy.org/en/v1.9.1/orbit.html). You can use the standard *galpy* Milky Way potential (MWPotential2014).  
#     - How does this cluster move with respect to the Galactic disc?
#     
#     
# If you'd like more practice, you can repeat the exercise for another globular cluster. NGC 3201 was a "nice one" because it has relatively low central density and is nearby and therefore the stars have well-measured parallaxes. Other clusters may suffer from bad Gaia data quality in the central regions because they are so dense, and are typically further away and therefore have larger uncertainties on measurements. 
#  
# -------
#  
# -------
#  
# -------

# **Some general hints:**
# 
# - Remember the typical observers astronomical coordinate system, with coordinates right ascension (RA) and declination (Dec). This is the ICRS coordinate system. A common way of dealing with coordinates in Python is with the *astropy* coordinates package (https://docs.astropy.org/en/stable/coordinates/index.html - *astropy* is a useful astronomy Python package in general with many functionalities). It can for example convert between different coordinate systems (e.g. ICRS and Galactic) or between different notations. 
#     - coordinates can be written in the hms/dms notation, e.g. (RA, Dec) = (18h 36m 56.336s +38° 47′ 01.2802″), or in the deg/deg notation, e.g. (RA, Dec) = (279.23473333 38.78368894). 
# - The available magnitudes in Gaia are BP (blue), RP (red) and G (overall brightness).  
# - Some useful Gaia columns will be (but are not necessarily limited to): source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, parallax, pmra, pmdec (and the parallax and proper motion uncertainties)  
# - All the Gaia column names are explained here: https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html
# 
# -----

# **Some dereddening hints:**
# 
# m = reddened magnitude  
# m0 = dereddened magnitude
# 
# m0 = m - R\*E(B-V)
# 
# ===================
# 
# Values of R for Gaia magnitudes (from Casagrande et al. 2021, MNRAS, 507, 2684C):  
# 
# R_G  = 2.609 - 0.475\*C + 0.053\*C\*\*2  
# R_BP = 2.998 - 0.140\*C - 0.175\*C\*\*2 + 0.062\*C\*\*3  
# R_RP = 1.689 - 0.059\*C  
# 
# where C = (BP - RP)  
# 
# ===================
# 
# E(B-V) is the reddening along the line of sight for a given coordinate. The standard source for getting E(B-V) is the SFD map (Schlegel, Finkbeiner & Davis 1998, ApJ, 500, 525). It can be queried using the Python *dustmaps* package (https://dustmaps.readthedocs.io/en/latest/). 
# 
# ------

# # Solution

# In[ ]:





# In[ ]:





# In[ ]:




