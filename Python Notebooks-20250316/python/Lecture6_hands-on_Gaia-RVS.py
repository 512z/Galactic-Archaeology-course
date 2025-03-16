#!/usr/bin/env python
# coding: utf-8

# # Imports

# In[30]:


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io.votable import parse
import glob

from matplotlib import rc
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 20}	
rc('font', **font)

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap            # umap-learn package


# ----------------------------------------
# **Information about the data**
# 
# Download from here: https://drive.google.com/file/d/1xuh69L_fO6G-h2rPsXdT0hs-f9bepVuZ/view?usp=sharing
# 
# The rvsspec folder contains the relevant data:
# 
# - individual RVS *.xml files 
# - specpar.csv (parameters from the Gaia team)
# - cnn_gaia_rvs_catalog_guiglion_et_al_2023_v2.csv (CNN parameters from Guiglion et al. 2023)
# - rave_linelist.csv
# 
# To obtain the XML files, I ran the query below on the Gaia archive
# and used the DataLink service to download the associated spectra.
# You could also adapt this to obtain additional parameters from Gaia.
# 
# ```
# SELECT TOP 5000 source_id
#    FROM gaiadr3.gaia_source
#   WHERE has_rvs = 'True'
#     AND random_index < 7000000
#     AND rvs_spec_sig_to_noise > 50 
# ```
# ----------------------------------------
# 

# In[2]:


### GET RAVE LINELIST 

rave = pd.read_csv('rvsspec/rave_linelist.csv', delim_whitespace=True)

## convert air to vacuum wavelengths
sigma2 = (10000./rave.wl)**2.
fact = 1.0 + (6.4328e-5 + 2.94981e-2 / (146. - sigma2) +
               2.5540e-4 / (41. - sigma2))
wl_vac = rave.wl*fact

rave = rave.assign(wl_vac = wl_vac)

## get lines for different elements
rave_mg = rave[rave.el == "MgI"]
rave_si = rave[rave.el == "SiI"]
rave_ti = rave[rave.el == "TiI"]
rave_ni = rave[rave.el == "NiI"]
rave_al = rave[rave.el == "AlI"]
rave_fe = rave[rave.el == "FeI"]


# In[3]:


### DEFINE PLOTTING FUNCTION 

def plot_spectra(spectra, wl_range, txt=None):

    nspectra = min(10, len(spectra))
    plt.figure(figsize=(14, 1.5 * nspectra))
    ax = plt.subplot(111)
    for i in range(nspectra):
        ax.plot(wl_range, spectra.iloc[i][:2401] + i * 0.7, 'k-')
    if txt is not None:
        for i, t in enumerate(txt):
            ax.text(wl_range[20], i * 0.7 + 1.1, t)
    ax.set_ylim(0, 0.7 + 0.7 * nspectra)
    ax.set_xlim(wl_range[0], wl_range[-1])
    y0, y1 = ax.get_ylim()
    ax.set_xlabel(r'Wavelength $\mathrm{[\AA]}$')
    ax.set_ylabel(r'Flux')

    ### PLOT RAVE LINELIST
    for wlist, kleur,lab in zip([rave_mg.wl_vac, rave_si.wl_vac, rave_ti.wl_vac, rave_ni.wl_vac, rave_al.wl_vac, rave_fe.wl_vac],
                            ['green', 'orange', 'cyan', 'magenta', 'blue', 'black'],
                            ['Mg', 'Si', 'Ti', 'Ni', 'Al', 'Fe']):
        for wli in wlist:
            ax.axvline(wli, linewidth=1, linestyle='--', color=kleur)
        
        ax.axvline(0.0, linewidth=1, linestyle='--', color=kleur, label=lab)

    for wl in (8500.35, 8544.44, 8664.52):
        ax.plot([wl, wl], [y0, y1], color='grey', alpha=0.3)
    ax.axvline(0.0, linewidth=1, linestyle='-', color='grey', alpha=0.3, label='Ca')

    ax.legend(ncol=7, fontsize=14, loc='lower center', framealpha=0.8)

    plt.show()


# # Data

# In[ ]:


par0 = pd.read_csv("rvsspec/cnn_gaia_rvs_catalog_guiglion_et_al_2023_v2.csv")
par = par0[(par0.flag_boundary == 0) & (par0.efeh < 0.25)]    # select stars with good parameters

### you can come back here after going through the notebook, try what happens without the "efeh" cut

list(par)


# In[ ]:


### could also try the Gaia team parameters (but first try the CNN ones above):
# par = pd.read_csv("rvsspec/specpar.csv")
# par = par.rename(columns={'teff_gspspec':'teff', 'logg_gspspec':'logg','mh_gspspec':'mh','alphafe_gspspec':'alpham'})


# In[ ]:


files = glob.glob("rvsspec/*.xml")

spectra = pd.DataFrame()
snrs = []
names = []
teffs = []
loggs = []
mhs = []
ams = []

i = 0
for file in files: 
    votable = parse(file)

    name = file.split("/")[-1].split(".")[0].split(" ")[-1]
    

    if i == 0:
        wl_range = votable.get_first_table().array['wavelength']*10
    fl = votable.get_first_table().array['flux']
    snr = np.median(votable.get_first_table().array['flux']/votable.get_first_table().array['flux_error'])

    
    if (snr > 10) & (int(name) in par.source_id.values):
        spectra = pd.concat([spectra, pd.DataFrame(fl).T], ignore_index=True)
        snrs.append(snr)
        names.append(name)

        teffs.append(par[par.source_id == int(name)].teff.values[0])
        loggs.append(par[par.source_id == int(name)].logg.values[0])
        mhs.append(par[par.source_id == int(name)].mh.values[0])
        ams.append(par[par.source_id == int(name)].alpham.values[0])
    i+=1

spectra = spectra.fillna(1)

spectra_std = StandardScaler().fit_transform(spectra)

spectra.shape


# # Plotting data

# ## Spectra sample

# In[ ]:


plot_spectra(spectra.loc[spectra.sample(n=10).index.values], wl_range)


# ## Kiel diagram

# In[ ]:





# ## Metallicity distribution

# In[ ]:





# ## alpha/Fe vs. Fe

# In[ ]:





# ## Explore spectra with different parameters

# In[ ]:





# # Dimensionality reduction

# **Note:** *reduced X & Y  = algorithm_results[:, 0] & algorithm_results[:, 1]*

# ## PCA (Principal Component Analysis)

# In[13]:


pca = PCA(n_components=2)
pca_results = pca.fit_transform(spectra_std)


# In[ ]:


## Plot results and colour-code by different parameters



# ## t-SNE (t-distributed Stochastic Neighbor Embedding)

# In[ ]:


tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(spectra_std)


# In[ ]:


## Plot results and colour-code by different parameters



# ## UMAP (Uniform Manifold Approximation and Projection)

# In[ ]:


reducer = umap.UMAP()
umap_results = reducer.fit_transform(spectra_std)


# In[ ]:


## Plot results and colour-code by different parameters


