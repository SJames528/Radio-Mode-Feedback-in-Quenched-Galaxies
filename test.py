from astropy.io import fits
from astropy.wcs import WCS
import astropy
import pandas as pd
import numpy as np
from library.utils import *

data_folder = "~/Cardiff/Summer Project/data/"

##read in matched catalogue
matched_catalogue = fits.open(data_folder + 'matched_lofar_sdss')

##read in mosaic
#mosaic = fits.open(data_folder+'/mosaics/p164+47-mosaic.fits')
#mosaic = fits.open(data_folder+'/mosaics/p176+60-mosaic.fits')
mosaic = fits.open(data_folder+'/mosaics/p169+55-mosaic.fits')

##matched catalogue to dataframe
mat_cat_df = pd.DataFrame(matched_catalogue[1].data)

##visualise all points in the DF which coincide with the selected mosaic
visualise(mat_cat_df, mosaic, s=50)
