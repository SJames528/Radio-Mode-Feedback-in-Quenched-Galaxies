from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import astropy
import pandas as pd
import numpy as np
from library.utils import *

data_folder = "~/Cardiff/Summer Project/data/"

##read in matched catalogue
matched_catalogue = fits.open(data_folder + 'matched_lofar_sdss')
mat_cat_df = Table(matched_catalogue[1].data).to_pandas()

##read in mosaic
#mosaic = fits.open(data_folder+'/mosaics/p164+47-mosaic.fits')
#mosaic = fits.open(data_folder+'/mosaics/p176+60-mosaic.fits')
mosaic = fits.open(data_folder+'/mosaics/p169+55-mosaic.fits')

##visualise all points in the DF which lie within the selected mosaic
#visualise(mat_cat_df, mosaic, s=50, coord_names=["RA_1","DEC_1"])

##Read in entire quenched catalogue
quenched_cat = fits.open(data_folder + 'SDSS_quenched_in_skyarea')
quenched_df = Table(quenched_cat[1].data).to_pandas().drop(columns = ["Source_S"])
quenched_df["combined_coord"] = [item for item in zip(quenched_df["RA"],quenched_df["DEC"])]

##Keep record of matched radio objects
SDSS_matched_coords = [item for item in zip(mat_cat_df["RA_2"],mat_cat_df["DEC_2"])]

##Retain only quenched sources which are not radio sources
quenched_no_radio_df = quenched_df[~quenched_df["combined_coord"].isin(SDSS_matched_coords)]

snap_info_quiet = visualise(quenched_no_radio_df[:1000], mosaic, s=25, ret_snap_info=True)
snap_info_loud = visualise(mat_cat_df, mosaic, s=25, coord_names=["RA_1","DEC_1"], ret_snap_info=True)

#plot histograms of snapshot pixel values
if False:
    fig, ax = plt.subplots(2,1,figsize=(10,10))
    ax[0].hist(snap_info_quiet[0].flatten(), bins=np.linspace(np.min(snap_info_loud[0]), np.max(snap_info_loud[0]), 50))
    ax[0].set_title("Histogram of pixel values for radio-quiet quenched galaxy"); ax[0].set_xlabel("Pixel value"); ax[0].set_ylabel("Count")
    ax[1].hist(snap_info_loud[0].flatten(), bins=np.linspace(np.min(snap_info_loud[0]), np.max(snap_info_loud[0]), 50))
    ax[1].set_title("Histogram of pixel values for radio-loud quenched galaxy"); ax[1].set_xlabel("Pixel value"); ax[1].set_ylabel("Count")
    plt.show()

for index, mat in enumerate(snap_info_quiet):
    if not index:
        avg_snap = mat
    else:
        avg_snap += mat
plt.imshow(avg_snap); plt.show()
