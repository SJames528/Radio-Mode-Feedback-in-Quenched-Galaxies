from astropy.io import fits
from astropy.wcs import WCS
import astropy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#takes a dataframe of points and returns an array of their pixel positions
def df_coords_to_pixels(df, coord_system):
    output = []
    for skypoint in df[["RA_1","DEC_1"]].values:
        output.append(coord_system.world_to_pixel(astropy.coordinates.SkyCoord(skypoint[0],skypoint[1],unit="deg")))
    return output

#takes a pixel position and returns a region of pixels around that central position (
def pixel_to_snapshot(coord, mosaic, size):
    nearest_pix = [int(np.round(a)) for a in coord]
    image_data = mosaic[0].data
    y_dim, x_dim = image_data.shape
    y_cut = image_data[int(nearest_pix[1]+1-(size/2)):int(nearest_pix[1]+1+(size/2))]
    x_cut = np.array([row[int(nearest_pix[0]+1-(size/2)):int(nearest_pix[0]+1+(size/2))] for row in y_cut])
    return x_cut

#calculate the RA/DEC range for the current mosaic
def mosaic_dim_limits(mos):
    head_ = mos[0].header
    ra_min = mos[0].header["CRVAL1"] + mos[0].header["CRPIX1"]*mos[0].header["CDELT1"]
    ra_max = mos[0].header["CRVAL1"] - mos[0].header["CRPIX1"]*mos[0].header["CDELT1"]
    dec_min = mos[0].header["CRVAL2"] - mos[0].header["CRPIX2"]*mos[0].header["CDELT2"]
    dec_max = mos[0].header["CRVAL2"] + mos[0].header["CRPIX2"]*mos[0].header["CDELT2"]
    return ra_min, ra_max, dec_min, dec_max

#plot the image data for all catalogue points in the current mosaic
def visualise(df, mos, s=10, fig_size=(10,10)):
    coord_sys = WCS(mos[0].header)
    ra_min, ra_max, dec_min, dec_max = mosaic_dim_limits(mos)
    pixs = df[(ra_min<=df["RA_1"]) & (df["RA_1"]<=ra_max) & (dec_min<=df["DEC_1"]) & (df["DEC_1"]<=dec_max)]
    points = df_coords_to_pixels(pixs, coord_sys)
    fig, ax = plt.subplots(int(0.5+len(pixs)/2),2,figsize=fig_size)
    for index, point in enumerate(points):
        ax[index//2,index%2].imshow(pixel_to_snapshot(point, mos, s),cmap='gray')
        ax[index//2,index%2].set_axis_off()
    if len(pixs)%2: ax[-1,-1].axis('off')
    plt.show()
        
