from astropy.io import fits
from astropy.wcs import WCS
import astropy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#takes a dataframe of points and returns an array of their pixel positions
def df_coords_to_pixels(df, coord_system, coord_names):
    output = []
    for skypoint in df[coord_names].values:
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
def visualise(df, mos, s=5, fig_size=(10,10), coord_names=["RA","DEC"], ret_snapshot_info=False):
    
    """Produce images from a given mosaic centred at points in a dataframe. The function will select only those objects within the mosaic's (assuemd square) field, so the dataframe need not be pre-processed to contain only relevant sources. Size of regions and fig sizes can be specified. Coord names are assumed to be "RA" and "DEC", but these can be specified. Can optionally return the image data in an array rather than plotting.

    :param df: dataframe - the dataframe containing catalogue points to display. Columns for right ascension (RA) and declination (DEC) are required, but other contents are not required.
    :param mos: mosaic object - usually loaded in from a fits mosaic with fits.open(file). Mosaic is assumed to be square (that is, all objects within the RA and DEC bounds will attemt to be plotted.
    :param s: int - size of the region (s x s square) to display, with the object position at the centre pixel, or above-left of centre in the case of even integer s value. Default value s=5.
    :param fig_size: (int, int) tuple - the size of the subplots (matplotlib.pyplot) for the returned images. Default value fig_size=(10, 10).
    :param coord_names: [str, str] list - the column names for RA and DEC respectively in the provided dataframe. Default value coord_names=["RA","DEC"].
    :param ret_snapshot_info: bool - if True, function will forgo displaying the image data, and instead return the data in a numpy array. Default value ret_snapshot_info=False.
    """
    
    coord_sys = WCS(mos[0].header)
    ra_min, ra_max, dec_min, dec_max = mosaic_dim_limits(mos)
    pixs = df[(ra_min<=df[coord_names[0]]) & (df[coord_names[0]]<=ra_max) & (dec_min<=df[coord_names[1]]) & (df[coord_names[1]]<=dec_max)]
    points = df_coords_to_pixels(pixs, coord_sys, coord_names)
    fig, ax = plt.subplots(int(0.5+len(pixs)/2),2,figsize=fig_size)
    ret = []
    for index, point in enumerate(points):
        ax[index//2,index%2].imshow(pixel_to_snapshot(point, mos, s),cmap='gray')
        if ret_snapshot_info:
            ret.append(pixel_to_snapshot(point, mos, s))
        ax[index//2,index%2].set_axis_off()
    if len(pixs)%2: ax[-1,-1].axis('off')
    if ret_snapshot_info:
        return np.array(ret)
    plt.show()
        
