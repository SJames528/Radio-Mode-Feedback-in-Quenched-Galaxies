from astropy.io import fits
from astropy.wcs import WCS
import astropy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def snap_data(df):
    return np.array([i for i in df["snapshot"]])

#takes a dataframe of points and returns an array of their pixel positions
def df_coords_to_pixels(df, coord_system, coord_names, new_names=["pix_x","pix_y"]):
    pix_coords = []
    for skypoint in df[coord_names].values:
        pix_coords.append(coord_system.world_to_pixel(astropy.coordinates.SkyCoord(skypoint[0],skypoint[1],unit="deg")))
    df.loc[:,new_names[0]] = [i[0] for i in pix_coords]
    df.loc[:,new_names[1]] = [i[1] for i in pix_coords] 
    return df

#takes an array and pads until a certain pixel is in the centre
# # Deprecated, as functionality has been incorporated into pixel_to_snapshot(). Keeping for other projects # #
def pad_square(array, central_pixel, pad_with=0):
    as_list = array.tolist()
    
    space_left = central_pixel[1]
    space_right = array.shape[1] - central_pixel[1] - 1
    space_up = central_pixel[0]
    space_down = array.shape[0] - central_pixel[0] - 1

    #horizontal padding
    left_padding = np.max(space_right-space_left,0) * [pad_with]
    right_padding = np.max(space_left-space_right,0) * [pad_with]
    for row_idx in range(len(as_list)):
        as_list[row_idx] = left_padding + as_list[row_idx] + right_padding

    #vertical padding
    up_padding = np.max(space_down-space_up,0) * [len(as_list[0]) * [pad_with]]
    down_padding = np.max(space_up-space_down,0) * [len(as_list[0]) * [pad_with]]
    as_list = up_padding + as_list + down_padding

    return np.array(as_list)

#takes a pixel position and returns a region of pixels around that central position
def pixel_to_snapshot(coord, mosaic, s, pad=True):
    if not (s % 2):
        raise Exception("Image size must be odd integer")
    nearest_pix = [int(np.round(a)) for a in coord]
    image_data = mosaic[0].data
    y_dim, x_dim = image_data.shape

    left_cutoff = int(nearest_pix[0]-((s-1)/2))
    right_cutoff = int(nearest_pix[0]+((s+1)/2))
    up_cutoff = int(nearest_pix[1]-((s-1)/2))
    down_cutoff = int(nearest_pix[1]+((s+1)/2))

    df_cut = np.array(image_data[max(up_cutoff,0):down_cutoff, max(left_cutoff,0):right_cutoff])

    if pad and df_cut.shape != (s,s):
        df_cut_list = df_cut.tolist()
        left_padding = max(0-left_cutoff,0) * [0]
        right_padding = max(right_cutoff-x_dim,0) * [0]
        up_padding = max(0-up_cutoff,0) * [s * [0]]
        down_padding = max(down_cutoff-y_dim,0) * [s * [0]]

        for row_idx in range(len(df_cut_list)):
            df_cut_list[row_idx] = left_padding + df_cut_list[row_idx] + right_padding
        df_cut_list = up_padding + df_cut_list + down_padding
        df_cut = np.array(df_cut_list)
    
    return df_cut

#calculate the RA/DEC range for the current mosaic
def mosaic_dim_limits(mos):
    head_ = mos[0].header
    ra_min = mos[0].header["CRVAL1"] + mos[0].header["CRPIX1"]*mos[0].header["CDELT1"]
    ra_max = mos[0].header["CRVAL1"] - mos[0].header["CRPIX1"]*mos[0].header["CDELT1"]
    dec_min = mos[0].header["CRVAL2"] - mos[0].header["CRPIX2"]*mos[0].header["CDELT2"]
    dec_max = mos[0].header["CRVAL2"] + mos[0].header["CRPIX2"]*mos[0].header["CDELT2"]
    return ra_min, ra_max, dec_min, dec_max

#produce plots of mosaic cutouts
def visualise(arr, fig_size=(10,10)):
    if isinstance(arr, pd.core.frame.DataFrame):
        arr = np.array(arr["snapshot"])
    if len(arr)==0 or not isinstance(arr, (list,np.ndarray)):
        raise Exception("Data empty, or not provided as list/array")
    while arr.shape[0]==1:
        arr = arr[0]
    if len(arr.shape)==2:
        plt.imshow(arr,cmap='grey')
        plt.axis('off')
    else:
        fig, ax = plt.subplots(int(0.5+len(arr)/2),2,figsize=fig_size)
        for index, cutout in enumerate(arr):
            subfig = ax[index//2,index%2] if len(arr)>2 else ax[index%2]
            subfig.imshow(cutout,cmap='gray')
            subfig.set_axis_off()
        if len(arr)%2:
            ax[-1,-1].axis('off')
    plt.show()

#plot the image data for all catalogue points in the current mosaic
def snapshots(df, mos, s=15, coord_names=["RA","DEC"], vis=False, vis_figsize=(10,10)):
    
    """Produces cutouts from a given mosaic centred at points given in a dataframe. The function will select only those objects within the mosaic's (assumed square) field, so the dataframe need not be pre-processed to contain only relevant sources. Size of cutout region can be specified. Coord names are assumed to be "RA" and "DEC", but these can be specified.

    :param df: dataframe OR list - the dataframe or list of tuples containing catalogue points to display. In the case of dataframe, columns for right ascension (RA) and declination (DEC) are required, but other contents are not required. If a list is entered, will be first converted to an appropriate dataframe.
    :param mos: mosaic object - usually loaded in from a fits mosaic with fits.open(file). Mosaic is assumed to be square (that is, all objects within the RA and DEC bounds will attemt to be plotted.
    :param s: int - size of the region (s x s square) to display, with the object position at the centre pixel, or above-left of centre in the case of even integer s value. Default value s=15.
    :param coord_names: [str, str] list - the column names for RA and DEC respectively in the provided dataframe. Default value coord_names=["RA","DEC"].
    """

    if isinstance(df, list):
        df = pd.DataFrame(df)
        df.columns = coord_names
        
    coord_sys = WCS(mos[0].header)
    ra_min, ra_max, dec_min, dec_max = mosaic_dim_limits(mos)
    points = df[(ra_min<=df[coord_names[0]]) & (df[coord_names[0]]<=ra_max) & (dec_min<=df[coord_names[1]]) & (df[coord_names[1]]<=dec_max)].copy()
    
    if len(points) == 0:
        raise Exception("Empty intersection with mosaic")
    points = df_coords_to_pixels(points, coord_sys, coord_names)
    
    ret = []
    for point in points[["pix_x", "pix_y"]].values:
        ret.append(pixel_to_snapshot(point, mos, s))
    points["snapshot"] = ret

    if vis: visualise(points, figsize = vis_figsize)
    return points

def stack(cutout_arr, weight_method):
    
    return
