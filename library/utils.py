from astropy.io import fits
from astropy.wcs import WCS
import astropy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def snap_data(df):
    return np.array([i for i in df["snapshot"]])

#takes a dataframe of points and adds new columns with their pixel positions
def df_skycoords_to_pixels(df, mos, coord_system):
    mos_pointer = mos
    pix_coords = [coord_system.world_to_pixel(skypoint) for skypoint in df["skycoord"].values]
    df.loc[:,"mosaic_coord"] = [(mos_pointer, coords) for coords in pix_coords]
    return df

def df_add_skycoords(df, coord_names):
    df.loc[:,"skycoord"] = [astropy.coordinates.SkyCoord(skypoint[0],skypoint[1],unit="deg") for skypoint in df[coord_names].values]
    return df

#takes an array and pads until a certain pixel is in the centre
# # DEPRECATED, as functionality has been incorporated into pixel_to_snapshot(). Keeping for other projects # #
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
def pixel_to_snapshot(mos_coord, s, pad=True):
    mosaic, coord = mos_coord
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

    df_cut = np.nan_to_num(df_cut, nan=0)

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

#produce plots of mosaic cutouts
def visualise(arr, title="", fig_size=(10,10)):
    if isinstance(arr, pd.core.frame.DataFrame):
        arr = snap_data(arr)
    if len(arr)==0 or not isinstance(arr, (list,np.ndarray,pd.core.series.Series)):
        raise Exception("Data empty, or not provided as list/array")
    while arr.shape[0]==1:
        arr = arr[0]
    if len(arr.shape)==2:
        arr = arr * 1000
        plt.imshow(arr,cmap='grey')
        plt.axis('off')
        plt.colorbar(label="mJy/beam")
        plt.title(title)
    else:
        fig, ax = plt.subplots(int(0.5+len(arr)/2),2,figsize=fig_size)
        for index, cutout in enumerate(arr):
            subfig = ax[index//2,index%2] if len(arr)>2 else ax[index%2]
            subfig.imshow(cutout,cmap='gray')
            subfig.set_axis_off()
        if len(arr)%2:
            ax[-1,-1].axis('off')
        fig.suptitle(title)
    plt.show()

#plot the image data for all catalogue points in the current mosaic
def snapshots(df, mosaics, s=25, coord_names=["RA","DEC"], vis=False, vis_figsize=(10,10)):
    
    """Produces cutouts from a given mosaic centred at points given in a dataframe. The function will select only those objects within the mosaic's (assumed square) field, so the dataframe need not be pre-processed to contain only relevant sources. Size of cutout region can be specified. Coord names are assumed to be "RA" and "DEC", but these can be specified.

    :param df: dataframe OR list - the dataframe or list of tuples containing catalogue points to display. In the case of dataframe, columns for right ascension (RA) and declination (DEC) are required, but other contents are not required. If a list is entered, will be first converted to an appropriate dataframe.
    :param mos: mosaic object or list of objects - usually loaded in from a fits mosaic with fits.open(file). Mosaic is assumed to be square (that is, all objects within the RA and DEC bounds will attemt to be plotted.
    :param s: int - size of the region (s x s square) to display, with the object position at the centre pixel, or above-left of centre in the case of even integer s value. Default value s=15.
    :param coord_names: [str, str] list - the column names for RA and DEC respectively in the provided dataframe. Default value coord_names=["RA","DEC"].
    """

    if isinstance(df, list):
        df = pd.DataFrame(df); df.columns = coord_names
    if isinstance(mosaics, astropy.io.fits.hdu.hdulist.HDUList):
        mosaics = [mosaics]

    df = df_add_skycoords(df, coord_names)
    ret = pd.DataFrame([])
    for mos in mosaics:
        coord_sys = WCS(mos[0].header)
        if not "mosaic_coord" in df.columns:
            points = df_skycoords_to_pixels(df.copy(), mos, coord_sys)
            points = points[((0<=np.array([i[1][0] for i in points["mosaic_coord"]]))&(np.array([i[1][0] for i in points["mosaic_coord"]])<coord_sys.pixel_shape[0]))&((0<=np.array([i[1][1] for i in points["mosaic_coord"]]))&(np.array([i[1][1] for i in points["mosaic_coord"]])<coord_sys.pixel_shape[1]))]
        else:
            points = df.copy()

        if len(points):
            points["snapshot"] = [pixel_to_snapshot(mos_coord, s) for mos_coord in points["mosaic_coord"].values]
        
            ret = pd.concat([ret,points])

    if not len(ret):
        raise Exception("Empty intersection with mosaic")
    ret = ret[np.array([not not np.count_nonzero(i) for i in ret["snapshot"]])]

    if vis: visualise(ret, figsize = vis_figsize)
    return ret

def stack(df, weight_method="IVW"):
    avg_snap = np.zeros(snap_data(df)[0].shape)
    inv_var_sum = 0
    for mat in snap_data(df):
        var = np.std(mat)**2
        if not var:
            if avg_snap.shape == (1,1):
                avg_snap += mat/len(df)
            continue
        weighted_cutout = mat/var
        avg_snap += weighted_cutout
        inv_var_sum += 1/var
    return avg_snap/inv_var_sum if inv_var_sum else avg_snap

def skycoord_in_survey_region(point):
    ra, dec = point
    #is coord in ra/dec bounds?
    if not (162.5<ra<237.5 and 46<dec<61.25): return False
    #is coord within survey shape?
    return ((170<ra<225 and 46<dec<61.25) or
            (225<ra<237.5 and 46<dec<56.5) or
            (225<ra<229.5 and 56.5<dec<58.8) or
            (162.5<ra<170 and 46<dec<48.5) or
            (167.5<ra<170 and 48.5<dec<58.75) or
            (ra<167.5 and dec>48.5 and (9*dec<(22*ra)-3199)))

def random_stack(rng, mos, len_, s=25):
    random_coord_dict = {'RA': [], 'DEC':[]}
    while len(random_coord_dict['RA'])<len_:
        new_rand_point = rng.random(2)*(75,15.25) + (162.5,46)
        if skycoord_in_survey_region(new_rand_point):
            random_coord_dict['RA'].append(new_rand_point[0])
            random_coord_dict['DEC'].append(new_rand_point[1])
    random_coord_df = snapshots(pd.DataFrame(random_coord_dict), mos, s)
    
    random_coord_df = random_coord_df[~random_coord_df.index.duplicated()]
    return stack(random_coord_df)

def prog_bar(i,n,bar_len=25):
    progress = int(np.floor((i/n)*bar_len))
    print("|"+"#"*progress+"_"*(bar_len-progress)+f"| Progress={i}/{n}", end='\r') if i!=n else print("|"+"#"*progress+"_"*(bar_len-progress)+f"| Progress={i}/{n}")

def verbose_iterate_to_array(n=100,bar_len=25):
    def decorator(func):
        list_ = []
        def wrapper(*args, **kwargs):
            for i in range(1,n+1):
                list_.append(func(*args, **kwargs))
                prog_bar(i,n,bar_len)
            return np.array(list_)
        return wrapper
    return decorator
