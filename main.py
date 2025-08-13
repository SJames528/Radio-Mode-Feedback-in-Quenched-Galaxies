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
mosaic_small = fits.open(data_folder+'/mosaics/p164+47-mosaic.fits')
#mosaic = fits.open(data_folder+'/mosaics/p176+60-mosaic.fits')
#mosaic = fits.open(data_folder+'/mosaics/p169+55-mosaic.fits')
mosaic = fits.open(data_folder+'/mosaics/mosaic-i_1.fits')
mosaic2 = fits.open(data_folder+'/mosaics/mosaic-i_2.fits')
mosaic_lr = fits.open(data_folder+'/mosaics/mosaic-lr.fits')

mosaic_active = [mosaic, mosaic2]
#mosaic_active = mosaic_lr
global_size = 25

##visualise all points in the DF which lie within the selected mosaic
#vis_data = snapshots(mat_cat_df, mosaic, s=15, coord_names=["RA_1","DEC_1"])
#visualise(vis_data)

##Read in entire quenched catalogue
quenched_cat = fits.open(data_folder + 'SDSS_quenched_in_skyarea')
quenched_df = Table(quenched_cat[1].data).to_pandas().drop(columns = ["Source_S","Z","Z_ERR","LGM_TOT_P50","SFR_TOT_P50","SSFR"])
quenched_df["combined_coord"] = [item for item in zip(quenched_df["RA"],quenched_df["DEC"])]

##Keep record of matched radio objects
SDSS_matched_coords = [item for item in zip(mat_cat_df["RA_2"],mat_cat_df["DEC_2"])]

##Retain only quenched sources which are not radio sources
quenched_no_radio_df = quenched_df[~quenched_df["combined_coord"].isin(SDSS_matched_coords)]
quenched_no_radio_df = quenched_no_radio_df.drop(columns=["combined_coord"])

snap_info_quiet = snapshots(quenched_no_radio_df, mosaic_active, s=global_size)
snap_info_loud = snapshots(mat_cat_df, mosaic_active, s=global_size, coord_names=["RA_1","DEC_1"])

##Mean/std per cutout
if False:
    std_list = np.array([np.mean(cutout) for cutout in snap_data(snap_info_quiet)])
    plt.hist(std_list, bins=bins, log=True); plt.title("Histogram of mean cutout values over entire catalogue")
    plt.xlabel("Mean flux (Jy/beam)"); plt.ylabel("Logged count")
    plt.show()


##See what stacking looks like on radio-loud galaxies, as an example of a positive result
if False:
    stack_loud = stack(snap_info_loud)
    visualise(stack_loud, title="stacking analysis on radio-loud data")

##Plot some histograms of typical snapshot pixel values
if False:
    fig, ax = plt.subplots(2,1,figsize=(10,10))
    ax[0].hist(snap_data(snap_info_quiet[:100]).flatten(), bins=np.linspace(np.min(snap_data(snap_info_loud)), np.max(snap_data(snap_info_loud)), 50), log=True)
    ax[0].set_title("Histogram of pixel values for first 100 radio-quiet quenched galaxies"); ax[0].set_xlabel("Pixel value"); ax[0].set_ylabel("Count")
    ax[1].hist(snap_data(snap_info_loud).flatten(), bins=np.linspace(np.min(snap_data(snap_info_loud)), np.max(snap_data(snap_info_loud)), 50), log=True)
    ax[1].set_title("Histogram of pixel values for radio-loud quenched galaxies"); ax[1].set_xlabel("Pixel value"); ax[1].set_ylabel("Count")
    plt.show()

##1 - Duplicated points (appear on both mosaics). TODO: decide whether to delete duplicates after, or put a condition in the snapshot() function which checks if that index already exists with a snapshot
duplicated_points = snap_info_quiet[snap_info_quiet.index.duplicated(keep=False)]
snap_info_quiet = snap_info_quiet[~snap_info_quiet.index.duplicated()]
if False:
    print(f'Are these duplicated points present exactly twice? {np.all(duplicated_points.index.value_counts()==2)}')
    print(f'Number of duplicated points: {int(len(duplicated_points)/2)}')
    plt.scatter(duplicated_points["RA"], duplicated_points["DEC"]); plt.title("Points appearing on both mosaics"); plt.show()

    duplicated_differences = pd.DataFrame([])
    duplicated_differences.index = list(set(duplicated_points.index))
    duplicated_differences["diffs"] = [np.linalg.norm(a-b) for a, b in [duplicated_points.loc[i]["snapshot"] for i in duplicated_differences.index]]

    plt.scatter([i for i in range(0,len(duplicated_differences))],duplicated_differences["diffs"]); plt.title("Euclidean distance of separation over the set of duplicated cutouts"); plt.show()
    visualise(duplicated_points.loc[duplicated_differences[duplicated_differences["diffs"]>2].index], title="Duplicated cutouts which differ the most from eachother")
    # -> Appears that there's just a rotational difference? Particuarly clear for brighter cutouts. Noisy cutouts appear less consistent, noise in both is independent?


##2 - Some snapshots contain bright radio sources. These should in theory be dealt with by the variance weighting when stacking, but better safe than sorry.
if False:
    max
    cutouts_with_large_vals = snap_info_quiet[[0.03<np.max(i)<0.04 for i in snap_data(snap_info_quiet)]]
    mid_vals = snap_info_quiet[[0.03<np.max(i)<0.04 for i in snap_data(snap_info_quiet)]]
    visualise(mid_vals[:10], title="")
if False:
    problematic_snaps = snap_info_quiet[[np.max(i) > 0.03 for i in snap_data(snap_info_quiet)]]
    snap_info_quiet = snap_info_quiet.drop(problematic_snaps.index)

##3 - Some snapshots have naturally higher baseline signal (due to being too close to bright sources or other noise artifacts). Need some way to subtract the baseline, or remove those with high median values
if False:
    low_median_unproblematic = snap_info_quiet[[np.median(i)<0.01 for i in snap_data(snap_info_quiet)]]

##Unused data - fixed with updates to utils.py (no longer pre-filtering, using all skycoords which have valid pixel positions)
if False:
    unused = quenched_no_radio_df[[index not in snap_info_quiet.index for index in quenched_no_radio_df.index]]
    plt.scatter(unused["RA"], unused["DEC"]); plt.show()

##Stacking analysis on filtered dataset
stack_data = stack(snap_info_quiet)
visualise(stack_data, title="Stacked cutout on processed dataset")


##Try a random stack to compare - need to rework to use random RA/DEC not pixel position

random_coord_dict = {'RA': [], 'DEC': []}
rng = np.random.default_rng(123)
while len(random_coord_dict['RA'])<len(snap_info_quiet):
    new_rand_point = rng.random(2)*(75,15.25) + (162.5,46)
    if skycoord_in_survey_region(new_rand_point):
        random_coord_dict['RA'].append(new_rand_point[0])
        random_coord_dict['DEC'].append(new_rand_point[1])
random_coord_df = pd.DataFrame(random_coord_dict)

random_coord_snaps = snapshots(random_coord_df, mosaic_active, s=global_size)
#Visual of duplicated region:
if True:
    dup_points = random_coord_snaps[random_coord_snaps.index.duplicated(keep=False)]
    plt.scatter(random_coord_snaps["RA"], random_coord_snaps["DEC"], s=2)
    plt.scatter(dup_points["RA"], dup_points["DEC"], s=2)
    plt.title("Mosaic overlap region - points in the intersection\nhave two distinct data cutouts available")
    plt.show()

random_coord_snaps = random_coord_snaps[~random_coord_snaps.index.duplicated()]
random_coord_stack = stack(random_coord_snaps)
visualise(random_coord_stack, title="Stacked random points in the mosaics")

##100 stacked cutouts of random populations:
@verbose_iterate_to_array(n=1000)
def random_iterated_stack(*args,**kwargs):
    return random_stack(*args,**kwargs)

rng = np.random.default_rng()
if True:
    random_snapshots = random_iterated_stack(rng, mosaic_active, len_=len(snap_info_quiet), s=1)
    with open("100_random_n_points.npy", 'wb') as f:
        np.save(f, random_snapshots)
else:  
    with open("100_random_n_stacks.npy", 'rb') as f:
        random_snapshots = np.load(f)

##Histograms of pixel values:
if False:
    plt.hist(stack_data.flatten()*1000, bins=np.linspace(-0.1,0.4,num=100))
    plt.hist(random_coord_stack.flatten()*1000, bins=np.linspace(-0.1,0.4,num=100))
    plt.title("Historgram of pixel values for stacked galaxies vs random cutouts")
    plt.show()
#look into why most snaps have positive mean but stacked the mean is negative
if False:
    plt.hist([random_coord_stack.flatten()*100,snap_data(random_coord_snaps)[:100].flatten()], log=True, alpha=0.2, bins=np.linspace(-0.02,0.05,num=100))
    plt.show()

##Angular separation analysis

if True:
    snap_info_quiet = snap_info_quiet.sort_index()
    ang_sep_l = []
    for i in snap_info_quiet.index:
        nearby = snap_info_quiet[(np.abs(snap_info_quiet["RA"]-snap_info_quiet["RA"][i])<1)&(np.abs(snap_info_quiet["DEC"]-snap_info_quiet["DEC"][i])<1)].drop([i])
        ang_sep = np.min(np.array([snap_info_quiet.loc[i]["skycoord"].separation(nearby_point).arcsecond for nearby_point in nearby["skycoord"]]))
        ang_sep_l.append(ang_sep)
    snap_info_quiet["closest_angsep_arc"] = ang_sep_l

    #!! there are 788 instances of exactly matching objects
    non_overlapping_points = snap_info_quiet[snap_info_quiet["closest_angsep_arc"]!=0]
    #!! there are 82 instances where the separation is less than 15" FWHM


    snap_info_quiet = snap_info_quiet[snap_info_quiet["closest_angsep_arc"]>3]

    #!! there are 355 instances where the separation is less than the size of a cutout


