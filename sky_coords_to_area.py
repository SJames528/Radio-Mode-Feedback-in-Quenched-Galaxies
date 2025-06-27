import numpy as np
import matplotlib.pyplot as plt
def area_calc(RA_range,DEC_range):
     RA_range_rad = [x*np.pi/180 for x in RA_range]
     DEC_range_rad = [x*np.pi/180 for x in DEC_range]
     area = (RA_range_rad[1]-RA_range_rad[0])*(np.sin(DEC_range_rad[1])-np.sin(DEC_range_rad[0]))
     area_deg_2 = area * (180/np.pi)**2
     return area_deg_2

coord_ranges = [([170,225],[46,61.25]),
                ([225,237.5],[46,56.5]),
                ([225,229.5],[56.5,58.8]),
                ([165.2,170],[46,48.5]),
                ([167.5,170],[48.5,58.75])]

tot_area = sum([area_calc(area[0],area[1]) for area in coord_ranges])


match_record = [(1.0,175),
                (2.0,451),
                (3.0,635),
                (4.0,736),
                (5.0,823),
                (6.0,891),
                (7.0,927),
                (8.0,958),
                (9.0,993),
                (10.0,1022)]

plt.scatter(x = [coord[0] for coord in match_record], y = [coord[1] for coord in match_record])
plt.xlabel("Max error for match (arcseconds)")
plt.ylabel("Number of matches")
plt.show()
