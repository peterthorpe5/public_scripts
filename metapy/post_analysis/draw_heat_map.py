#!/usr/bin/env python3
from scipy.io import netcdf
#https://github.com/widdowquinn/Teaching-Data-Visualisation/blob/master/exercises/colormaps_surfaces_netcdf/colormaps_surfaces_netcdf.ipynb
import warnings
warnings.filterwarnings('ignore')








in order to generate the data.fram.

use pandas
pandas.dataframe
look at this function

https://github.com/widdowquinn/pyani/blob/master/pyani/anib.py
def process_blast



to darw the datfram.

https://github.com/widdowquinn/pyani/blob/master/pyani/pyani_graphics.py


def heatmap_seaborn(df, 
















#############################################################################

#to run the script       

usage = """usage :

script to graphically represent clustering output.


python draw_heatmap.py -i clustering_file -o summarise_clusters.out

requires:
use pip install ...
 seaborn
 matplotlib
 numpy

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")
parser.add_option("--heatmap", dest="heatmap", default=True,
                  help="draw a heat map of the species clustering")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
heatmap = options.heatmap
out_file = options.out_file

#run the program

parse_tab_file_get_clusters(in_file, out_file)

#if heatmap:
    #from draw_heat_map import *

print "done"
=======
<<<<<<< HEAD
#!/usr/bin/env python
from scipy.io import netcdf
#https://github.com/widdowquinn/Teaching-Data-Visualisation/blob/master/exercises/colormaps_surfaces_netcdf/colormaps_surfaces_netcdf.ipynb
import warnings
warnings.filterwarnings('ignore')








in order to generate the data.fram.

use pandas
pandas.dataframe
look at this function

https://github.com/widdowquinn/pyani/blob/master/pyani/anib.py
def process_blast



to darw the datfram.

https://github.com/widdowquinn/pyani/blob/master/pyani/pyani_graphics.py


def heatmap_seaborn(df, 
















#############################################################################

#to run the script       

usage = """usage :

script to graphically represent clustering output.


python draw_heatmap.py -i clustering_file -o summarise_clusters.out

requires:
use pip install ...
 seaborn
 matplotlib
 numpy

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")
parser.add_option("--heatmap", dest="heatmap", default=True,
                  help="draw a heat map of the species clustering")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
heatmap = options.heatmap
out_file = options.out_file

#run the program

parse_tab_file_get_clusters(in_file, out_file)

#if heatmap:
    #from draw_heat_map import *

print "done"
>>>>>>> b5a3cd9f461b5e65b0b0748712c44268bd2236a1
