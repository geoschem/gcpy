'''get gamap's WhGrYlRd color scheme from file'''
import numpy as np
from matplotlib.colors import ListedColormap
import os

current_dir = os.path.dirname(__file__)

rgb_WhGrYlRd = np.genfromtxt(current_dir+'/colormap_data/WhGrYlRd.txt',
                             delimiter=' ')
WhGrYlRd = ListedColormap(rgb_WhGrYlRd/255.0)
