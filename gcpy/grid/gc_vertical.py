import numpy as np
import scipy.sparse

class vert_grid:
    def __init__(self,AP=None,BP=None,p_sfc=1013.25):
        if (AP.size != BP.size) or (AP is None):
            # Throw error?
            print('Inconsistent vertical grid specification')
        self.AP = np.array(AP)
        self.BP = np.array(BP)
        self.p_sfc = p_sfc
    def p_edge(self):
        # Calculate pressure edges using eta coordinate
        return self.AP + self.BP * self.p_sfc
    def p_mid(self):
        p_edge = self.p_edge()
        return (p_edge[1:]+p_edge[:-1])/2.0

# Standard vertical grids
GEOS_72L_AP = np.array([ 0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
         1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
         4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
         7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
         1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
         1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
         2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
         2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
         1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
         7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01,
         4.017541e+01, 3.381001e+01, 2.836781e+01, 2.373041e+01,
         1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
         9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00,
         4.076571e+00, 3.276431e+00, 2.620211e+00, 2.084970e+00,
         1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
         6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01,
         2.113490e-01, 1.594950e-01, 1.197030e-01, 8.934502e-02,
         6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
         1.000000e-02 ])

GEOS_72L_BP = np.array([ 1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
         9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
         8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
         7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
         6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
         4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
         2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
         6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00 ])

GEOS_72L_grid = vert_grid(GEOS_72L_AP, GEOS_72L_BP)

# Reduced grid
GEOS_47L_AP = np.zeros(48)
GEOS_47L_BP = np.zeros(48)

# Fill in the values for the surface
GEOS_47L_AP[0] = GEOS_72L_AP[0]
GEOS_47L_BP[0] = GEOS_72L_BP[0]

# Build the GEOS 72-layer to 47-layer mapping matrix at the same time
xmat_i = np.zeros((72))
xmat_j = np.zeros((72))
xmat_s = np.zeros((72))

# Index here is the 1-indexed layer number
for i_lev in range(1,37):
    # Map from 1-indexing to 0-indexing
    x_lev = i_lev - 1
    # Sparse matrix for regridding
    # Below layer 37, it's 1:1
    xct = x_lev
    xmat_i[xct] = x_lev
    xmat_j[xct] = x_lev
    xmat_s[xct] = 1.0
    # Copy over the pressure edge for the top of the grid cell
    GEOS_47L_AP[i_lev] = GEOS_72L_AP[i_lev]
    GEOS_47L_BP[i_lev] = GEOS_72L_BP[i_lev]

# Now deal with the lumped layers
skip_size_vec = [2,4]
number_lumped = [4,7]

# Initialize
i_lev = 36
i_lev_72 = 36
for lump_seg in range(2):
    skip_size = skip_size_vec[lump_seg]
    # 1-indexed starting point in the 47-layer grid
    first_lev_47 = i_lev + 1
    first_lev_72 = i_lev_72 + 1
    
    # Loop over the coarse vertical levels (47-layer grid)
    for i_lev_offset in range(number_lumped[lump_seg]):
        # i_lev is the index for the current level on the 47-level grid
        i_lev = first_lev_47 + i_lev_offset
        # Map from 1-indexing to 0-indexing
        x_lev = i_lev - 1
        
        # Get the 1-indexed location of the last layer in the 72-layer grid 
        # which is below the start of the current lumping region
        i_lev_72_base = first_lev_72 + (i_lev_offset*skip_size) - 1
        
        # Get the 1-indexed location of the uppermost level in the 72-layer
        # grid which is within the target layer on the 47-layer grid
        i_lev_72 = i_lev_72_base + skip_size
        
        # Do the pressure edges first
        # These are the 0-indexed locations of the upper edge for the 
        # target layers in 47- and 72-layer grids
        GEOS_47L_AP[i_lev] = GEOS_72L_AP[i_lev_72]
        GEOS_47L_BP[i_lev] = GEOS_72L_BP[i_lev_72]
        
        # Get the total pressure delta across the layer on the lumped grid
        # We are within the fixed pressure levels so don't need to account
        # for variations in surface pressure
        dp_total = GEOS_47L_AP[i_lev-1] - GEOS_47L_AP[i_lev]
        
        # Now figure out the mapping
        for i_lev_offset_72 in range(skip_size):
            # Source layer in the 72 layer grid (0-indexed)
            x_lev_72 = i_lev_72_base + i_lev_offset_72
            xct = x_lev_72
            xmat_i[xct] = x_lev_72
            # Target in the 47 layer grid
            xmat_j[xct] = x_lev
            
            # Proportion of 72-layer grid cell, by pressure, within expanded layer
            xmat_s[xct] = (GEOS_72L_AP[x_lev_72] - GEOS_72L_AP[x_lev_72+1])/dp_total
    start_pt = i_lev

# Do last entry separately (no layer to go with it)
xmat_72to47 = scipy.sparse.coo_matrix((xmat_s,(xmat_i,xmat_j)),shape=(72,47))

GEOS_47L_grid = vert_grid(GEOS_47L_AP, GEOS_47L_BP)

# CAM 26-layer grid
CAM_26L_AP = np.flip(np.array([ 219.4067,   489.5209,   988.2418,   1805.201,        
                                2983.724,   4462.334,   6160.587,   7851.243,        
                                7731.271,   7590.131,   7424.086,   7228.744,        
                                6998.933,   6728.574,   6410.509,   6036.322,        
                                5596.111,   5078.225,   4468.96,    3752.191,        
                                2908.949,   2084.739,   1334.443,   708.499,         
                                252.136,    0.,         0.  ]),axis=0)*0.01
CAM_26L_BP = np.flip(np.array([ 0.,         0.,         0.,         0.,             
                                0.,         0.,         0.,         0.,             
                                0.01505309, 0.03276228, 0.05359622, 0.07810627,     
                                0.1069411,  0.14086370, 0.180772,   0.227722,       
                                0.2829562,  0.3479364,  0.4243822,  0.5143168,      
                                0.6201202,  0.7235355,  0.8176768,  0.8962153,      
                                0.9534761,  0.9851122,  1.        ]),axis=0)

CAM_26L_grid = vert_grid(CAM_26L_AP, CAM_26L_BP)
