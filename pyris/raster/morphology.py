from __future__ import division
import numpy as np
from skimage import morphology as mm

def CleanIslands( bw, size, conn=1 ):
    '''CleanIslands( bw, size, conn=1 ) - Remove Islands inside channel'''
    return ( mm.remove_small_objects( bw.astype( bool )==False, size, conn ) == False ).astype( bw.dtype )


def RemoveSmallObjects( bw, size, conn=1 ):
    '''RemoveSmallObjects( bw, size, conn=1 )'''
    return bw * mm.remove_small_objects( bw.astype(bool), size, conn ).astype( bw.dtype )

    
def Skeletonize( bw ):
    '''Combine Skeletonization with Isles Removal.'''
    skel, dist = mm.medial_axis( bw, return_distance=True )
    skel = mm.skeletonize( bw ) # Overwrite Medial Axis with Skeletonization (works better)
    return skel, dist
