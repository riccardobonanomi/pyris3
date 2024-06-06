from __future__ import division
import os, sys
import numpy as np
from skimage import morphology as mm
from skimage.util import img_as_ubyte
from skimage.filters import threshold_otsu, rank

def l7smooth( u ):
    # Landsat7 Correction if SLC-off
    return np.where( np.isnan(u), mm.closing(u, mm.disk(12)), u )


def Thresholding( rgb, band=None ):
    '''Thresholding(rgb) - Apply Otsu's Thresholding Method'''
    # Assign band
    if band is None: idx = 0 # Defaut band is R
    elif isinstance(band, str):
        if band.lower()[0] == 'r': idx = 0
        elif band.lower()[0] == 'g': idx = 1
        elif band.lower()[0] == 'b': idx = 2
    # Apply Threshold to selected Band
    img = rgb[:,:,idx] # Band Index
    thresh = threshold_otsu( img ) # Compute Otsu's Threshold
    bw = img < thresh # Apply Threshold
    return bw


def SegmentationIndex( *args, **kwargs ):
    '''Apply Index'''

    R = kwargs['R'].astype( float )
    G = kwargs['G'].astype( float )
    B = kwargs['B'].astype( float )
    NIR = kwargs.pop( 'NIR', np.full(R.shape,np.nan) ).astype( float )
    MIR = kwargs.pop( 'MIR', np.full(R.shape,np.nan) ).astype( float )
    SWIR = kwargs.pop( 'SWIR', np.full(R.shape,np.nan) ).astype( float )
    index = kwargs.pop( 'index', None )
    rad = kwargs.pop( 'radius', 20 )
    method = kwargs.pop( 'method', 'local' )
    L7correction = kwargs.pop( 'L7correction', True )

    if index == 'NDVI':
        IDX =  (R - NIR) / (R + NIR)
    elif index == 'MNDWI':
        IDX =  (G - MIR) / (G + MIR)
    elif index == 'BAR':
        IDX = SWIR
    elif index == 'MIX':
        IDX = (R - NIR) / (R + NIR)
        IDXX = (G - MIR) / (G + MIR)
        IDXXX = SWIR
    elif index == 'AWEI':
        raise NotImplementedError
    else:
        err = 'Index %s not recognized' % IDX
        raise ValueError(err)

    if L7correction:
        IDX = l7smooth( IDX )
        if index == 'MIX':
            IDXX = l7smooth( IDXX )
            IDXXX = l7smooth( IDXXX )

    # Apply Local Otsu's Method
    selem = mm.disk( rad )
    globthresh = threshold_otsu( IDX[np.isfinite(IDX)] )

    if index=='MIX':
        globthreshX = threshold_otsu( IDXX[np.isfinite(IDXX)] )
        globthreshXX = threshold_otsu( IDXXX[np.isfinite(IDXXX)] )

    if method == 'local':
        print("applying local Otsu method - this may require some time... ", )
        thresh = rank.otsu( img_as_ubyte(IDX), selem ).astype(float)
        if index=='MIX':
            threshX = rank.otsu( img_as_ubyte(IDXX), selem ).astype(float)
            threshXX = rank.otsu( img_as_ubyte(IDXXX), selem ).astype(float)
        print('done')
    else:
        thresh = globthresh
        if index=='MIX':
            threshX = globthreshX
            threshXX = globthreshXX

    if index == 'MIX':
        MASK = np.logical_or( ( IDX>=thresh ) * ( mm.binary_dilation(IDXX>=threshX,mm.disk(0.3*rad)) ), IDXXX>threshXX)
    else: MASK = IDX >= thresh

    return IDX, MASK.astype( int ), globthresh

