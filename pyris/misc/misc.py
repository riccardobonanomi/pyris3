# ==================================================
# Module: misc
# Package: PyRIS
# Author: Federico Monegaglia
# Date: April 2016
# Description: Collection of miscellaneous functions
# ==================================================

from __future__ import division
import os
import time
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from skimage.io import imread
from skimage.morphology import closing, disk
from osgeo import gdal
import sys
import warnings

# Matplotlib Setting
mpl.rc( 'font', size=20 )
mpl.rcParams['text.latex.preamble'] = r"\usepackage{siunitx} \siteup{detect-all} \usepackage{helvet} \usepackage{sansmath} \usepackage{amsmath} \sansmath"
mpl.rcParams['axes.formatter.limits'] = [-3,3]
plt.close( 'all' )


# Functions
# =========
def ediff1d0( x ):
    '''
    Compute the element-wise different of an array starting from 0
    '''
    if len( x ) == 0:
        return np.ediff1d( x )
    return np.ediff1d( x, to_begin=0 )


def NaNs( N ):
    '''
    Create a 1d array of NaNs of size N
    '''
    return np.full( N, np.nan )


def crossprod2( v1, v2 ):
    '''
    3rd component of 2d vectors cross product
    '''
    return v1[0]*v2[1] - v1[1]*v2[0]


def Intersection( P, Q, R, S, return_point=True ):
    '''Check if two segments PQ and RS intersect and where'''
    QP = Q - P
    CRS = crossprod2( R, S ) 
    t = ( QP[0]*S[1] - S[0]*QP[1] ) / CRS
    u = ( QP[0]*R[1] - R[0]*QP[1] ) / CRS
    if ( abs(CRS) > 0 and 0 <= abs(t) <= 1 and 0 <= abs(u) <= 1 ):
        # Segments Intersect!
        if return_point: return True, P + t*R
        else: return True
    if return_point: return False, NaNs(2)
    else: return False


def LoadLandsatData( dirname ):
    '''
    Load Relevant Bands for the Current Landsat Data
    '''
    if any( [os.path.split(dirname)[-1].startswith( s ) for s in ['LC8', 'LC08', 'LC9', 'LC09']] ): bidx = range( 2, 8 )
    else: bidx = [1,2,3,4,5,7]
    base = os.path.join( dirname, os.path.basename(dirname) )
    ext = '.TIF'
    bnames = [ ('_B'.join(( base, '%d' % i )))+ext for i in bidx ]
    [ B, G, R, NIR, MIR, SWIR ] = [ imread( band ) for band in bnames ]
    bands = [ R, G, B, NIR, MIR, SWIR ]

    geo = gdal.Open( bnames[0] )
    GeoTransf = {    
        'PixelSize' : abs( geo.GetGeoTransform()[1] ),
        'X' : geo.GetGeoTransform()[0],
        'Y' : geo.GetGeoTransform()[3],
        'Lx' : bands[0].shape[1],
        'Ly' : bands[0].shape[0]
        }
    return bands, GeoTransf

def LoadRawMask( maskfile ):
    '''
    Load raw .tif mask for the current raw mask data
    '''
    mask = imread( maskfile )

    geo = gdal.Open( maskfile )
    GeoTransf = {    
        'PixelSize' : abs( geo.GetGeoTransform()[1] ),
        'X' : geo.GetGeoTransform()[0],
        'Y' : geo.GetGeoTransform()[3],
        'Lx' : mask.shape[1],
        'Ly' : mask.shape[0]
        }
    return mask, GeoTransf


# Classes
# =======
class Line2D( object ): # This must be put down better!
    '''
    Line2D - a class to store and manage channel centerline coordinates and width
    '''
    def __init__( self, x=np.array([]), y=np.array([]), B=np.array([]) ):
        '''
        Read centerline and width sequences or initialize empty ones
        '''
        dx = ediff1d0( x )
        dy = ediff1d0( y )
        ds = np.sqrt( dx**2 + dy**2 )
        s = np.cumsum( ds )
        try: L = s[-1]
        except IndexError: L = 0
        self.x, self.y, self.B, self.s, self.L = x, y, B, s, L
        return None

    def join( self, line2d ):
        '''
        Join two subsequent centerlines
        '''        
        if len(self.x) == 0:
            ds = 0
        else:
            ds = np.sqrt( (self.x[-1]-line2d.x[-1])**2 + (self.y[-1]-line2d.y[-1])**2 )
        self.x = np.concatenate( (self.x, line2d.x) )
        self.y = np.concatenate( (self.y, line2d.y) )
        self.B = np.concatenate( (self.B, line2d.B) )
        self.s = np.concatenate( (self.s, line2d.s+ds) )
        self.L = self.s[-1]
        return None

    
class GeoReference( object ):
    '''Provide Georeferenced Coordinates for an Image Object'''

    def __init__( self, GeoTransf ):
        '''
        Store a GeoTransf and define the extent of the EO image for viewing
        '''
        self.GeoTransf = GeoTransf
        self.extent = [ GeoTransf['X'], # xmin
                        GeoTransf['X'] + GeoTransf['PixelSize']*GeoTransf['Lx'], # xmax
                        GeoTransf['Y'] - GeoTransf['PixelSize']*GeoTransf['Ly'], # ymin
                        GeoTransf['Y'] # ymax
                    ]
        return None

    def RefCurve( self, X, Y, inverse=False ):
        '''Compute pixel-to-georeferenced coordinate conversion (and inverse)'''
        X, Y = np.asarray(X), np.asarray(Y)
        if inverse:
            Cx = ( X - self.GeoTransf['X'] ) / self.GeoTransf['PixelSize']
            Cy = -( Y - self.GeoTransf['Y'] ) / self.GeoTransf['PixelSize']
            return Cx, Cy
        else:
            self.Cx, self.Cy = X, Y
            self.CX = self.GeoTransf['X'] + self.Cx*self.GeoTransf['PixelSize']
            self.CY = self.GeoTransf['Y'] - self.Cy*self.GeoTransf['PixelSize']
            return self.CX, self.CY


class interactive_mask( object ):

    '''
    A callable class where the user can view an False Color composite
    and select areas to include/exclude from the analysis
    '''
    
    def __init__( self, path ):
        '''
        Path to Landsat directory
        '''
        self.path = os.path.normpath( path )
        if os.path.isdir(self.path):
            self.name = self.path.split( os.sep )[-1]
            self.landsat = True
        else:
            self.name = self.path.split( os.sep )[-2]
            self.landsat = False

    def build_real_color( self ):
        '''
        Create an RGB  False Color composite
        '''
        if any( [self.name.startswith( s ) for s in ['LC8', 'LC08', 'LC9', 'LC09']] ):
            w = 'Landsat 8 and 9 may return distorted images as real color.'
            warnings.warn( w, UserWarning )
            b1, b2, b3 = 'B6', 'B5', 'B4'
        else:
            b1, b2, b3 = 'B5', 'B4', 'B3'
        B1, B2, B3 = map( imread, [ os.path.join( self.path, '_'.join((self.name,bname))+'.TIF' ) for bname in [ b1, b2, b3 ] ] )
        return np.dstack( ( B1, B2, B3 ) )

    def _set_mask( self ):
        '''
        Interactively select areas through drag-and-drop
        '''
        white_masks = []
        plt.ioff()
        fig = plt.figure()
        plt.title( 'Press-drag a rectangle for your mask. Close when you are finish.' )
        if self.landsat:
            real_color = self.build_real_color()
            plt.imshow( real_color, cmap='binary_r' )
        else:
            raw_bw = imread (self.path)
            cmap = ListedColormap(['black', 'white'])
            plt.imshow(raw_bw, cmap = cmap)
        plt.axis('equal')
        x_press = None
        y_press = None
        def onpress(event):
            global x_press, y_press
            x_press = int(event.xdata) if (event.xdata != None) else None
            y_press = int(event.ydata) if (event.ydata != None) else None
        def onrelease(event):
            global x_press, y_press
            x_release = int(event.xdata) if (event.xdata != None) else None
            y_release = int(event.ydata) if (event.ydata != None) else None
            if (x_press != None and y_press != None
                and x_release != None and y_release != None):
                (xs, xe) = (x_press, x_release+1) if (x_press <= x_release) \
                  else (x_release, x_press+1)
                (ys, ye) = (y_press, y_release+1) if (y_press <= y_release) \
                  else (y_release, y_press+1)
                print( "The mask you selected is [{0}:{1},{2}:{3}]".format(
                    xs, xe, ys, ye) )
                white_masks.append( [ ys, ye, xs, xe ] )
                plt.fill( [xs,xe,xe,xs,xs], [ys,ys,ye,ye,ys], 'r', alpha=0.25 )
                event.canvas.draw()
            x_press = None
            y_press = None
        cid_press   = fig.canvas.mpl_connect('button_press_event'  , onpress  )
        cid_release = fig.canvas.mpl_connect('button_release_event', onrelease)
        plt.show()
        return white_masks

    def get_georef( self ):
        if self.landsat:
            bands, GeoTransf = LoadLandsatData( self.path )
        else:
            bands, GeoTransf = LoadRawMask( self.path )
        return GeoReference( GeoTransf )

    def _georeference_masks( self, masks, inverse=False ):
        GRF = self.get_georef()
        gmask = []
        for mask in masks:
            [ys, ye, xs, xe] = mask
            X, Y = GRF.RefCurve( np.array([xs,xe]), np.array([ys,ye]), inverse=inverse )
            gmask.append( [ Y[0], Y[1], X[0], X[1] ] )
        return gmask
    
    def georeference( self, masks ):
        gmasks = self._georeference_masks( masks )
        return gmasks

    def dereference( self, gmasks ):
        masks = self._georeference_masks( gmasks, inverse=True )
        return masks

    def __call__( self, *args, **kwargs ):
        inverse = kwargs.pop( 'inverse', False )
        if inverse:
            gmasks = list(sys.argv[1] )
            return self.dereference( gmasks )
        masks = self._set_mask()
        return self.georeference( masks )

    
class MaskClean( object ):

    '''
    A callable class where the user can select areas to include/exclude from the analysis
    directly from a mask
    '''
    
    def __init__( self, bw, bg=None ):
        self.bw = bw.astype( int )
        self.bg = np.zeros(bw.shape) if bg is None else bg
        self.press = False
        
    def __call__( self ):
        #real_color = self.build_real_color()
        white_masks = []
        plt.ioff()
        fig = plt.figure()
        plt.title( 'Press-drag a rectangle for your mask. Close when you are finish.' )
        plt.imshow( self.bg, cmap='binary_r' )
        plt.imshow( self.bw, cmap='jet', alpha=0.5 )
        plt.axis('equal')
        x_press = None
        y_press = None
        def onpress(event):
            global x_press, y_press
            x_press = int(event.xdata) if (event.xdata != None) else None
            y_press = int(event.ydata) if (event.ydata != None) else None
        def onrelease(event):
            global x_press, y_press
            x_release = int(event.xdata) if (event.xdata != None) else None
            y_release = int(event.ydata) if (event.ydata != None) else None
            if (x_press != None and y_press != None
                and x_release != None and y_release != None):
                (xs, xe) = (x_press, x_release+1) if (x_press <= x_release) \
                  else (x_release, x_press+1)
                (ys, ye) = (y_press, y_release+1) if (y_press <= y_release) \
                  else (y_release, y_press+1)
                print( "Slice [{0}:{1},{2}:{3}] will be set to {4}".format(
                    xs, xe, ys, ye, 0) )
                self.bw[ ys:ye,xs:xe ] = 0
                self.press = True
                plt.fill( [xs,xe,xe,xs,xs], [ys,ys,ye,ye,ys], 'r', alpha=0.25 )
                event.canvas.draw()
            x_press = None
            y_press = None
        cid_press   = fig.canvas.mpl_connect('button_press_event'  , onpress  )
        cid_release = fig.canvas.mpl_connect('button_release_event', onrelease)
        plt.show()
        return self.bw, self.press


class ToC( object ):
    '''
    Compute and print the computation time
    '''
    def __init__(self, message = 'Computation'):
        self.start = time.time()
        self.message = message
        return None

    def print_elapsed_time(self):
        self.end = time.time()
        # convert time in seconds to hours, minutes and seconds
        hours, rem = divmod(self.end - self.start, 3600)
        minutes, seconds = divmod(rem, 60)
        if hours == 0:
            if minutes == 0:
                end_message = '%.2f seconds' % seconds
            else:
                end_message = '%d minutes and %04.2f seconds' % (minutes, seconds)
        else:
            end_message = '%d hours, %02d minutes and %04.2f seconds' % (hours, minutes, seconds)
        print(self.message + ' ended in ' + end_message)
        return None