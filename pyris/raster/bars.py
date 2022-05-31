# ===========================================================
# Module: vector
# File: axis.py
# Package: PyRIS
# Author: Federico Monegaglia
# Date: April 2016
# Description: skeleton centerline extractor and vectorizator
# ===========================================================

from __future__ import division
import os, sys
import numpy as np
from scipy import ndimage, interpolate as scipy_interp
from skimage import morphology as mm
from skimage.feature import peak_local_max
from skimage import measure as sm
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from ..raster.segmentation import SegmentationIndex
from ..misc.misc import GeoReference, NaNs


class Unwrapper( object ):

    '''
    Unwrapper( object )
    ===================

    Unwrap the river channel by performing coordinate transform
    (x,y) --> (s,n)

    Arguments
    -----
    data        Output for a single river planform
    mig         Output of the migration result for the current river planform
    GeoTransf   Geospatial Coordinate Transform

    Methods
    -------
    unwrap( shape, Npts=100 )   Regrids the river planform over intrinsic and cartesian coordinates
    interpolate( band )         Interpolates the band values over the regridded reference system    
    '''

    def __init__( self, data, mig, GeoTransf ):
        '''Read the Georeferenced River Planform'''
        self.data = data
        self.GeoTransf = GeoTransf
        self.x = data[0]
        self.y = data[1]
        self.s = data[2]
        self.theta = data[3]
        self.Cs = data[4]
        self.b = data[5]
        self.Bend = mig[4].astype( int )
        self.NextBend = mig[5].astype( int )
        self.Ipoints = ( np.where( data[6].astype( int )==2 )[0] ).astype( int )
        self.BendIndexes = np.unique( self.Bend[ self.Bend>=0 ] )
        return None

    def unwrap( self, shape, Npts=100 ):

        '''
        Returns the Cartesian Coordinates (XC,YC) together with the associated Intrinsic Coordinates (s, n)
        as MeshGrids
        '''
        x, y, s, b, GeoTransf = self.x, self.y, self.s, self.b, self.GeoTransf
        # Pixel units river planform    
        X = ( x -  ( GeoTransf['X'] ) ) / GeoTransf['PixelSize']
        Y = ( y -  ( (1-shape[0])*GeoTransf['PixelSize'] + GeoTransf['Y'] ) ) / GeoTransf['PixelSize']
        S = s / GeoTransf['PixelSize']
        B = b / GeoTransf['PixelSize']

        # Transverse Axis
        N = np.linspace( -1.0, +1.0, Npts, endpoint=True )
        self.N = N
        # Knots
        self.Xc, self.Yc = np.zeros((S.size, N.size)), np.zeros((S.size, N.size))
        angle = np.arctan2( np.gradient( Y ),  np.gradient( X ) )

        # Create Cartesian Coorinates Array for Intrinsic Coordinate Grid
        for i in xrange( S.size ):
            n = N * B[i] # Pixel Units Transverse Coordinate
            self.Xc[i,:] = X[i] + n[:]*np.cos( angle[i]-np.pi/2 )
            self.Yc[i,:] = Y[i] + n[:]*np.sin( angle[i]-np.pi/2 )

        self.I, self.J = np.arange( S.size ), np.arange( N.size )

        self.XC = self.Xc*GeoTransf['PixelSize'] + GeoTransf['X']
        self.YC = self.Yc*GeoTransf['PixelSize'] + GeoTransf['Y'] + (1-shape[0])*GeoTransf['PixelSize']
        self.Sc, self.Nc = np.meshgrid( s/b.mean(), N )

        return [self.XC, self.YC],  [self.Sc, self.Nc]


    def interpolate( self, band ):
        '''Interpolate a band or a bands combination over the grid'''
        return ndimage.interpolation.map_coordinates( band[::-1,:], [self.Yc, self.Xc] )



class BarFinder( object ):

    '''
    BarFinder( object )
    ===================

    A callable class to detect bare sediment bars along channels

    Arguments
    ---------
    unwrapper         An 'Unwrapper' instance
    bands             A sequence of EO bands as extracted by LoadLandsatData
    close             Apply binary_closing to the bar mask (default True)
    remove_small      Remove very small detected bars (default True)

    '''


    def __init__( self, unwrapper ):
        '''
        BarFinder(unwrapper) - Read Coordinate Unwrapper
        '''
        self.unwrapper = unwrapper
        self.BendIndexes = self.unwrapper.BendIndexes
        return None


    def FindBars( self, bands, close=True, remove_small=True ):
        '''
        BarFinder.FindBars(bands) - Uses Otsu's global thresholding
                                    with both MNDWI and NDVI indexes
                                    in order to find channel bars.
                                    A convex hull is applied to the
                                    unwrapper bar objects.

        Args
        ====
        bands: dictionary - must contain R, G, B, NIR, MIR bands

        Kwargs
        ======
        close: boolean - perform binary closing on bar mask

        Returns
        =======
        FindBars: Labeled Channel Bars
        '''
        
        Wbands = { } # Water bands
        Vbands = { } # Vegetation bands

        for nband in [ 'R', 'G', 'B', 'NIR', 'MIR', 'SWIR' ]:
            Wbands[nband] = ndimage.interpolation.map_coordinates(
                bands[nband][::-1,:], [self.unwrapper.Yc, self.unwrapper.Xc] )

        Idx, Bars, otsu_glob = SegmentationIndex(
            R=Wbands['R'], G=Wbands['G'], B=Wbands['B'],
            NIR=Wbands['NIR'], MIR=Wbands['MIR'], SWIR=Wbands['SWIR'],
            index='BAR', method='global' ) ## This must be made locally, otherwise we dont see bars eventually
        
        # Apply a Convex Hull to Channel Bars
        #Bars = mm.convex_hull_object( Bars )

        if close:
            # 1/8 of the total average channel width is used ( 1/4*mean(b) )
            rad = 1 #max( 0, self.unwrapper.b.mean()/self.unwrapper.GeoTransf['PixelSize'] )
            Bars = mm.binary_closing( Bars, mm.disk(rad) )

        if remove_small:
            Amin = 0.1*self.unwrapper.N.size/2 * ( self.unwrapper.b.mean() /  (self.unwrapper.s[1]-self.unwrapper.s[0]) )
            mm.remove_small_objects( Bars, 1, in_place=True ) # Remove small Bars
            mm.remove_small_holes(   Bars, Amin, in_place=True ) # Remove Internal Spots

        # Apply a Convex Hull to Channel Bars
        Bars = mm.convex_hull_object( Bars )

        # Identify the Largest Bar on the Bend
        labeled_array, num_features = ndimage.measurements.label( Bars )
        self.Bars = labeled_array
        self.BarIdx = np.arange( num_features, dtype=int )+1
        return self.Bars

    def BarCentroid( self ):
        '''Compute Centroid of Each Bar'''

        num_features = self.BarIdx.max() # Number of Bar Features
        IC, JC = np.zeros(num_features,dtype=int), np.zeros(num_features,dtype=int) # Indexes of Bar Centroids
        # Compute Bar Centroids
        for n in xrange( 1,num_features ): [ IC[n-1], JC[n-1] ] = ndimage.measurements.center_of_mass( self.Bars==(n) )
        # Apply a Correction to the Transverse Position of the Centroid: we take the position of the maximum longitudinal extension (main s axis)
        for n in xrange( 1,num_features ):
            s_lengths = (self.Bars==n).sum(axis=0)
            pos = np.where( s_lengths == s_lengths.max() )[0]
            JC[n-1] = int( ( pos.min() + pos.max() ) / 2 )
        self.Centroid = np.vstack((IC, JC))
        return IC, JC


    def BarArea( self ):
        ''''Measure the number of Pixels of each Bar'''
        self.Area = np.zeros( self.BarIdx.size, dtype=int )
        for n in self.BarIdx: self.Area[n-1] = np.sum( self.Bars==n )
        return self.Area

    def BarType( self ):
        '''
        Define if a Bar is either Central, Lateral or Mixed

        ---------- - +1     ^ N
        2222222222 _        |
        1111111111   +0.5   |
        0000000000 -  0    -|--------> S
        1111111111 _ -0.5   |
        2222222222          |
        ---------- - -1     |
        '''
        JC = self.Centroid[1]
        NC = self.unwrapper.N[ JC ]
        TYPE = np.ones( NC.size, dtype=int ) # MiXeD
        TYPE[ np.abs(NC)<0.25 ] = 0 # MidChannel
        TYPE[ np.abs(NC)>0.5 ] = 2 # Lateral        
        self.TYPE = TYPE
        self.Central = self.BarIdx[ TYPE==0 ]
        self.Mixed = self.BarIdx[ TYPE==1 ]
        self.Lateral = self.BarIdx[ TYPE==2 ]
        self.TYPES = [ self.Central, self.Mixed, self.Lateral ]
        return TYPE

    def BarBend( self ):
        '''Return the Bar-Bend Index (to which bend the bar belongs)'''
        self.BBIdx = np.full( self.Centroid[0].size, -1, dtype=int )
        for j, c in enumerate( self.Centroid[0] ):
            for i in self.BendIndexes:
                mask = self.unwrapper.Bend==i
                s = self.unwrapper.s[mask]
                if s[0] < self.unwrapper.s[c] <= s[-1]:
                    self.BBIdx[j] = i
                    break
        return self.BBIdx

    def MainBarTypeBend( self, TYPE=None ):
        '''
        Identify the main bar for a specific bar type for each bend
        None : Any Kind of Bar
        0    : Central
        1    : Mixed
        2    : Lateral
        '''
        if TYPE is None: BBIdx = self.BBIdx
        else: BBIdx = self.TYPES[ TYPE ]

        BIdx = np.unique( BBIdx[ BBIdx>=0 ] )
        BarBendTypeIdx = -np.ones( BIdx.size, dtype=int )
        for i, n in enumerate( BIdx ):
            Areas =  np.asarray( [ (self.Bars==(k+1)).sum() for ki, k in enumerate(self.BarIdx[ BBIdx==n ]) ] )
            Indexes = self.BarIdx[ BBIdx==n ]
            BarBendTypeIdx[i] = Indexes[ Areas.argmax() ]
        return BarBendTypeIdx


    def MainBarBend( self ):
        '''Identify the main bar for each bend'''
        self.BarBendIdx = self.MainBarTypeBend()
        return self.BarBendIdx


    def BarContour( self ):
        '''Compute the contour line of each emerging bar'''
        self.Contours = []
        for n in self.BarIdx:
            contour = sm.find_contours( self.Bars==n, 0.5 )
            I, J = contour[0][:,0], contour[0][:,1]

            if np.abs(I[-1]-I[0])>1:
                _I = np.arange( I[-1], I[0], np.sign(I[0]-I[-1]) )
                _J = np.ones( _I.size ) * J[-1]
                I, J = np.append( I, _I ), np.append( J, _J )

            self.Contours.append( [I, J] )
        return self.Contours            


    def BarProps( self ):
        '''Compute Properties of Channel Bars'''
        self.BarCentroid()
        self.BarArea()
        self.BarType()
        self.BarBend()    
        self.BarContour()
        self.MainBarBend()
        return None


    def Show( self, bands ):

        plt.figure()
        gs = gridspec.GridSpec( 12, 2 )
        ax1 = plt.subplot( gs[:1, :] ) # Colorbar
        ax2 = plt.subplot( gs[2:, 0] ) # Topography
        ax3 = plt.subplot( gs[2:, 1], sharex=ax2, sharey=ax2 ) # RGB        
        GR = GeoReference( self.unwrapper.GeoTransf )
        XC, YC = self.unwrapper.XC, self.unwrapper.YC
        Zc = self.unwrapper.interpolate( bands['B'] )
        ax2.imshow( bands['G'], cmap='binary', extent=GR.extent )
        pcm = ax2.pcolormesh( XC, YC, Zc, cmap='Spectral_r', alpha=0.75 )
        ax2.contour( XC, YC, Zc )
        plt.plot( XC[:,0], YC[:,0], 'k', lw=2 )
        plt.plot( XC[:,-1], YC[:,-1], 'k', lw=2 )    
        Zi = np.dstack( (bands['MIR'], bands['NIR'], bands['R']) )
        ax3.imshow( Zi, extent=GR.extent )
        cb = plt.colorbar( pcm, cax=ax1, orientation='horizontal' )
        plt.axis('equal')

        plt.figure()
        gs = gridspec.GridSpec( 12, 2 )
        ax1 = plt.subplot( gs[:1, :] ) # Colorbar
        ax2 = plt.subplot( gs[2:, 0] ) # Topography
        ax3 = plt.subplot( gs[2:, 1], sharex=ax2, sharey=ax2 ) # RGB        
        GR = GeoReference( self.unwrapper.GeoTransf )
        XC, YC = self.unwrapper.XC, self.unwrapper.YC
        Zc = self.unwrapper.interpolate( bands['B'] )
        ax2.imshow( bands['G'], cmap='binary', extent=GR.extent )
        pcm = ax2.pcolormesh( XC, YC, Zc, cmap='Spectral_r', alpha=0.75 )
        ax2.contour( XC, YC, Zc )
        plt.plot( XC[:,0], YC[:,0], 'k', lw=2 )
        plt.plot( XC[:,-1], YC[:,-1], 'k', lw=2 )    
        # Draw Cross Sections
        for i in xrange(1, self.unwrapper.s.size-1, 10):
            ax2.plot( XC[i,:], YC[i,:], 'k' )
            ax2.text( XC[i,-1], YC[i,-1], 's=%s' % int(self.unwrapper.s[i]/self.unwrapper.b.mean()) )
            ax3.plot( XC[i,:], YC[i,:], 'k' )
            ax3.text( XC[i,-1], YC[i,-1], 's=%s' % int(self.unwrapper.s[i]/self.unwrapper.b.mean()) )
        Zi = np.dstack( (bands['MIR'], bands['NIR'], bands['R']) )
        ax3.imshow( Zi, extent=GR.extent )
        cb = plt.colorbar( pcm, cax=ax1, orientation='horizontal' )
        plt.axis('equal')
    
        plt.figure()
        gs = gridspec.GridSpec( 60, 60 )
        ax1 = plt.subplot( gs[7:11, :] ) # Colorbar
        ax2 = plt.subplot( gs[15:25,:] ) # Surface
        ax3 = plt.subplot( gs[30:40, :], sharex=ax2 ) # Width
        ax4 = plt.subplot( gs[48:58, :], sharex=ax2 ) # Curvature    
        # Surface
        pcm = ax2.pcolormesh( self.unwrapper.s/self.unwrapper.b.mean(), self.unwrapper.N, Zc.T, cmap='Spectral_r' )
        ax2.contour( self.unwrapper.Sc, self.unwrapper.Nc, Zc.T )
        ax2.set_ylabel( r'$n^*/B^*$' )    
        # Colorbar
        cb = plt.colorbar( pcm, cax=ax1, orientation='horizontal' )        
        # Width
        ax3.plot( self.unwrapper.s/self.unwrapper.b.mean(), self.unwrapper.b/self.unwrapper.b.mean() )
        ax3.set_ylabel( r'$B^*/B_0^*$' )    
        # Curvature
        ax4.plot( self.unwrapper.s/self.unwrapper.b.mean(), self.unwrapper.Cs )
        ax4.set_ylabel( r'$\mathcal{C^*}$' )
        ax4.set_xlabel( r'$s^*/B_0^*$' )
        plt.axis('tight')

        plt.show()

    
    def __call__( self, bands, close=True, remove_small=True ):
        self.FindBars( bands, close, remove_small )
        self.BarProps()
        return None



class TemporalBars( object ):

    '''
    TemporalBars( object )
    ======================

    A class to reconstruct the temporal dynamics of bare sediment bars along rivers
    
    '''
    def __init__( self ):
        '''Perform a Temporal Analysis of Channel Bars'''
        self.T = []
        self.Bars = []
        return None

    def GetFinder( self, T, Finder ):
        '''Read a BarFinder Instance with its time and store it'''
        self.T.append( T )
        self.Bars.append( Finder )
        return None

    def IterData( self ):
        '''Bar Iterator'''
        for ti, ( T, Bar ) in enumerate( zip( self.T, self.Bars ) ):
            yield ti, ( T, Bar )

    def BendIndexes( self, ti ):
        '''Return the sequence of Bend Indexes for a given time index'''
        return self.Bars[ti].BendIndexes

    def IterBends( self, ti=0 ):
        '''Iterate over individual meander bends at a given time index ti'''
        for idx in self.BendIndexes( ti ):
            yield idx

    def CentroidsEvol( self, bend_idx, normalize=True ):
        '''Follow the evolution of the centroid of main bar of an individual meander bend'''

        bend = bend_idx

        centroids_IJ = []
        centroids_XY = []
        centroids_SN = []
        
        for ti, ( T, Bars ) in self.IterData():

            if bend == -1: break # CutOff is Reached
            
            # Find the Current Bend
            bend_indexes = Bars.unwrapper.Bend
            bend_indexes_next = Bars.unwrapper.NextBend
            mask = bend_indexes==bend

            # Find the Dominant Bar of the Current Bend
            try: ibar = Bars.BarBendIdx[ bend ]
            except IndexError: continue ## ??? what's wrong here ???
            
            # Get Bar Properties
            IC, JC = Bars.Centroid[ :,ibar ]
            X = Bars.unwrapper.XC[ IC, JC ]
            Y = Bars.unwrapper.YC[ IC, JC ]
            S = Bars.unwrapper.s[ IC ]
            N = Bars.unwrapper.N[ -JC ]

            if normalize:
                Sbend = Bars.unwrapper.s[ mask ]
                SS = S # Bend-Normalized S
                SS -= Sbend[ int( Sbend.size / 2 ) ] # Relative to "Bend Apex" (FIXME: use proper apex)
                SS /= (0.5*( Sbend[-1] - Sbend[0] )) # Normalize to Bend Half-Length

            centroids_IJ.append( [ IC, JC ] )
            centroids_XY.append( [ X, Y ] )
            centroids_SN.append( [ SS, N ] )

            # Get the Bend Index for the Next Time Step
            bend = bend_indexes_next[ mask ][0]

        self.centroids_IJ, self.centroids_XY, self.centroids_SN = centroids_IJ, centroids_XY, centroids_SN

        return centroids_IJ, centroids_SN, centroids_XY

    def MainBarEvol( self, bend_idx, normalize=True ):
        '''Follow the evolution of the main bar of an individual meander bend'''

        bend = bend_idx

        contours_XY = []
        contours_SN = []
        contours_IJ = []
        
        for ti, ( T, Bars ) in self.IterData():

            if bend == -1: break # CutOff is Reached
            
            # Find the Current Bend
            bend_indexes = Bars.unwrapper.Bend
            bend_indexes_next = Bars.unwrapper.NextBend
            mask = bend_indexes==bend

            # Find the Dominant Bar of the Current Bend
            try: ibar = Bars.BarBendIdx[ bend ]
            except IndexError: continue ## ??? what's wrong here ???
            
            # Get Bar Properties
            contour = Bars.Contours[ ibar ]
            X = Bars.unwrapper.XC[ contour[0].astype(int),contour[1].astype(int) ]
            Y = Bars.unwrapper.YC[ contour[0].astype(int),contour[1].astype(int) ]
            S = Bars.unwrapper.s[ contour[0].astype(int) ]
            N = -Bars.unwrapper.N[ contour[1].astype(int) ]

            if normalize:
                Sbend = Bars.unwrapper.s[mask]
                S -= Sbend[ int( Sbend.size / 2 ) ] # Relative to "Bend Apex" (FIXME: use proper apex)
                S /= ( 0.5*(Sbend[-1] - Sbend[0]) ) # Normalize to Bend Half-Length
            
            contours_IJ.append( [ contour[0], contour[1] ] )
            contours_XY.append( [ X, Y ] )
            contours_SN.append( [ S, N ] )

            # Get the Bend Index for the Next Time Step
            bend = bend_indexes_next[ mask ][0]

        self.contours_IJ, self.contours_XY, self.contours_SN = contours_IJ, contours_XY, contours_SN

        return contours_IJ, contours_SN, contours_XY


    def Show( self, landsat_dirs, geodir, bend=None ):

        for BEND in self.IterBends():
            if bend is not None:
                if not BEND==bend: continue

            if hasattr( self, 'contours_IJ' ): [ contours_IJ, contours_SN, contours_XY ] = [ self.contours_IJ, self.contours_SN, self.contours_XY ]
            else: [ contours_IJ, contours_SN, contours_XY ] = self.MainBarEvol( BEND )

            if hasattr( self, 'centroids_IJ' ): [ centroids_IJ, centroids_SN, centroids_XY ] = [ self.centroids_IJ, self.centroids_SN, self.centroids_XY ]
            else: [ centroids_IJ, centroids_SN, centroids_XY ] = self.CentroidsEvol( BEND )
    
            colors = [ cm.Spectral_r(k) for k in np.linspace(0, 1, len(contours_SN)) ]
            lws = np.linspace( 0.5, 2.5, len(contours_SN) )
        
            smin, smax = 0, 0
            xmin, xmax = contours_XY[0][0].min(), contours_XY[0][0].max()
            ymin, ymax = contours_XY[0][1].min(), contours_XY[0][1].max()
            for i, (Csn, Cxy) in enumerate( zip( contours_SN, contours_XY ) ):
                s = np.append( Csn[0], Csn[0][-1] )
                n = np.append( Csn[1], Csn[1][-1] )
                x = np.append( Cxy[0], Cxy[0][-1] )
                y = np.append( Cxy[1], Cxy[1][-1] )
                smin, smax = min( smin, s.min() ), max( smax, s.max() )
                xmin, xmax = min( xmin, x.min() ), max( xmax, x.max() )
                ymin, ymax = min( ymin, y.min() ), max( ymax, y.max() )

            for i, (Csn, Cxy) in enumerate( zip( contours_SN, contours_XY ) ):
                plt.figure()
                gs = gridspec.GridSpec( 60, 60 )
                ax1 = plt.subplot( gs[5:55,:23] ) # SN
                ax2 = plt.subplot( gs[5:55,28:51] ) # XY

                geofile = sorted(os.listdir(geodir))[i]
                name = ''.join((os.path.splitext( os.path.basename( geofile ) )[0].split('_')))
                found = False
                for landsat_dir in landsat_dirs:
                    lname = os.path.splitext(os.path.split(landsat_dir)[-1])[0][9:16]
                    if name == lname:
                        found = True
                        break
                if found == True:
                    from ..misc import LoadLandsatData
                    [ R, G, B, NIR, MIR ], GeoTransf = LoadLandsatData( landsat_dir, L7correction=False )
                    bands = { 'R' : R, 'G' : G, 'B' : B, 'NIR' : NIR, 'MIR' : MIR }                    
                    ax2.imshow( np.dstack( (bands['MIR'], bands['NIR'], bands['R']) ), extent=GeoReference(GeoTransf).extent )

                if i>0 and (centroids_SN[i][0]-centroids_SN[i-1][0]) > 4: continue
                for j, (csn, cxy) in enumerate( zip( contours_SN[:i+1], contours_XY[:i+1] ) ):
                    s = np.append( csn[0], csn[0][-1] )
                    n = np.append( csn[1], csn[1][-1] )
                    x = np.append( cxy[0], cxy[0][-1] )
                    y = np.append( cxy[1], cxy[1][-1] )
                    ax1.plot( s, n, label=r'%s' % int(self.T[j]), lw=4, c=colors[j], alpha=0.75 )
                    ax1.fill( s, n, color=colors[j], alpha=0.5 )
                    ax2.plot( x, y, label=r'%s' % int(self.T[j]), lw=4, c=colors[j], alpha=0.75 )
                    ax2.fill( x, y, color=colors[j], alpha=0.5 )
                ax1.set_xlabel( r'$s/B_0\, [-]$' )
                ax1.set_ylabel( r'$n/B_0\, [-]$' )
                ax2.set_xlabel( r'$x_{\mathrm{UTM}} [\mathrm{m}]$' )
                ax2.set_ylabel( r'$y_{\mathrm{UTM}} [\mathrm{m}]$' )
                dx, dy = xmax-xmin, ymax-ymin
                smax = max( abs(smax), abs(smin) )
                smin = -smax
                ds = smax-smin
                ax1.set_xlim( [ smin-0.2*ds, smax+0.2*ds ] )
                ax1.set_ylim( [ -1, 1 ] )
                ax2.set_xlim( [ xmin-dx, xmax+dx ] )
                ax2.set_ylim( [ ymin-dy, ymax+dx ] )
                ax1.axvline( 0, color='gray', lw=2 )
                ax1.text( 0.05, 0.5, 'bend apex', rotation=90 )
                ax1.text( smin-0.1*ds, 0.93, 'flow' )
                ax1.arrow( smin-0.1*ds, 0.9, 0.2*ds, 0 )
                ax2.legend( loc='center left', bbox_to_anchor=(1,0.5), ncol=2 )
                plt.show()


class FreeTemporalBars( TemporalBars ):

    '''
    FreeTemporalBars( TemporalBars )
    ================================

    A class to reconstruct the temporal dynamics of free migrating bare sediment bars
    '''
    
    def AccumulateBends( self ):
        
        '''For each bend we follow its history through indices'''

        BendIdx = self.Bars[0].unwrapper.Bend.astype(int)
        Bends = np.unique(BendIdx[BendIdx>=0]).astype(int)
        NextBend = self.Bars[0].unwrapper.NextBend.astype(int)
        self.BendAccumulator = -np.ones( (Bends.size,len(self.T)), dtype=int )
        for iBend, Bend in enumerate( Bends ):
            self.BendAccumulator[iBend,0] = Bend
            b = NextBend[ BendIdx==Bend ][0]
            for iFinder, Finder in enumerate( self.Bars[1:], 1 ):
                if b == -1: break
                self.BendAccumulator[iBend,iFinder] = int(b)
                b = Finder.unwrapper.NextBend[ Finder.unwrapper.Bend==b ][0]
        return self.BendAccumulator
        

    def CorrelateBars( self ):
        '''For each BarIdx(t) compute BarIdx(t+dt)'''

        accumulator = self.AccumulateBends()
        self.BarAccumulator = -np.ones( (self.Bars[0].BarIdx.size,len(self.T)), dtype=int )
        self.BarAccumulator[:,0] = self.Bars[0].BarIdx
        self.BarsCorr = []
        
        for iBars, (BarsL, BarsR, TL, TR) in enumerate( zip( self.Bars[:-1], self.Bars[1:], self.T[:-1], self.T[1:] ) ):

            BarCorr = []
            dT = TR - TL

            for iBarL, BarL in enumerate( BarsL.BarIdx ):

                # Bar Centroid (t)
                bL = BarsL.unwrapper.b.mean()
                IL, JL = BarsL.Centroid[0,iBarL], BarsL.Centroid[1,iBarL]
                XcL, YcL = BarsL.unwrapper.XC[IL,JL], BarsL.unwrapper.YC[IL,JL]
                NL = BarsL.unwrapper.N[JL] # Transverse Coordinate
                SL = BarsL.unwrapper.s[IL] # Longitudinal Coordinate
                BendL = BarsL.BBIdx[iBarL] # Bend to which Bar(t) belongs
                SL_0 = BarsL.unwrapper.s[ BarsL.unwrapper.Bend==BendL ][0] # Coordinate of Upstream Inflection Point
                RL = SL - SL_0 # Relative Longitudinal Coordinate

                # Bar Centroid (t+dt)
                bR = BarsR.unwrapper.b.mean()
                xR, yR = BarsR.unwrapper.XC[BarsR.Centroid[0,:], BarsR.Centroid[1,:]], BarsR.unwrapper.YC[BarsR.Centroid[0,:], BarsR.Centroid[1,:]]
                NR = BarsR.unwrapper.N[ BarsR.Centroid[1,:] ] # Transverse Coordinate
                SR = BarsR.unwrapper.s[ BarsR.Centroid[0,:] ] # Longitudinal Coordinate
                BendR = BarsL.unwrapper.NextBend[BarsL.unwrapper.Bend==BendL][0] # Relative to the same Bend of BarsL
                SR_0 = BarsR.unwrapper.s[ BarsR.unwrapper.Bend==BendR ][0] # Coordinate of Upstream Inflection Point
                RR = SR - SR_0 # Relative Longitudinal Coordinate

                if any(( BendL<0, BendR<0 )):
                    BarCorr.append( [iBarL, IL, JL, XcL, YcL, -1, -1, -1, np.nan, np.nan, np.nan, np.nan, np.nan] )
                    continue

                # Reference System (0)
                try:
                    b0 = self.Bars[0].unwrapper.b.mean()
                    Bend0 = self.Bars[0].unwrapper.Bend[ self.Bars[0].unwrapper.Bend==accumulator[accumulator[:,iBars]==BendL,0] ][0]
                    S_0 = self.Bars[0].unwrapper.s[ self.Bars[0].unwrapper.Bend==Bend0 ][0]
                    L0 = self.Bars[0].unwrapper.s[self.Bars[0].unwrapper.Bend==Bend0][-1] - S_0 # Bends Length (0)
                except IndexError:
                    BarCorr.append( [iBarL, IL, JL, XcL, YcL, -1, -1, -1, np.nan, np.nan, np.nan, np.nan, np.nan] )
                    continue

                lL = BarsL.unwrapper.s[BarsL.unwrapper.Bend==BendL][-1] - BarsL.unwrapper.s[BarsL.unwrapper.Bend==BendL][0] # Bend Length (t)
                lR = BarsR.unwrapper.s[BarsR.unwrapper.Bend==BendR][-1] - BarsR.unwrapper.s[BarsR.unwrapper.Bend==BendR][0] # Bend Length (t+dt)

                mask = np.logical_or.reduce(( np.abs(NR-NL)>0.25,
                                              (RR-RL)>2*bL*dT,
                                              (RR-RL)<-0.5*bL*dT,
                                              np.sqrt( (xR-XcL)**2 + (yR-YcL)**2 )>2*bL*dT
                                          ))
                RR[ mask ] = np.nan

                try:
                    iBarR = np.nanargmin( np.abs(RR-RL) )
                    BarR = BarsR.BarIdx[iBarR]
                    IR, JR = BarsR.Centroid[0,iBarR], BarsR.Centroid[1,iBarR]
                    if iBarL>0 and BarCorr[-1][5] == iBarR:
                        if (RR[iBarR]-RL) >= BarCorr[-1][10]:
                            raise ValueError
                        else:
                            BarCorr[-1][5:12] = [-1, -1, -1, np.nan, np.nan, np.nan, np.nan]

                    # Scale on Bend Elongation (if more than one bend is involved, we account for all of them)
                    if SR[iBarR] > BarsR.unwrapper.s[BarsR.unwrapper.Bend==BendR][-1]: ibend = 1 # The bar has moved to the next bend
                    else: ibend = 0
                    LL0 = self.Bars[0].unwrapper.s[self.Bars[0].unwrapper.Bend==(Bend0+ibend)][-1] - self.Bars[0].unwrapper.s[self.Bars[0].unwrapper.Bend==Bend0][0] # Bends Length (0)
                    LL = BarsL.unwrapper.s[BarsL.unwrapper.Bend==(BendL+ibend)][-1] - BarsL.unwrapper.s[BarsL.unwrapper.Bend==BendL][0] # Bends Length (t)
                    LR = BarsR.unwrapper.s[BarsR.unwrapper.Bend==(BendR+ibend)][-1] - BarsR.unwrapper.s[BarsR.unwrapper.Bend==BendR][0] # Bends Length (t+dt)
                    rR = (RR[iBarR]-RL) * LL/LR / BarsL.unwrapper.b.mean()
                    nR = (NR[iBarR]-NL) * LL/LR
                    BarCorr.append( [iBarL, IL, JL, XcL, YcL, iBarR, IR, JR, xR[iBarR], yR[iBarR], RR[iBarR]*L0/lR-RL*L0/lL, nR, RL*L0/lL+S_0] )

                    self.BarAccumulator[self.BarAccumulator[:,iBars]==BarL,iBars+1] = iBarR
                except ValueError:
                    BarCorr.append( [iBarL, IL, JL, XcL, YcL, -1, -1, -1, np.nan, np.nan, np.nan, np.nan, RL*L0/lL+S_0] )
                    continue

                if False: # Centroids Correlation
                    f = plt.figure()
                    ax1 = f.add_subplot( 211 )
                    ax2 = f.add_subplot( 212, sharex=ax1, sharey=ax1 )
                    ax1.set_title('%d' % (iBarL))
                    ax2.set_title('%d' % (iBarR))
                    p1 = ax1.pcolormesh( BarsL.unwrapper.XC, BarsL.unwrapper.YC, BarsL.Bars, cmap='spectral' )
                    ax1.plot( BarsL.unwrapper.XC[IL,JL], BarsL.unwrapper.YC[IL,JL], 'wo', markersize=16 )
                    p2 = ax2.pcolormesh( BarsR.unwrapper.XC, BarsR.unwrapper.YC, BarsR.Bars, cmap='spectral' )
                    ax2.plot( BarsR.unwrapper.XC[IR,JR], BarsR.unwrapper.YC[IR,JR], 'wo', markersize=16 )
                    ax1.set_aspect('equal')
                    ax2.set_aspect('equal')
                    plt.show()

            # Test Plots!!
            if False: # Bars Arrows
                f = plt.figure()
                plt.pcolormesh( BarsL.unwrapper.XC, BarsL.unwrapper.YC, BarsL.Bars, cmap='Spectral', alpha=0.5 )
                plt.contour( BarsR.unwrapper.XC, BarsR.unwrapper.YC, BarsR.Bars>0, 1, cmap='binary_r', linewidths=3 )
                for i in xrange( len(BarCorr) ):
                    if BarCorr[i][5]<0: continue
                    [ x0, y0 ] = BarCorr[i][3:5]
                    #[ x1, y1 ] = BarCorr[i][8:10]
                    I = int(BarCorr[i][1])
                    alpha = BarsL.unwrapper.theta[np.abs(BarsL.unwrapper.s-BarsL.unwrapper.s[I]).argmin()]
                    ds_i, dn_i = BarCorr[i][10], BarCorr[i][11]*BarsL.unwrapper.b.mean()
                    x1 = x0 + ds_i*np.cos(alpha) - dn_i*np.sin(alpha) #np.sqrt(ds_i**2+dn_i**2)*np.cos(alpha)
                    y1 = y0 - ds_i*np.sin(alpha) + dn_i*np.cos(alpha) #np.sqrt(ds_i**2+dn_i**2)*np.sin(alpha)
                    plt.plot( x0, y0, 'yo' )
                    plt.arrow( x0, y0, x1-x0, y1-y0, facecolor='k', edgecolor='k', head_width=50, head_length=50, width=30 )
                    plt.text( x0, y0, '%d' % int(BarCorr[i][0]), color='w' )
                plt.axis('equal')
                plt.show()                            

            self.BarsCorr.append( BarCorr )
        return self.BarsCorr

    
    def CentroidsEvol( self, bend_idx, normalize=True ):
        '''Follow the evolution of the centroid of main bar of an individual meander bend'''

        # We compute the Local Migration Rates and Wavenumbers of Channel Bars

        self.CorrelateBars() # Create Accumulator Matrix of Bars Correlation

        hwidths = NaNs( len(self.Bars)-1 ) # Reference Half-Width of the channel for each TimeFrame

        # Iterate over TimeFrames
        for iFinder, (T1, T2, Finder, BarCorr) in enumerate( zip( self.T[:-1], self.T[1:], self.Bars[:-1], self.BarsCorr ) ):
            dT = T2 - T1
            s, n = Finder.unwrapper.s, Finder.unwrapper.N
            NMAX = len( BarCorr )
            I, J = np.zeros( NMAX, dtype=int ), np.zeros( NMAX, dtype=int )
            dsi, dni, dxi, dyi, dmi, dzi = map( NaNs, [NMAX]*6 )
            hwidths[iFinder] = Finder.unwrapper.b.mean()
            for i in xrange( NMAX ): # For each Bar in the Time Frame
                I[i] = BarCorr[i][1] # Longitudinal Position of Bar Centroid
                J[i] = BarCorr[i][2] # Transversal Position of Bar Centroid
                dsi[i] = BarCorr[i][10] # Longitudinal Distance to which the Bar Centroid has Migrated
                dni[i] = BarCorr[i][11] # Transversal Distance to which the Bar Centroid has Migrated
                dzi[i] = np.sqrt( (dsi[i])**2 + (dni[i])**2 ) # Total Intrinsic Distance Migrated
                dxi[i] = BarCorr[i][8] - BarCorr[i][3] # Cartesian x Distance of Migration
                dyi[i] = BarCorr[i][9] - BarCorr[i][4] # Cartesian y Distance of Migration
                dmi[i] = np.sqrt( (dxi[i])**2 + (dyi[i])**2 ) # Total Cartesian Distance Migrated

            # Position of Channel Bars in Intrinsic and Cartesian Reference Systems
            si, ni, xi, yi = s[I]/Finder.unwrapper.b.mean(), n[J], Finder.unwrapper.XC[I,J], Finder.unwrapper.YC[I,J]

            if False: #True: 
                # Plot (X,Y) Bars Correlation with Arrows for the current TimeFrame
                plt.figure()
                mask = np.isfinite(dsi)
                Z = scipy_interp.griddata( (si[mask],ni[mask]), dsi[mask], (Finder.unwrapper.Sc, Finder.unwrapper.Nc), method='nearest' ).T /hwidths[iFinder]/dT
                plt.pcolor( Finder.unwrapper.XC, Finder.unwrapper.YC, np.ma.array(Z,mask=np.isnan(Z)), cmap='viridis', alpha=0.75 )
                plt.colorbar()
                plt.contour( Finder.unwrapper.XC, Finder.unwrapper.YC, Finder.Bars>0, 1, colors='g', linewidths=2 )
                plt.contour( self.Bars[iFinder+1].unwrapper.XC, self.Bars[iFinder+1].unwrapper.YC, self.Bars[iFinder+1].Bars>0, 1, colors='r', linewidths=2 )
                for i in xrange( len(BarCorr) ):
                    if BarCorr[i][5]<0: continue
                    [ x0, y0 ] = BarCorr[i][3:5]
                    #[ x1, y1 ] = BarCorr[i][8:10]
                    alpha = Finder.unwrapper.theta[np.abs(s-s[I[i]]).argmin()]
                    ds_i, dn_i = BarCorr[i][10], BarCorr[i][11]*hwidths[iFinder]
                    x1 = x0 + ds_i*np.cos(alpha) - dn_i*np.sin(alpha) #np.sqrt(ds_i**2+dn_i**2)*np.cos(alpha)
                    y1 = y0 - ds_i*np.sin(alpha) + dn_i*np.cos(alpha) #np.sqrt(ds_i**2+dn_i**2)*np.sin(alpha)
                    #plt.plot( x0, y0, 'ko', markersize=2 )
                    plt.arrow( x0, y0, x1-x0, y1-y0, facecolor='b', edgecolor='b', head_width=50, head_length=50, width=30 )
                    #plt.annotate('', xytext=(x0,y0), xy=(x1,y1), arrowprops=dict(arrowstyle='simple', edgecolor='w', facecolor='w') )
                    bidx = np.where(np.abs(np.ediff1d(Finder.unwrapper.Bend, to_begin=0)).astype(int) == 1)[0]
                    plt.plot( [ Finder.unwrapper.XC[bidx,0], Finder.unwrapper.XC[bidx,-1] ], [ Finder.unwrapper.YC[bidx,0],Finder.unwrapper.YC[bidx,-1] ], 'g', lw=2 )
                plt.axis('equal')
                #plt.show()

                # Plot (S,N) Bars Correlation with Arrows for the current TimeFrame
                plt.figure(figsize=(10.24, 2.56))
                plt.pcolor( Finder.unwrapper.Sc, Finder.unwrapper.Nc, np.ma.array(Z,mask=np.isnan(Z)).T, cmap='viridis' )
                plt.colorbar()
                plt.contour( Finder.unwrapper.Sc, Finder.unwrapper.Nc, Finder.Bars.T>0, 1, colors='r' )
                for i in xrange( len(BarCorr) ):
                    if BarCorr[i][5]<0: continue
                    [ s0, n0 ] = si[i], ni[i]
                    [ ds, dn ] = dsi[i]/hwidths[iFinder], dni[i]
                    #plt.plot( s0, n0, 'ko' )
                    plt.arrow( s0, n0, ds, dn, facecolor='b', edgecolor='b' )
                plt.show()

        Bavg = np.nanmean( hwidths ) # Time-Space Averaged Channel Width

        # Locate Bars and Define Migration Vectors on the Reference River Structure for ALL the TimeFrames
        Xi, Yi, Si, Ni, DSi, DNi, YYi, Wi, Bi = [], [], [], [], [], [], [], [], [] # Pos(S), Pos(N), Migr(S), Migr(N), Year Index, Wavenumber

        for cnt, (BarCorr, Finder, T1, T2) in enumerate( zip(self.BarsCorr, self.Bars[:-1], self.T[:-1], self.T[1:]) ): # For all the TimeFrames
            dT = (T2-T1)
            icnt = 0
            for i in xrange( len(BarCorr) ): # For all the Channel Bars in the Current TimeFrame
                if BarCorr[i][5]<0: continue
                icnt += 1
                Xi.append( BarCorr[i][3] ), Yi.append( BarCorr[i][4] )
                [ s0, n0 ] = BarCorr[i][12]/Bavg, Finder.unwrapper.N[BarCorr[i][2]] # Position (s,n)
                [ ds, dn ] = BarCorr[i][10]/Bavg/dT, BarCorr[i][11]/dT # Migration Vector (ds,dn)
                Si.append( s0 ), Ni.append( n0 )
                DSi.append( ds ), DNi.append( dn )
                YYi.append(cnt)
                Bi.append( self.Bars[0].unwrapper.b[np.abs(self.Bars[0].unwrapper.s-s0).argmin()] )
            # Channel Bar Wavenumbers for the Current TimeFrame
            b0 = hwidths[ cnt ] # Channel Width for the current TimeFrame
            wi = NaNs( icnt ) # Wavenumbers of the Current TimeFrame
            try:
                maskl, maskr = np.asarray(Ni)[-icnt:]>0, np.asarray(Ni)[-icnt:]<0 # Left/Right-Handed Channel Bars
                #wi[maskl] = 2*np.pi / np.gradient( np.asarray(Si)[-icnt:][maskl] ) # Wavenumbers of the Left Channel Bars
                #wi[maskr] = 2*np.pi / np.gradient( np.asarray(Si)[-icnt:][maskr] ) # Wavenumbers of the Right Channel Bars
                wi = np.pi / np.gradient( np.asarray(Si)[-icnt:] ) # Wavenumbers of the Left Channel Bars
                Wi += ( wi*np.asarray(Bi)[-icnt:]/Bavg ).tolist() # Correct Width - it may change significantly, therefore use local width!!!
            except ValueError:
                Wi += NaNs( icnt ).tolist()

        Xi, Yi, Si, Ni, DSi, DNi, YYi, Wi = map( np.asarray, (Xi, Yi, Si, Ni, DSi, DNi, YYi, Wi) )

        X, Y = self.Bars[0].unwrapper.XC, self.Bars[0].unwrapper.YC # Reference 
        S, N = self.Bars[0].unwrapper.Sc*self.Bars[0].unwrapper.b.mean()/Bavg, self.Bars[0].unwrapper.Nc
        
        # In order to apply the criteria to an evolving meandering river, we must 'freeze its shape',
        # i.e. we define a reference river structure according to its initial structure
        X, Y = self.Bars[0].unwrapper.XC, self.Bars[0].unwrapper.YC # Reference 
        S, N = self.Bars[0].unwrapper.Sc*self.Bars[0].unwrapper.b.mean()/Bavg, self.Bars[0].unwrapper.Nc
        Z = scipy_interp.griddata( (Si, Ni), DSi, (S,N), method='linear' ).T

        si, dsi, wi = [np.asarray(xx) for xx in zip(*sorted(zip(Si, DSi, Wi), key=lambda pair: pair[0]))]
        DDS = np.gradient( si )

        if False: #True: #False:
            plt.figure(figsize=(16,2))
            plt.title( 'Average annual longitudinal migration rate of channel bars (%d-%d)' % (int(self.T[0]), int(self.T[-1])) )
            plt.pcolormesh( S, N, np.ma.array(Z,mask=np.isnan(Z)).T, cmap='jet' )
            cs = [cm.RdYlBu_r(xx) for xx in np.linspace(0,1,len(self.BarsCorr))]
            for cnt, (BarCorr, Finder, T1, T2) in enumerate( zip(self.BarsCorr, self.Bars[:-1], self.T[:-1], self.T[1:]) ):
                dT = (T2 - T1)
                for i in xrange( len(BarCorr) ):
                    if BarCorr[i][5]<0: continue
                    #BarCorr.append( [iBarL, IL, JL, XcL, YcL, iBarR, IR, JR, xR[iBarR], yR[iBarR], rR, nR, RL*L0/lL+S_0] )
                    [ s0, n0 ] = BarCorr[i][12]/Bavg, Finder.unwrapper.N[BarCorr[i][2]]
                    [ ds, dn ] = BarCorr[i][10]/Bavg/dT, BarCorr[i][11]/dT
                    plt.plot( s0, n0, 'o', color=cs[cnt] )
                    plt.arrow( s0, n0, ds, dn, facecolor=cs[cnt], edgecolor=cs[cnt] )
            plt.xlabel( r'$x$' )
            plt.xlabel( r'$y$' )
            plt.colorbar()
            plt.show()
        if False: #True: #False:
            plt.figure()
            plt.title( 'Average annual longitudinal migration rate of channel bars (%d-%d)' % (int(self.T[0]), int(self.T[-1])) )
            plt.pcolormesh( X, Y, np.ma.array(Z,mask=np.isnan(Z)), cmap='jet' )
            plt.xlabel( r'$x$' )
            plt.xlabel( r'$y$' )
            plt.colorbar()
            plt.axis( 'equal' )
            plt.show()

        if False: #True: # False:
            plt.figure()
            si, ni, dsi, dni = [np.array(x) for x in zip(*sorted(zip(Si.tolist(), Ni.tolist(), DSi.tolist(), DNi.tolist()), key=lambda pair: pair[0]))]
            for ifilter in xrange(10): dsi[1:-1] * 0.25 * (dsi[:-2] + dsi[2:] + 2*dsi[1:-1])
            plt.plot(si, dsi, 'o')

            plt.show()

        return S, N, X, Y, Z, Xi, Yi, Si, Ni, DSi, DNi, YYi, Wi, BarCorr
