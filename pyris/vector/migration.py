# ==================================================
# Module: vector
# File: migration.py
# Package: PyRIS
# Author: Federico Monegaglia
# Date: April 2016
# Description: Migration vectors and bend separation
# ==================================================

from __future__ import division
import numpy as np
from scipy import interpolate
from ..misc import NaNs, Intersection
from .interpolation import InterpPCS, CurvaturePCS
# import matplotlib.pyplot as plt # check if it works, if not add also cm

class AxisMigration( object ):

    '''
    AxisMigration( object )
    =======================

    A callable class for computing Bend Separation and Migration Vectors

    Arguments
    ---------
    Xseries      Temporal sequence (list) of x coordinates of the river centerline
    Yseries      Temporal sequence (list) of y coordinates of the river centerline

    Returns
    -------
    dx           Temporal sequence (list) of x component of the centerline migration vectors
    dy           Temporal sequence (list) of y component of the centerline migration vectors
    dz           Temporal sequence (list) of the magnitude of the centerline migration vectors
    Css          Temporal sequence (list) of the smoothed planform curvature
    BI           Temporal sequence (list) of the bed indexes
    B12          Temporal sequence (list) of the next bed indexes (where the current bend goes)
    BUD          Temporal sequence (list) of the bend upstream-downstream index
    '''

    Css = []
    method = 'curvature' # distance must be taken away
    
    def __init__( self, Xseries, Yseries ):
        '''
        Get Subsequent Channel Planforms
        '''
        self.data = []
        for x, y in zip( Xseries, Yseries ):
            x, y = np.asarray(x), np.asarray(y)
            dx = np.ediff1d( x, to_begin=0 )
            dy = np.ediff1d( y, to_begin=0 )
            ds = np.sqrt( dx**2 + dy**2 )
            s = np.cumsum( ds )
            theta = np.arctan2( dy, dx )
            rtheta = theta.copy()
            for i in range( 1, theta.size ): # Set theta continuous
                if (theta[i]-theta[i-1])/np.pi > +1.9: theta[i] -=2*np.pi
                if (theta[i]-theta[i-1])/np.pi < -1.9: theta[i] +=2*np.pi
            c = -np.gradient( theta ) / np.gradient( s )
            self.data.append( { 'x': x, 'y': y, 's': s, 'c':c, 't':theta, 'r':rtheta } )
        return None

    def IterData( self ):
        '''Planform data iterator'''
        for i, d in enumerate( self.data ): yield i, d

    def IterData2( self ):
        '''Planform data pair iterator'''
        for i, (d1, d2) in enumerate( zip( self.data[:-1], self.data[1:]) ): yield i, (d1, d2)

    def RevIterData2( self ):
        '''Reverse planform data pair iterator'''
        for i, (d1, d2) in enumerate( zip( self.data[:-1][::-1], self.data[1:][::-1] ) ): yield ( len(self.data)-2-i ), (d2, d1)

    def Iterbends( self, Idx ):
        '''Bend Iterator'''
        for i, (il,ir) in enumerate( zip( Idx[:-1], Idx[1:] ) ): yield i, (il, ir)

    def Iterbends2( self, Idx1, Idx2 ):
        '''Bend Pair Iterator'''
        for i, (i1l,i1r,i2l,i2r) in enumerate( zip( Idx1[:-1], Idx1[1:], Idx2[:-1], Idx2[1:] ) ): yield i, (i1l,i1r,i2l,i2r)

    def FilterAll( self, pfreq=10 ):
        '''Perform ICWT Filtering on all the Data'''
        for i, d in self.IterData():
            if i==0:
                self.Css.append( self.FilterCs( d['c'], d['s'], d['x'], d['y'], pfreq=int(2.5*pfreq) ) )
            else:
                self.Css.append( self.FilterCs( d['c'], d['s'], d['x'], d['y'], pfreq=int(pfreq) ) )
        return None

    def FilterCs( self, s, Cs, x, y, pfreq ):
        '''Use PCS interpolation to filter out curvature signal'''
        xp_PCS, yp_PCS, d1xp_PCS, d1yp_PCS, d2xp_PCS, d2yp_PCS = InterpPCS(
            x[::pfreq], y[::pfreq], N=x.size, s=x.size+np.sqrt(2*x.size),
            with_derivatives=True, k=2 )
        sp_PCS, thetap_PCS, Cs_PCS = CurvaturePCS(
            xp_PCS, yp_PCS, d1xp_PCS, d1yp_PCS, d2xp_PCS,
            d2yp_PCS, method=2, return_diff=False )
        return Cs_PCS
        
    def GetInflections( self, Cs ):
        '''Compute 0-crossings of channel curvature'''
        IDX = np.where( Cs[1:]*Cs[:-1] < 0 )[0]
        return IDX

    def DistanceInflections( self, data, prev_data, prev_I ):
        '''Compute Inflection points by moving orthogonally from previous inflection points'''
        I = NaNs( prev_I.size )
        for j, prev_i in enumerate(prev_I):
            if np.isnan( prev_i ): continue
            i = self.FindOrthogonalPoint( data, prev_data, prev_i )
            if j > 0 and i <= I[j-1]: continue
            I[j] = i if i is not None else np.nan
        return np.asarray( I )            

    def GetAllInflections( self ):
        '''Get Inflection points on Filtered Curvature.'''
        self.I = []
        if self.method == 'curvature':
            for i, d in self.IterData():
                self.I.append( self.GetInflections( self.Css[i] ) )
        elif self.method == 'distance':
            for i, d in self.IterData():
                if i == 0: self.I.append( self.GetInflections( self.Css[i] ) )
                else:
                    self.I.append( self.DistanceInflections( d, self.data[i-1], self.I[i-1] ) )
        return None

    def CorrelateInflections( self, *args, **kwargs ):
        '''Find the Closest Inflection Points'''

        self.CI1 = [] # Points on the Current Planform
        self.CI12 = [] # Points to which the First Planform Points Converge to the Second Planform

        if self.method == 'distance':            
            self.CI1 = [ [] for _ in range( len( self.data ) ) ]
            self.CI12 = [ [] for _ in range( len( self.data ) ) ]
            for i, (d1, d2) in self.IterData2():
                mask = np.isfinite( self.I[i+1] )
                self.CI1[i+1] = self.I[i+1][ mask ].astype( int )
                self.CI1[i] = self.I[i][ mask ].astype( int )
                self.CI12[i] = self.CI1[i+1]                
            self.CI12[-1] = self.CI1[-1]
            return None

        elif self.method == 'curvature':
            C1 = self.I[0] # Initial Reference Planform
            # Go Forward
            for i, (d1, d2) in self.IterData2():
                C2 = self.I[i+1]
                C12 = np.zeros_like( C1, dtype=int )
                x1, y1, R1, T1 = d1['x'], d1['y'], d1['r'], d1['t']
                x2, y2, R2, T2 = d2['x'], d2['y'], d2['r'], d2['t']

                for ipoint, Ipoint in enumerate( C1 ):
                    xi1, yi1, ti1 = x1[Ipoint], y1[Ipoint], R1[Ipoint]
                    xC2, yC2, tC2 = x2[C2], y2[C2], R2[C2] # Do not care about sign
                    mask = np.logical_and( abs( (R2[C2+1])-(R1[Ipoint+1]) )>0.5*np.pi, abs( (R2[C2+1])-(R1[Ipoint+1]) )<1.5*np.pi ) # Data to be masked with NaNs
                    xC2 = np.where( mask, np.nan, x2[C2] )
                    yC2 = np.where( mask, np.nan, y2[C2] )
                    tC2 = np.where( mask, np.nan, R2[C2] )
                    # Find the Closest
                    try:
                        C12[ipoint] = C2[ np.nanargmin( np.sqrt( (xC2-xi1)**2 + (yC2-yi1)**2 ) ) ]
                    except ValueError:
                        raise ValueError('not able to compute inflection correlation for planform n. %d. Please check your axis' % (i+2))
                    
                # There are some duplicated points - we need to get rid of them
                unique, counts = np.unique(C12, return_counts=True)
                duplic = unique[ counts>1 ]
                cduplic = counts[ counts > 1 ]                
                for idup, (dup, cdup) in enumerate( zip( duplic, cduplic ) ):
                    idxs = np.where( C12==dup )[0]
                    idx = np.argmin( np.sqrt( (x2[dup]-x1[C1][idxs])**2 + (y2[dup]-y1[C1][idxs])**2 ) )                    
                    idxs = np.delete( idxs, idx )
                    C1 = np.delete( C1, idxs )
                    C12 = np.delete( C12, idxs )
                    
                while np.any(C12[:-1]>C12[1:]):
                    for j in range( 1, C1.size ):
                        l, r = j-1, j
                        cl, cr = C12[l], C12[r]
                        if cr < cl:
                            if np.sqrt( (x2[cl]-x1[C1[l]])**2 + (y2[cl]-y1[C1[l]])**2 ) < np.sqrt( (x2[cr]-x1[C1[l]])**2 + (y2[cr]-y1[C1[l]])**2 ):
                                C1, C12 = np.delete( C1, r ), np.delete( C12, r )
                            else:
                                C1, C12 = np.delete( C1, l ), np.delete( C12, l )
                            break

                self.CI1.append(C1)
                self.CI12.append(C12)
                C1 = C12
        self.CI1.append(C12)
        return None

    def BendUpstreamDownstream( self, I, Cs ):
        '''Bend Upstream-Downstream Indexes'''
        BUD = NaNs( Cs.size )
        for i, (il,ir) in self.Iterbends( I ):
            iapex = il + np.abs( Cs[ il:ir ] ).argmax()
            BUD[ il ] = 2 # Inflection Point
            BUD[ ir ] = 2 # Inflection Point
            BUD[ iapex ] = 0 # Bend Apex
            BUD[ il+1:iapex ] = -1  # Bend Upstream
            BUD[ iapex+1:ir ] = +1 # Bend Downstream
        return BUD

    def AllBUDs( self ):
        '''Bend Upstream-Downstream Indexes for All Planforms'''
        self.BUD = []
        for i, d in self.IterData():
            self.BUD.append( self.BendUpstreamDownstream( self.CI1[i], self.Css[i] ) )
        return None

    def GetBends( self, c ):
        '''Returns Inflection Points, Bend Indexes'''
        Idx = self.GetInflections( c )
        BIDX = self.LabelBends( c.size, Idx )
        return BIDX, Idx

    def LabelBends( self, *args, **kwargs ):
        'Bend label for each point of the planform'
        N = args[0] if isinstance( args[0], int ) else args[0].size
        Idx = args[1]
        labels = -np.ones( N, dtype=int )
        for i, (il, ir) in self.Iterbends( Idx ):
            labels[il:ir] = i
        return labels

    def LabelAllBends( self ):
        '''Apply Bend Labels to Each Planform'''
        self.BI = []
        for i, d in self.IterData():
            self.BI.append( self.LabelBends( d['s'].size, self.CI1[i] ) )
        return None

    def CorrelateBends( self ):
        '''Once Bends are Separated and Labeled, Correlate Them'''
        self.B12 = []
        self.ignore = []
        for di, (d1, d2) in self.IterData2():
            B1 = self.BI[di]
            B2 = self.BI[di+1]
            B12 = -np.ones( B1.size, dtype=int )
            I1 = self.CI1[di]
            I2 = self.CI1[di+1]
            I12 = self.CI12[di]
            x1, y1, c1 = d1['x'], d1['y'], d1['c']
            x2, y2, c2 = d2['x'], d2['y'], d2['c']

            # X il momento tengo la correlazione tra gli inflections
            for i, (i1l, i1r, i2l, i2r) in self.Iterbends2( I1, I12 ):
                vals, cnts = np.unique( B2[i2l:i2r], return_counts=True )
                if len( vals ) == 0:
                    B12[i1l:i1r] = -1
                else:
                    B12[i1l:i1r] = vals[ cnts.argmax() ]

            self.B12.append( B12 )
        self.B12.append( -np.ones( x2.size ) ) # Add a Convenience -1 Array for the Last Planform
        return None        

    def FindOrthogonalPoint( self, data1, data2, i2, L=None ):
        '''Find the orthogonal point to second line on the first one'''
        [ x1, y1, s1 ] = data1['x'], data1['y'], data1['s']
        [ x2, y2, s2 ] = data2['x'], data2['y'], data2['s']
        if L is None: L = 10*np.gradient( s1 ).mean()
        a0 = np.arctan2( ( y2[i2+1] - y2[i2-1] ), ( x2[i2+1] - x2[i2-1] ) )
        a = a0 - np.pi/2 # Local Perpendicular Angle
        P = np.array( [ x2[i2], y2[i2] ] )
        R = np.array( [ np.cos(a), np.sin(a) ] ) * L
        hits = []
        for i in range( 1, x1.size ):
            Q = np.array( [ x1[i-1], y1[i-1] ] )
            S = np.array( [ x1[i], y1[i] ] ) - Q
            # Bound angle
            a1 = np.arctan2( (y1[i]-y1[i-1]), (x1[i]-x1[i-1]) )
            if ( a0 > +np.pi/2 and a1 < -np.pi/2 ): a1 += 2*np.pi
            if ( a0 < -np.pi/2 and a1 > +np.pi/2 ): a1 -= 2*np.pi
            if a1 > a0+np.pi/4 or a1 < a0-np.pi/4: continue

            segments_intersect, (xi, yi) = Intersection( P, Q, R, S )
            if segments_intersect: hits.append( np.sqrt( (x1-xi)**2 + (y1-yi)**2 ).argmin() )
        
        if hits == []: return None
        return np.min( hits )

    def MigrationRates( self, data1, data2, I1, I12, B1, B2, B12 ):
        '''Compute Local Migration Rates by connected individual bends'''

        [ x1, y1, s1 ] = data1['x'], data1['y'], data1['s']
        [ x2, y2, s2 ] = data2['x'], data2['y'], data2['s']
        [ dx, dy, dz]  = [ NaNs( x1.size ), NaNs( x1.size ), NaNs( x1.size ) ]

        for i, (il,ir) in self.Iterbends( I1 ):
            # Isolate Bend

            if B12[il] < 0: continue # Bend Is not Correlated
            if B12[il-1] == B12[il]: continue # More bends go into one
            if B12[ir+1] == B12[il]: continue # More bends go into one

            mask1 = np.full( s1.size, False, dtype=bool ); mask1[il:ir]=True
            mask2 = B2==B12[il]
            bx1, by1, bs1, N1 = x1[mask1], y1[mask1], s1[mask1], mask1.sum() # Bend in First Planform
            bx2, by2, bs2, N2 = x2[mask2], y2[mask2], s2[mask2], mask2.sum() # Bend in Second Planform
            if N1<=1 or N2<=1: continue
            if N2 > N1: # Remove Random Points from Second Bend in order to interpolate
                idx = np.full( N2, True, bool )
                idx[ np.random.choice( np.arange(2,N2-2), N2-N1, replace=False ) ] = False
                bx2 = bx2[ idx ]
                by2 = by2[ idx ]
                N2 = bx2.size
            # ReInterpolate Second Planform (Parametric Cubic Spline)
            if N1 <= 1 or N2 <= 1: continue
            if N1 <= 3 or N2 <= 3: kpcs=1 # If we have too few points, use linear interpolation
            elif N1 <= 5 or N2 <= 5: kpcs=3 # If we have too few points, use linear interpolation
            else: kpcs=5
            bx2, by2 = InterpPCS( bx2, by2, N=N1, s=N2, k=kpcs, with_derivatives=False )
            # Compute Migration Rates for the whole bend
            dxb = bx2 - bx1
            dyb = by2 - by1
            dzb = np.sqrt( dxb**2 + dyb**2 )

            # Sinuosity Control
            sigma2 = ( bs2[-1] - bs2[0] ) / np.sqrt( (by2[-1]-by2[0])**2 + (bx2[-1]-bx2[0])**2 )
            sigma1 = ( bs1[-1] - bs1[0] ) / np.sqrt( (by1[-1]-by1[0])**2 + (bx1[-1]-bx1[0])**2 )
            # If Sinuosity has decreased significantly, assume a CutOff occurred
            if sigma1/sigma2 > 1.2: dxb, dyb, dzb = NaNs( N1 ), NaNs( N1 ), NaNs( N1 )
            # Set Migration Rate into Main Arrays
            dx[ mask1 ] = dxb
            dy[ mask1 ] = dyb
            dz[ mask1 ] = dzb
        return dx, dy, dz

    def AllMigrationRates( self ):
        '''Apply Migration Rates Algorithm to the whole set of planforms'''
        self.dx = []
        self.dy = []
        self.dz = []
        
        for i, (d1, d2) in self.IterData2():
            I1, I12 = self.CI1[i], self.CI12[i]
            B1, B2, B12 = self.BI[i], self.BI[i+1], self.B12[i]
            dxi, dyi, dzi = self.MigrationRates( d1, d2, I1, I12, B1, B2, B12 )
            self.dx.append( dxi ), self.dy.append( dyi ), self.dz.append( dzi )
        N = ( d2['s'] ).size
        self.dx.append( NaNs( N ) ), self.dy.append( NaNs( N ) ), self.dz.append( NaNs( N ) )
        return None
    
    def __call__( self, pfreq=10 ):
        self.FilterAll( pfreq=pfreq )
        self.GetAllInflections()
        self.CorrelateInflections()
        self.LabelAllBends()
        self.AllBUDs()
        self.CorrelateBends()
        self.AllMigrationRates()
        return self.dx, self.dy, self.dz, self.Css, self.BI, self.B12, self.BUD

