from __future__ import division
import numpy as np
from numpy.lib import stride_tricks
from skimage.measure import regionprops

class Pruner( object ):
    '''
    Spurs Removal Class    
    '''

    # Primitives
    primitives = ( np.array([[0,0,0],
                             [1,1,0],
                             [0,0,0]]),
                   np.array([[1,0,0],
                             [0,1,0],
                             [0,0,0]]),
                   np.array([[1,1,1],
                             [0,1,0],
                             [0,0,0]]),
                   np.array([[1,1,0],
                             [0,1,0],
                             [0,0,0]]),
                   np.array([[1,0,0],
                             [1,1,0],
                             [0,0,0]]) )

    wrong_blocks = ( np.array([[0,0,0],
                              [1,0,1],
                              [0,0,0]]),
                     np.array([[0,1,0],
                               [1,0,1],
                               [0,1,0]]) )
    
    smooth_block = np.array([[0,1,0],
                             [1,0,1],
                             [0,0,0]])
                        
    def __init__( self, img ):
        # Borders are automatically removed
        self.values = img
        self.img = np.where( img>0, 1, 0 )
        self.img[:2,:] = 0
        self.img[-2:,:] = 0
        self.img[:,:2] = 0
        self.img[:,-2:] = 0
        self.original = self.img.copy()
        # Create Bounding Box
        rp = regionprops( self.img )
        [xl, yl, xr, yr] = [int(b) for b in rp[0].bbox]
        xl -= 1
        yl -= 1
        xr += 1
        yr += 1
        self.bbox = ( slice( xl,xr ), slice( yl,yr ) )
        self.img = self.original[self.bbox]
        self.values = self.values[self.bbox]
                
    def BuildStrides( self ):
        '''Build cache-friendly Strides Array'''
        n = 3
        i = 1 + self.img.shape[0] - 3
        j = 1 + self.img.shape[1] - 3
        self.strides = stride_tricks.as_strided( self.img, (i,j,n,n),
                                                 strides=2*self.img.strides )

    def ValuesStrides( self ):
        '''Build cache-friendly Strides Array'''
        n = 3
        i = 1 + self.values.shape[0] - 3
        j = 1 + self.values.shape[1] - 3
        self.values_strides = stride_tricks.as_strided( self.values, (i,j,n,n),
                                                 strides=2*self.values.strides )
        
    def Update( self, seed_elems ):
        '''Update(i) - Prune for the current iteration i'''
        remainder = 0 # Convergence check
        for primitive in self.primitives:
            for side in range( 4 ):
                seed = np.rot90( primitive, side )
                mask = ( self.strides == seed ).all( axis=(2,3) )
                remainder += np.count_nonzero( mask )
                self.img[1:-1,1:-1] = np.where( mask, 0, self.img[1:-1,1:-1] )
        return remainder
    
    def FixWrongSkelBlocks( self ):
        '''Here we Fix Skeletonization Process'''
        self.BuildStrides()
        self.ValuesStrides()
        values_means = np.nanmean( np.where( self.values_strides==0, np.nan, self.values_strides ), axis=(2,3) )
        for iw, wrong_block in enumerate( self.wrong_blocks ):
            if iw == 1: Nsides=1 # Only one is required since this block is symmetric
            else: Nsides=4
            for side in range( Nsides ):
                seed = np.rot90( wrong_block, side )
                mask = ( self.strides == seed ).all( axis=(2,3) )
                self.img[1:-1,1:-1] = np.where( mask, 1, self.img[1:-1,1:-1] )
                self.values[1:-1,1:-1] = np.where( mask, values_means, self.values[1:-1,1:-1] )
            
    def Smooth( self, NMAX ):
        '''Smooth() - Remove Very High Local Curvature Gradients'''
        self.BuildStrides()
        print('smoothing centerline...')
        for side in range( 4 ):
            seed = np.rot90( self.smooth_block, side )
            mask = ( self.strides == seed ).all( axis=(2,3) )
            self.img[1:-1,1:-1] = np.where( mask, 1, self.img[1:-1,1:-1] )
        self.Prune( NMAX, True )        
                
    def Convergence( self, remainder ):
        '''
        Convergence( remainder ) - check remainder in order to state convergence
        '''
        if remainder <= 2: return True
        return False
                        
    def Prune( self, NMAX, verbose=True ):
        '''Prune( NMAX ) - Pruning Iterations'''
        for i in range( NMAX ):
            self.BuildStrides()
            if verbose: print('pruning iteration %03d on %3d' % ( i+1, NMAX ))
            remainder = self.Update( self.primitives )
            if self.Convergence( remainder ):
                if verbose:
                    print('pruning convergence reached after %3d iterations' % (i+1))
                break
            if i+1 == NMAX:
                print('maximum number of pruning iterations reached')

    def __call__( self, NMAX, verbose=True, fix_skel=True, smooth=True ):
        if fix_skel: self.FixWrongSkelBlocks()
        self.Prune( NMAX, verbose )
        if smooth: self.Smooth( NMAX )
        pruned_values = self.img * self.values
        ret = np.zeros( self.original.shape )
        ret[self.bbox] = pruned_values
        return ret.astype( int )



def Pruning( raster, NMAX=100, fix_skel=True, smooth=True, verbose=True ):
    
    '''
    Pruning( raster, maxiter=100 ) - Remove Spurs from skeletonized image.
    NMAX = maximum number of pruning iterations.
    Convenience Function for Class Pruner.
    '''
    img = raster.astype( int )
    pruner = Pruner( img )
    return pruner( NMAX, verbose, fix_skel, smooth )
