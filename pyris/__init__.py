# -*- coding: utf-8 -*-
# ===================================================================================
#
# Package: PyRIS
# Author: Federico Monegaglia
# Date: 2016
# Description: A Toolbox for extracting river planform features from satellite images
#
# ===================================================================================

'''
PyRIS :: Python - RIvers from Satellite
====================================
'''

# Imports
from __future__ import division
import os, sys
import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib import cm
from skimage import morphology as mm
from skimage.util import img_as_ubyte
from skimage.io import imread
from scipy import ndimage
from skimage.measure import regionprops
import warnings
import datetime

# Suppress Warnings
warnings.filterwarnings("ignore")

__author__ = 'Federico Monegaglia'
__email__ = 'f.monegaglia@gmail.com'
__version__ = '1.0'
__year__ = '2016'

__all__ = [
    # modules
    'raster', 'vector', 'misc',
    # pkg dataFcm
    '__author__', '__version__', '__email__', '__year__',
    # from standard packages
    'img_as_ubyte', 'imread', 'ndimage',
    # misc
    'GeoReference', 'NaNs', 'MaskClean',
    # raster
    'CleanIslands', 'RemoveSmallObjects', 'Skeletonize',
    'Pruner', 'Pruning',
    'Thresholding', 'SegmentationIndex',
    'Unwrapper', 'BarFinder', 'TemporalBars', 'FreeTemporalBars',
    # vector
    'AxisReader', 'ReadAxisLine',
    'InterpPCS', 'CurvaturePCS', 'WidthPCS',
    'AxisMigration', 'LoadLandsatData',
    ]

# Import Everything from SubModules
from .raster import *
from .vector import *
from .misc import *
from .config import *

def load( fname, *args, **kwargs ):
    '''Load file depending on the extension'''
    ext = os.path.splitext( fname )[-1]
    if ext == '.txt':
        return np.loadtxt( fname, *args, **kwargs )
    elif ext == '.npy':
        return np.load( fname, *args, **kwargs )
    else:
        e = 'Format %s not supported for file %s. Use either "npy" or "txt"' % ( ext, fname )
        raise TypeError(e)

def save( fname, *args, **kwargs ):
    '''Save file depending on the extension'''
    ext = os.path.splitext( fname )[-1]
    if ext == '.txt':
        return np.savetxt( fname, *args, **kwargs )
    elif ext == '.npy':
        return np.save( fname, *args, **kwargs )
    else:
        e = 'Format %s not supported for file %s. Use either "npy" or "txt"' % ( ext, fname )
        raise TypeError(e)

def get_year_jday( landsatname ):
    '''Get year and julian day from Landsat name'''
    collection = len(landsatname.split('_'))==7
    if collection: # Landsat Collection
        # Must Convert into Julian Calendar
        lname_arr = landsatname.split( '_' )
        year = lname_arr[3][:4]
        month = lname_arr[3][4:6]
        day = lname_arr[3][6:8]
        date_fmt = '%Y.%m.%d'
        date_s = '%s.%s.%s' % ( year, month, day )
        date = datetime.datetime.strptime( date_s, date_fmt )
        jday = '%03d' % date.timetuple().tm_yday
    else: # PreCollection
        # Already in Julian Calendar
        year = landsatname[9:13]
        jday = landsatname[13:16]
    # Return As Strings
    return year, jday

    
# ============================
# Functions of the Main Script 
# ============================

def segment_all( landsat_dirs, geodir, config, maskdir, auto_label=None ):
    '''
    segment_all( landsat_dirs, geodir, config, maskdir, auto_label=None )
    =====================================================================

    Iterate over all the Landsat Directories and perform image segmentation

    Arguments
    ---------
    landsat_dirs      directory containing all the landsat directories
    geodir            directory where GeoTransf instances are stored
    config            PyRIS' RawConfigParser instance
    maskdir           directory where channel masks are stored
    auto_label        mask selection method if more than one object occur (default None)

    Returns
    -------
    None
    '''
    to_skip = []
    # Iterate over Landsat Directories
    yearday = [ ''.join( get_year_jday(os.path.split(landsat_dir)[-1]) ) for landsat_dir in landsat_dirs ]
    yd, ldirs = [list(x) for x in zip(*sorted(zip(yearday, landsat_dirs), key=lambda pair: pair[0]))] # Sort by time
    for landsat in ldirs:
        # input
        landsatname = os.path.basename( landsat )
        year, jday = get_year_jday( landsatname )
        name = '_'.join( ( year, jday ) )
        # output
        maskfile = os.path.join( maskdir, '.'.join( (name, 'npy') ) )
        geofile = os.path.join( geodir, '.'.join( (name, 'p') ) )

        # skip the files which have already been processes
        if all( map( os.path.isfile, [ maskfile, geofile ] ) ):
            print
            print('data found for file %s - skipping '  % ( landsatname ))
            to_skip.append( name )
            continue

        print
        print('Processing file %s' % ( landsatname ))

        bands, GeoTransf = LoadLandsatData( landsat )

        print('applying BW masks...')

        # GeoReferencing of White and Black masks
        bw_masks_georef = GeoReference( GeoTransf )

        # Apply Black Mask
        black = np.ones( bands[0].shape, dtype=int )
        black_masks = eval( config.get( 'Segmentation', 'black_masks' ) )
        for s in black_masks:
            xx, yy = bw_masks_georef.RefCurve( np.asarray(s[2:]), np.asarray(s[:2]), inverse=True )
            sy, sx = slice( max(0, int(xx[0])), min(black.shape[1],int(xx[1])) ), slice( max(0,int(yy[0])), min(int(yy[1]),black.shape[0]) )
            black[ sx, sy ] = 0

        [ R, G, B, NIR, MIR, SWIR ] =  [ np.where(np.logical_or(band==0,black==0), np.nan, band) for band in bands ]

        # Set Dimensions
        pixel_width = config.getfloat('Data', 'channel_width') / GeoTransf['PixelSize'] # Channel width in Pixels
        radius = 2 * pixel_width # Circle Radius for Local Thresholding
        
        # Compute Mask
        print('computing mask...')
        _, mask, _ = SegmentationIndex( R=R, G=G, B=B, NIR=NIR, MIR=MIR, SWIR=SWIR, index=config.get('Segmentation', 'method'), 
                                        rad=radius, method=config.get('Segmentation', 'thresholding') )
        # Mask Landsat NoData
        print('nodata dilation...')
        nanmask = np.where( bands[0]==0, 1, 0 )
        nanmask = mm.binary_dilation( nanmask, mm.disk( 30 ) )
        mask = np.where( nanmask==1, 0, mask*black )

        # Image Cleaning
        print('cleaning mask...')
        mask = RemoveSmallObjects( mask, 100*pixel_width**2 ) # One Hundred Widths of Channel at Least is Required
        radius = max( np.floor( 0.5 * ( pixel_width ) ) - 3, 0 )
        mask = mm.binary_opening( mask, mm.disk( radius ) ) # Break Small Connectins
        mask = CleanIslands( mask, 10*pixel_width**2 ) # Clean Holes Inside the Planform
        mask = RemoveSmallObjects( mask, 100*pixel_width**2 ) # Remove New Small Objects
        mask = mask.astype( int )

        # Label Masks - we need to perform a rotation in order to have labels going from the largest to the smallest
        rot_angle = { 'b': 0, 'l': 1, 't': 2, 'r': 3 } # Counter-clockwise rotationan angle
        mask = np.rot90( mask, rot_angle[config.get( 'Data', 'flow_from' )] )
        mask_lab, num_features = ndimage.measurements.label( mask )
        # Rotate back to the original
        mask = np.rot90( mask, -rot_angle[config.get( 'Data', 'flow_from' )] )
        mask_lab = np.rot90( mask_lab, -rot_angle[config.get( 'Data', 'flow_from' )] )
        print('labelling feature in channel mask...')
        print('found %s features in river mask %s...' % ( num_features, os.path.basename(maskfile) ))

        if auto_label is None:
            plt.figure()
            plt.imshow( mask_lab, cmap=cm.nipy_spectral, interpolation='none' )
            plt.title( 'Indentify the label(s) corresponding to the river planform.' )
            for ifeat in range( 1, num_features+1 ):
                c0 = np.column_stack( np.where( mask_lab==ifeat ) )[0]
                plt.text( c0[1], c0[0], '%s' % ifeat, fontsize=30, bbox=dict( facecolor='white' ) )
            plt.show()
            labs = input( 'Please enter the label(s) do you want to use? (if more than one, separate them with a space): ' ).split(' ')
            mask *= 0
            for ilab, lab in enumerate( labs ): mask += np.where( mask_lab==int(lab), ilab+1, 0 )
        else:
            # The largest element in the image will be used.
            warnings.warn( 'automated labels may lead to erroneous planforms! please check your results!' )
            if auto_label == 'auto':
                if config.get( 'Segmentation', 'thresholding' ) == 'local': auto_label = 'max'
                elif config.get( 'Segmentation', 'thresholding' ) == 'global': auto_label = 'all'
            if auto_label == 'all':
                mask = mask_lab
            elif auto_label == 'max':
                labs = np.unique( mask_lab[mask_lab>0] )
                areas = np.zeros( labs.size )
                for ilab, lab in enumerate( labs ):
                    rp = regionprops( (mask_lab==lab).astype(int) )
                    [xl, yl, xr, yr] = [ int(b) for b in rp[0].bbox ]
                    areas[ilab] = abs( (xr-xl)*(yr-yl) )
                mask = mask_lab==labs[ areas.argmax() ]
            else:
                e = "labelling method '%s' not known. choose either 'auto', 'max', 'all' or None" % auto_label
                raise ValueError(e)

        print('saving  mask and GeoTransf data...')
        np.save( maskfile, mask )
        with open( geofile, 'wb' ) as gf: pickle.dump( GeoTransf, gf )
    return None


def import_gee_mask(config, geedir, geodir, maskdir ):
    '''
    import_gee_mask(geedir, geodir, maskdir )
    ===========================================

    Import .tif masks generated with gee outside PyRIS
     
    Arguments
    ---------
    geedir            directory containing all the mask files
    geodir            directory where GeoTransf instances are stored (default None)
    maskdir           directory containing all the mask files

    
    '''
    to_skip = []
    # Iterate over gee directory
    geemasks = sorted(os.listdir(geedir))
    for geemask in geemasks:
        # input
        geename, __ = os.path.splitext(geemask)
        year = geename[-12:-8]
        jday = str(213) # TODO CHANGE THIS WITH REAL DAY
        name = '_'.join( ( year, jday ) )
        # output
        maskfile = os.path.join( maskdir, '.'.join( (name, 'npy') ) )
        geofile = os.path.join( geodir, '.'.join( (name, 'p') ) )

        # skip the files which have already been processes
        if all( map( os.path.isfile, [ maskfile, geofile ] ) ):
            print
            print('data found for file %s - skipping '  % ( geename ))
            to_skip.append( name )
            continue

        print
        print('Loading %s' % ( geename ))

        # Loading geemask and georeferencing data
        mask, GeoTransf = LoadGeeMask( os.path.join( geedir, geemask ) )

        print('applying BW masks...')

        # # GeoReferencing of White and Black masks
        # bw_masks_georef = GeoReference( GeoTransf )

        # # Apply Black Mask
        # black = np.ones( mask.shape, dtype=int )
        # black_masks = eval( config.get( 'Segmentation', 'black_masks' ) )
        # for s in black_masks:
        #     xx, yy = bw_masks_georef.RefCurve( np.asarray(s[2:]), np.asarray(s[:2]), inverse=True )
        #     sy, sx = slice( max(0, int(xx[0])), min(black.shape[1],int(xx[1])) ), slice( max(0,int(yy[0])), min(int(yy[1]),black.shape[0]) )
        #     black[ sx, sy ] = 0

        # Set Dimensions
        pixel_width = config.getfloat('Data', 'channel_width') / GeoTransf['PixelSize'] # Channel width in Pixels
        radius = 2 * pixel_width # Circle Radius for Local Thresholding

        # # Mask Landsat NoData
        # print('nodata dilation...')
        # nanmask = np.where( mask==0, 1, 0 )
        # nanmask = mm.binary_dilation( nanmask, mm.disk( 30 ) )
        # mask = np.where( nanmask==1, 0, mask*black )

        # Image Cleaning
        print('cleaning mask...')
        mask = RemoveSmallObjects( mask, 100*pixel_width**2 ) # One Hundred Widths of Channel at Least is Required
        radius = max( np.floor( 0.5 * ( pixel_width ) ) - 3, 0 )
        mask = mm.binary_opening( mask, mm.disk( radius ) ) # Break Small Connectins
        mask = CleanIslands( mask, 10*pixel_width**2 ) # Clean Holes Inside the Planform
        mask = RemoveSmallObjects( mask, 100*pixel_width**2 ) # Remove New Small Objects
        mask = mask.astype( int )

        print('saving  mask and GeoTransf data...')
        np.save( maskfile, mask )
        with open( geofile, 'wb' ) as gf: pickle.dump( GeoTransf, gf )

    return None


def clean_masks( maskdir, geodir=None, config=None, file_only=False ):
    '''
    clean_masks( maskdir, geodir=None, config=None, file_only=False )
    =================================================================

    Manually remove branches from masks

    Arguments
    ---------
    maskdir           directory containing all the mask files
    geodir            directory where GeoTransf instances are stored (default None)
    config            PyRIS' RawConfigParser instance (default None)
    file_only         run on a single file only (default None)

    Returns
    -------
    None
    '''

    if not file_only:
        maskfiles = sorted( [ os.path.join(maskdir, f) for f in os.listdir(maskdir) ] )
        if geodir is not None: geofiles = sorted( [ os.path.join(geodir, f) for f in os.listdir(geodir) ] )
        else: gofiles = [ None for i in range( len(maskfiles) ) ]
    else:
        maskfiles = [ maskdir ]
        if geodir is not None:
            geofiles = [ os.path.join( geodir, os.path.splitext(os.path.split( maskdir )[-1] )[0] + '.p' )  ]
            geofiles = geofiles if os.path.isfile(geofiles[0]) else [ None ]
        else: geofiles = [ None ]
    for ifile, (maskfile,geofile) in enumerate( zip( maskfiles, geofiles ) ):
        print('cleaning file %s' % maskfile)
        # Look for the corresponding landsat image
        bg = None
        if config is not None:
            yearday = os.path.splitext( os.path.basename(maskfile) )[0]
            ldir = config.get('Data', 'input')
            landsats = [os.path.join(ldir,d) for d in os.listdir(ldir)]
            for l in landsats:
                lyearday = '_'.join( os.path.basename(l) )
                if yearday==lyearday:
                    bg = imread(os.path.join( l, os.path.basename(l).strip()+'_B1.TIF' ))
                    break
        GeoTransf = pickle.load( open(geofile, 'rb') ) if geofile is not None else None
        mask = np.load( maskfile )
        bw = np.where( mask>0, 1, 0 )
        bw = MaskClean( bw, bg )()
        if GeoTransf is not None and config is not None:
            pixel_width = int(config.get('Data', 'channel_width')) / GeoTransf[ 'PixelSize' ]
        else:
            pixel_width = 3
        bw = RemoveSmallObjects( bw, 100*pixel_width**2 ) # Remove New Small Objects

        ans = None
        while ans not in ['y', 'n']: ans = input('overwrite mask file?[y/n] ')
        if ans == 'y':
            print('saving mask file')
            np.save( maskfile, mask*bw )
        else:
            print('skipping')
    return None

def skeletonize_all( maskdir, skeldir, config ):
    '''
    skeletonize_all( maskdir, skeldir, config )
    ===========================================

    Apply skeletonization and distance transform on all the masks

    Arguments
    ---------
    maskdir           directory containing all the mask files
    skeldir           directory where skeleton and distance files will be stored
    config            PyRIS' RawConfigParser instance

    Returns
    -------
    None
    '''

    maskfiles = sorted( [ os.path.join(maskdir, f) for f in os.listdir(maskdir) ] )

    for ifile, maskfile in enumerate( maskfiles ):
        # input
        name = os.path.splitext( os.path.basename( maskfile ) )[0]
        # output
        skelfile = os.path.join( skeldir, '.'.join(( name, 'npy' )) )

        # skip the files which have already been processes
        if os.path.isfile( skelfile ):
            print('data found for file %s - skipping '  % ( skelfile ))
            continue
        print
        print('Processing file %s' % ( maskfile ))

        mask = np.load( maskfile ).astype( int )
        num_features = mask.max()

        # Skeletonization
        print('skeletonizing...')
        skel, dist = Skeletonize( np.where(mask>0,1,0).astype( int ) ) # Compute Axis and Distance
        labelled_skel = skel.astype(int) * mask.astype(int)

        # Apply Brute-Force cleaning
        skel = skel.astype( int )
        skel[ dist==1 ] = 0
        skelabs, nsl = ndimage.measurements.label( skel, structure=np.ones((3,3)) )
        for sl in range(1,nsl+1):
            if (skelabs==sl).sum() <= 500: skel[ skelabs==sl ] = 0

        # Pruning
        print('pruning n=%d labelled elements...' % num_features)
        pruned = np.zeros( mask.shape, dtype=int )
        if ( skelabs>0 ).sum() == 0:
            print('Something missing when pruning current label. Skipping...')
            continue
        for lab in range( 1, num_features+1 ):
            print('pruning label %d...' % lab)
            pruned += Pruning( labelled_skel==lab, int(config.get('Pruning', 'prune_iter')), smooth=False ) # Remove Spurs
        pruned *= dist.astype( int )

        # Check the number of junctions
        p = Pruner( np.where(pruned>0,1,0) )
        p.BuildStrides()
        Njunctions = 0
        for i in range( p.strides.shape[0] ):
            for j in range( p.strides.shape[1] ):
                if p.strides[i,j,1,1] == 1:
                    s = p.strides[i,j].sum()
                    if s == 4:
                        matrix = p.strides[i,j].copy()
                        matrix[1,1] = 0
                        pos = list(zip( *np.where( matrix > 0 ) ))
                        dists = np.zeros( len(pos) )
                        for ip in range( len(pos) ):
                            ipp = ip+1 if ip+1 < len(pos) else 0
                            dists[ip] = np.sqrt( (pos[ip][0]-pos[ipp][0])**2 + (pos[ip][1]-pos[ipp][1])**2 )
                        if np.any( dists<=1.001 ): continue
                        Njunctions += 1
        if Njunctions > 15:
            print('''
            'Warning!'
            'Expected at least %s recursion level.'
            'Consider increasing the pruning iteration or calling pyris with the --clean-mask flag'
            ''' % Njunctions)
        np.save( skelfile, pruned )
    return None

def vectorize_all( geodir, maskdir, skeldir, config, axisdir, use_geo=True ):
    '''
    vectorize_all( geodir, maskdir, skeldir, config, axisdir, use_geo=True )
    ===========================================

    Extract planform centerline coordinates and metrics from each skeleton and distance map

    Arguments
    ---------
    geodir            directory containing all the GeoTransf files
    maskdir           directory containing all the mask files
    skeldir           directory containing all the skeleton and distance files
    config            PyRIS' RawConfigParser instance
    axisdir           directory where the planform data files will be stored
    use_geo           use georeference provided by GeoTrensf files (default True)

    Returns
    -------
    None
    '''

    skelfiles = sorted( [ os.path.join(skeldir, f) for f in os.listdir(skeldir) ] )
    if use_geo: geofiles = sorted( [ os.path.join(geodir, f) for f in os.listdir(geodir) ] )

    for ifile, skelfile in enumerate( skelfiles ):
        # input
        name = os.path.splitext( os.path.basename( skelfile ) )[0]
        maskfile = os.path.join( maskdir, '.'.join(( name, 'npy' )) )
        geofile = os.path.join( geodir, '.'.join(( name, 'p' )) )
        # output
        axisfile = os.path.join( axisdir, '.'.join(( name, 'npy' )) )

        # skip the files which have already been processes
        if os.path.isfile( axisfile ):
            print('data found for file %s - skipping '  % ( axisfile ))
            continue
        print
        print('Processing file %s' % ( maskfile ))

        # Load mask, skeleton and GeoFile
        if use_geo: GeoTransf = pickle.load( open( geofile, 'rb' ) )
        mask = np.load( maskfile ).astype( int )
        skel = np.load( skelfile ).astype( int )
        num_features = mask.max()
        
        # Centerline Extraction
        print('extracting centerline of n=%d labelled elements...' % num_features)
        axis = Line2D()
        for lab in range( num_features, 0, -1 ):
            print('extracting label %d...' % lab)
            pdist = skel*(mask==lab)
            curr_axis = ReadAxisLine( pdist, flow_from=config.get('Data', 'flow_from'),
                                      method=config.get('Axis', 'reconstruction_method'),
                                      MAXITER=int(config.get('Axis', 'maxiter')) )
            axis.join( curr_axis )

        # Interpolation
        print('parametric cublic spline interpolation of the centerline...')
        step = int( max( 1, 0.5*int( axis.B.mean() ) ) ) # Discard points if too close
        # Npoints = axis.L / (0.25*axis.B.mean()) # Spacing = width/4
        Npoints = int(axis.L / (0.25*axis.B.mean())) # Spacing = width/4
        PCSs = 0.25*axis.x[::step].size # Degree of smoothness = n. of data points
        
        # Pixelled PCS
        xp_PCS, yp_PCS, d1xp_PCS, d1yp_PCS, d2xp_PCS, d2yp_PCS = InterpPCS( # Build Spline
            axis.x[::step], axis.y[::step], s=PCSs, N=Npoints
            )
        sp_PCS, thetap_PCS, Csp_PCS = CurvaturePCS( xp_PCS, yp_PCS, d1xp_PCS, d1yp_PCS, d2xp_PCS, d2yp_PCS,
                                                    method=2, return_diff=False )

        Bp_PCS = WidthPCS( axis.s/axis.s[-1]*sp_PCS[-1], axis.B, sp_PCS )
       
        # GeoReferenced PCS
        if use_geo:
            s_PCS, theta_PCS, Cs_PCS, B_PCS = \
                sp_PCS*GeoTransf['PixelSize'], thetap_PCS, Csp_PCS/GeoTransf['PixelSize'], Bp_PCS*GeoTransf['PixelSize']
            GR = GeoReference( GeoTransf )
            x_PCS, y_PCS = GR.RefCurve( xp_PCS, yp_PCS )
        else:
            x_PCS, y_PCS, s_PCS, theta_PCS, Cs_PCS, B_PCS = xp_PCS, yp_PCS, sp_PCS, thetap_PCS, Csp_PCS, Bp_PCS

        # Save Axis Properties
        print('saving main channel data...')
        save( axisfile, ( x_PCS, y_PCS, s_PCS, theta_PCS, Cs_PCS, B_PCS, xp_PCS, yp_PCS ) )
    return None
        

def migration_rates( axisfiles, migdir, columns=(0,1), show=False, pfreq=1 ):
    '''
    migration_rates( axisfiles, migdir, columns=(0,1), show=False, pfreq=10 ):
    ===========================================

    Compute migration vectors and bend separation for all the planform centerlines

    Arguments
    ---------
    axisfiles         temporally ordered list of files containing centerline data
    migdir            directory where migration rates and bend separation data will be stored
    columns           columns representing x and y planform coordinates in the axisfiles
    show              show all the migrations in a window at the end of the computation (default False)
    pfreq             smoothing coefficient (>=1) for the channel curvature for bend separation (default 1: no smoothing)

    Returns
    -------
    None
    '''
    
    migfiles = [ os.path.join( migdir, os.path.basename( axisfile ) ) for axisfile in axisfiles ]
    X, Y = [], []
    for axisfile in axisfiles:
        axis = load( axisfile )
        x, y = axis[ columns[0] ], axis[ columns[1] ]
        X.append( x ), Y.append( y )
    migrations = AxisMigration( X, Y )( pfreq=pfreq )
    for i, migfile in enumerate( migfiles ):
        [ dx, dy, dz, ICWTC, BI, B12, BUD ] = [ m[i] for m in migrations ]
        data = np.vstack( (dx, dy, dz, ICWTC, BI, B12, BUD) )
        save( migfile, data )

    colors = [ cm.jet(x) for x in np.linspace( 0, 1, len(axisfiles) ) ]
    lws = np.linspace( 1, 5, len(axisfiles) )

    if show:
        plt.figure()
        for i, f1 in enumerate( axisfiles ):
            a1 = load(f1)
            m1 = load( os.path.join(migdir, os.path.basename( f1 )) )
            name = '/'.join( (os.path.splitext(os.path.basename(f1))[0].split('_')[::-1]) )
            x, y, s = a1[columns[0]], a1[columns[1]], a1[2]
            dx, dy = m1[0], m1[1]
            b = m1[4]
            db = np.ediff1d(b, to_begin=0)
            idx = np.where(db>0)[0]
            plt.plot( x, y, c=colors[i], lw=lws[i], label=name )
            plt.plot( x[idx], y[idx], 'o', c=colors[i] )
            for j in range(idx.size): plt.text( x[idx[j]], y[idx[j]], str(int(b[idx[j]])) )
            for j in range(0,a1.shape[1],5): plt.arrow( x[j], y[j], dx[j], dy[j], fc='k', ec='k' )
            I = m1[6]==2
            plt.plot( [x[I], x[I]+dx[I]], [y[I], y[I]+dy[I]], 'k', lw=4 )
        plt.axis( 'equal' )
        plt.legend( loc='best' )
        plt.title( 'Migration rates' )
        plt.show()


def bars_detection( landsat_dirs, geodir, axisdir, migdir, bardir, show=False ):
    '''
    bars_detection( landsat_dirs, geodir, axisdir, migdir, bardir, show=False ):
    ===========================================

    Compute bare sediment bars morphodynamics

    Arguments
    ---------
    landsat_dirs      temporally ordered list of directories containing the Landsat data
    geodir            Directory containing GeoTransf files
    axisdir           directory containing planform data files
    migdir            directory containing migration vectors and bend separation files
    bardir            directory where bare sediment bar morphodynamic files will be stored
    show              show results in a window at the end of the computation (default False)

    Returns
    -------
    None
    '''
    
    axis_files = [ os.path.join( axisdir, f ) for f in sorted( os.listdir( axisdir ) ) ]
    found_files = []
    bars = FreeTemporalBars()
    
    for i_file, axis_file in enumerate( axis_files ):
        basename = os.path.splitext( os.path.split( axis_file )[-1] )[0]
        geo_file = os.path.join( geodir, '.'.join((basename,'p')) )
        mig_file = os.path.join( migdir, os.path.split( axis_file )[-1] )
        year, jday = [ v for v in basename.split('_') ]
        time = ( float(year) + float(jday)/365 )
        name = '%s%s' % ( year, jday )
        landsat_found = False
        for landsat_dir in landsat_dirs:
            lname = ''.join( get_year_jday( os.path.split(landsat_dir)[-1] ) )
            if name == lname:
                landsat_found = True
                break
        if not landsat_found:
            print('Landsat data not found for %s. Skipping...' % basename)
            continue
            found_files.append( axis_file )
        
        close=True
        remove_small=True

        print('Processing file %s' % basename)

        GeoTransf = pickle.load( open(geo_file, 'rb') )
        axis = load( axis_file )
        mig = load( mig_file )

         # ReLoad Landsat Data
        [ R, G, B, NIR, MIR, SWIR ], GeoTransf = LoadLandsatData( landsat_dir )
        bands = { 'R' : R, 'G' : G, 'B' : B, 'NIR' : NIR, 'MIR' : MIR, 'SWIR' : SWIR }

        # Compute Transformed Coordinates and Interpolate Band
        unwrapper = Unwrapper( axis, mig, GeoTransf )
        unwrapper.unwrap( R.shape, Npts=100 )
    
        # Find and Classify Bars
        barfinder = BarFinder( unwrapper )
        
        barfinder( bands, close=close, remove_small=remove_small )
        #barfinder.Show( bands )
        bars.GetFinder( time, barfinder )

    S, N, X, Y, Z, Xi, Yi, Si, Ni, DSi, DNi, YYi, Wi, BarCorr = bars.CentroidsEvol( 0 )

    np.save( os.path.join(bardir, 'interpolation.npy'), (S,N,X,Y,Z) )
    np.save( os.path.join(bardir, 'bendsep.npy'), (bars.Bars[0].unwrapper.Bend, bars.Bars[0].unwrapper.b) )
    np.save( os.path.join(bardir, 'values.npy'), (Xi, Yi, Si, Ni, DSi, DNi, YYi, Wi) )
    np.save( os.path.join(bardir, 'barcorr.npy'), BarCorr )

    if show: bars.Show( landsat_dirs, geodir )
    return None
