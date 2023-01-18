# =======================================
# Module: config
# Package: PyRIS
# Author: Federico Monegaglia
# Date: April 2016
# Description: Local Configuration Module
# =======================================

import configparser as ConfigParser

def default_config():
    '''Set up a configuration file for PyRIS'''
    
    
    cf = ConfigParser.RawConfigParser()
    
    cf.add_section( 'Data' )
    cf.set( 'Data', 'input', '' )
    cf.set( 'Data', 'output', '' )
    cf.set( 'Data', 'flow_from', '' )
    cf.set( 'Data', 'channel_width', None )

    cf.add_section( 'Segmentation' )
    cf.set( 'Segmentation', 'method', 'MIX' )
    cf.set( 'Segmentation', 'thresholding', 'global' )
    cf.set( 'Segmentation', 'black_masks', [] )

    cf.add_section( 'Pruning' )
    cf.set( 'Pruning', 'prune_iter', 25 )

    cf.add_section( 'Axis' )
    cf.set( 'Axis', 'reconstruction_method', 'std' )
    cf.set( 'Axis', 'maxiter', 100000 )

    return cf


def create_cfg_file( cf, fname='' ):
    '''Dump configuration cf on file fname'''
    with open( fname, 'w' ) as cfile: cf.write( cfile )
    return None


def get_cfg( fname ):
    '''Get Configuration'''
    cf = ConfigParser.RawConfigParser()
    cf.read( fname )
    return cf


def set_cfg( cf, section, name, value ):
    '''Set up a value'''
    cf.set( section, name, value )
    return None
