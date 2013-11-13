#! /usr/bin/env python

'''
Created on January 29, 2013 by Gail Schmidt, USGS/EROS

License:
  "NASA Open Source Agreement 1.3"

Description:
  Run the scene-based snow cover algorithm, including checking processing
  status.

History:
  Updated on ??? by Gail Schmidt, USGS/EROS
    - {specify modifications}
  Updated on October/2013 by Ron Dilley, USGS/EROS
    - Modified to use espa-common.
    - Modified to require the DEM file to use on the command line.

Usage:
  do_snow_cover.py --help prints the help message
'''

import sys
import os
import re
import commands
import datetime
from optparse import OptionParser

from espa_constants import *
from espa_logging import log

class Sca():
    '''
    Defines the class object for performing snow covered area processing
    '''

    def __init__(self):
        pass


    ########################################################################
    # Description: runSca will use the parameters passed for the input/output
    # files, logfile, and usebin.  If input/output files are None (i.e. not
    # specified) then the command-line parameters will be parsed for this
    # information.  The scene-based DEM and snow cover applications are then
    # executed on the specified input files.  If a log file was specified,
    # then the output from each application will be logged to that file.
    #
    # Inputs:
    #   metafile - name of the Landsat metadata file to be processed
    #   toa_infile - name of the input TOA reflectance HDF file to be processed
    #   dem - name of an ENVI formatted DEM file
    #   btemp_infile - name of the input brightness temp HDF file to
    #       be processed
    #   sca_outfile - name of the output snow cover HDF file
    #   logfile - name of the logfile for logging information; if None then
    #       the output will be written to stdout
    #   usebin - this specifies if the DEM and SCE exes reside in the $BIN
    #       directory; if None then the DEM and SCA exes are expected to be in
    #       the PATH
    #
    # Returns:
    #     ERROR - error running the DEM and SCA applications
    #     SUCCESS - successful processing
    #
    # Notes:
    #   1. The script obtains the path of the metadata file and changes
    #      directory to that path for running the scene-based DEM and SCA
    #      code.  If the metafile directory is not writable, then this script
    #      exits with an error.
    #   2. If the metadata file is not specified and the information is going
    #      to be grabbed from the command line, then it's assumed all the
    #      parameters will be pulled from the command line.
    #   3. Currently the raw binary products are automatically written in
    #      addition to the HDF-EOS product for the snow cover algorithm.
    #      The raw binary products will be used for downstream processing,
    #      but only the output HDF-EOS product should be delivered to the
    #      general public.
    #######################################################################
    def runSca (self, metafile=None, toa_infile=None, dem=None, \
        btemp_infile=None, sca_outfile=None, logfile=None, usebin=None):
        # if no parameters were passed then get the info from the
        # command line
        if metafile == None:
            # get the command line argument for the metadata file
            parser = OptionParser()
            parser.add_option ("-f", "--metafile", type="string",
                dest="metafile",
                help="name of Landsat MTL file", metavar="FILE")
            parser.add_option ("-t", "--toa_infile", type="string",
                dest="toa_infile",
                help="name of TOA reflectance HDF file", metavar="FILE")
            parser.add_option ('--dem', type='string',
                dest='dem',
                help="name of an ENVI formatted DEM file", metavar='FILE')
            parser.add_option ("-b", "--btemp_infile", type="string",
                dest="btemp_infile",
                help="name of brightness temperature HDF file", metavar="FILE")
            parser.add_option ("-s", "--sca_outfile", type="string",
                dest="sca_outfile",
                help="name of output snow cover HDF file", metavar="FILE")
            parser.add_option ("--usebin", dest="usebin", default=False,
                action="store_true",
                help="use BIN environment variable as the location of DEM and SCA apps")
            parser.add_option ("-l", "--logfile", type="string", dest="logfile",
                help="name of optional log file", metavar="FILE")
            (options, args) = parser.parse_args()
    
            # validate the command-line options
            usebin = options.usebin          # should $BIN directory be used
            logfile = options.logfile        # name of the log file

            # metadata file
            metafile = options.metafile
            if metafile == None:
                parser.error ("missing metafile command-line argument");
                return ERROR

            # TOA reflectance file
            toa_infile = options.toa_infile
            if toa_infile == None:
                parser.error ("missing TOA input file command-line argument");
                return ERROR

            # ENVI DEM file
            if dem == None:
                dem = options.dem
                if dem == None:
                    parser.error ("missing input dem file command-line" \
                        " argument");
                    return ERROR

            # brightness temp file
            btemp_infile = options.btemp_infile
            if btemp_infile == None:
                parser.error ("missing brightness temperature input file command-line argument");
                return ERROR
        
            # snow cover output file
            sca_outfile = options.sca_outfile
            if sca_outfile == None:
                parser.error ("missing snow cover output file command-line argument");
                return ERROR
        
        # let the world know what we are processing
        msg = 'SCA processing of Landsat metadata file: %s' % metafile
        log (msg)
        
        # should we expect the DEM and SCA applications to be in the PATH or
        # in the BIN directory?
        if usebin:
            # get the BIN dir environment variable
            bin_dir = os.environ.get('BIN')
            bin_dir = bin_dir + '/'
            msg = 'BIN environment variable: %s' % bin_dir
            log (msg)
        else:
            # don't use a path to the DEM/SCA applications
            bin_dir = ""
            msg = 'DEM and SCA executables expected to be in the PATH'
            log (msg)
        
        # make sure the metadata file exists
        if not os.path.isfile(metafile):
            msg = "Error: metadata file does not exist or is not accessible: " + metafile
            log (msg)
            return ERROR

        # use the base metadata filename and not the full path.
        base_metafile = os.path.basename (metafile)
        msg = 'Processing metadata file: %s' % base_metafile
        log (msg)
        
        # get the path of the MTL file and change directory to that location
        # for running this script.  save the current working directory for
        # return to upon error or when processing is complete.  Note: use
        # abspath to handle the case when the filepath is just the filename
        # and doesn't really include a file path (i.e. the current working
        # directory).
        mydir = os.getcwd()
        metadir = os.path.dirname (os.path.abspath (metafile))
        if not os.access(metadir, os.W_OK):
            msg = 'Path of metadata file is not writable: %s.  Script needs write access to the metadata directory.' % metadir
            log (msg)
            return ERROR
        msg = 'Changing directories for snow cover processing: %s' % metadir
        log (msg)
        os.chdir (metadir)

        # run snow cover algorithm, checking the return status of each module.
        # exit if any errors occur.
        cmdstr = "%sscene_based_sca --toa=%s --btemp=%s --dem=%s --snow_cover=%s --write_binary --verbose" % (bin_dir, toa_infile, btemp_infile, dem, sca_outfile)
#        print 'DEBUG: scene_based_sca command: %s' % cmdstr
        (status, output) = commands.getstatusoutput (cmdstr)
        log (output)
        exit_code = status >> 8
        if exit_code != 0:
            msg = 'Error running scene_based_sca.  Processing will terminate.'
            log (msg)
            os.chdir (mydir)
            return ERROR
        
        # successful completion.  return to the original directory.
        os.chdir (mydir)
        msg = 'Completion of scene based snow cover.'
        log (msg)

        return SUCCESS

######end of Sca class######

if __name__ == "__main__":
    return_value = Sca().runSca()

    if (return_value != SUCCESS):
        sys.exit (EXIT_FAILURE)
    else:
        sys.exit (EXIT_SUCCESS)

