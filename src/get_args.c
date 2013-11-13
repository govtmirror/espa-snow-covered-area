#include <getopt.h>
#include "sca.h"

/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/2/2013    Gail Schmidt     Original Development
2/15/2013   Gail Schmidt     Added support for write raw binary flag

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
******************************************************************************/
short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **toa_infile,    /* O: address of input TOA filename */
    char **btemp_infile,  /* O: address of input TOA filename */
    char **dem_infile,    /* O: address of input DEM filename */
    char **sc_outfile,    /* O: address of output snow cover filename */
    bool *write_binary,   /* O: write raw binary flag */
    bool *verbose         /* O: verbose flag */
)
{
    int c;                           /* current argument index */
    int option_index;                /* index for the command-line option */
    static int verbose_flag=0;       /* verbose flag */
    static int binary_flag=0;        /* write binary flag */
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "get_args";   /* function name */
    static struct option long_options[] =
    {
        {"verbose", no_argument, &verbose_flag, 1},
        {"write_binary", no_argument, &binary_flag, 1},
        {"toa", required_argument, 0, 't'},
        {"btemp", required_argument, 0, 'b'},
        {"dem", required_argument, 0, 'd'},
        {"snow_cover", required_argument, 0, 's'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Loop through all the cmd-line options */
    opterr = 0;   /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {   /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
     
            case 'h':  /* help */
                usage ();
                return (ERROR);
                break;

            case 't':  /* toa infile */
                *toa_infile = strdup (optarg);
                break;
     
            case 'b':  /* btemp infile */
                *btemp_infile = strdup (optarg);
                break;
     
            case 'd':  /* dem infile */
                *dem_infile = strdup (optarg);
                break;
     
            case 's':  /* snow cover outfile */
                *sc_outfile = strdup (optarg);
                break;
     
            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind-1]);
                error_handler (true, FUNC_NAME, errmsg);
                usage ();
                return (ERROR);
                break;
        }
    }

    /* Make sure the infiles and outfiles were specified */
    if (*toa_infile == NULL)
    {
        sprintf (errmsg, "TOA input file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    if (*btemp_infile == NULL)
    {
        sprintf (errmsg, "Brightness temperature input file is a required "
            "argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    if (*dem_infile == NULL)
    {
        sprintf (errmsg, "DEM input file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    if (*sc_outfile == NULL)
    {
        sprintf (errmsg, "Snow cover output file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Check the write binary flag */
    if (binary_flag)
        *write_binary = true;

    /* Check the verbose flag */
    if (verbose_flag)
        *verbose = true;

    return (SUCCESS);
}
