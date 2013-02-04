#include "sca.h"

/******************************************************************************
MODULE:  scene_based_sca

PURPOSE:  Calculate the snow cover mask, cloud cover mask, deep shadow mask,
and shaded relieve for the current scene, provided the TOA reflectance,
brightness temperature, and elevation.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the snow cover
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/31/2012  Gail Schmidt     Original Development

NOTES:
  1. The scene-based snow cover mask is based on an algorithm developed by
     Dave Selkowitz, Research Geographer, USGS Alaska Science Center.
  2. Processing will occur on a subset of lines at a time, however, the snow
     related buffers for snow cover and probability will need to be full-scene
     buffers to handle the 7x7 window post-processing.  The other option is
     to write out the NPROC_LINES buffer to a file, then read that file back
     for post-processing.  The later seems inefficient, and the former should
     be doable since the masks are 8-bit unsigned integers.
******************************************************************************/
int main (int argc, char *argv[])
{
    bool verbose;            /* verbose flag for printing messages */
    bool dem_top;            /* are we at the top of the dem for shaded
                                relief processing */
    bool dem_bottom;         /* are we at the bottom of the dem for shaded
                                relief processing */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char *toa_infile=NULL;   /* input TOA filename */
    char *btemp_infile=NULL; /* input brightness temperature filename */
    char *dem_infile=NULL;   /* input DEM filename */
    char *sc_outfile=NULL;   /* output snow cover filename */
    int retval;              /* return status */
    int k;                   /* variable to keep track of the % complete */
    int band;                /* current band to be processed */
    int line;                /* current line to be processed */
    int nlines_proc;         /* number of lines to process at one time */
    int start_line;          /* line of DEM to start reading */
    int extra_lines_start;   /* number of extra lines at the start of the DEM
                                to be read as part of the 3x3 window */
    int extra_lines_end;     /* number of extra lines at the end of the DEM
                                to be read as part of the 3x3 window */
    int nvals_to_read;       /* actual number of values to be read from the DEM
                                to include padding for the 3x3 */
    int offset;              /* offset in the raw binary DEM file to seek to
                                to begin reading the window in the DEM */
    int curr_snow_pix;       /* starting location/pixel of the current line
                                in the snow-cover related arrays, which are
                                full scene buffers */
    uint8 *refl_qa_mask=NULL;/* quality mask for the TOA reflectance products,
                                where fill and saturated values are flagged */
    uint8 *btemp_qa_mask=NULL;/* quality mask for the brightness temp products,
                                where fill and saturated values are flagged */
    uint8 *cloud_mask=NULL;  /* cloud mask */
    uint8 *snow_mask=NULL;   /* snow cover mask */
    uint8 *snow_prob=NULL;   /* snow cover probability score (percentage) */
    uint8 *deep_shad_mask=NULL; /* deep shadow mask */
    uint8 *shaded_relief=NULL;  /* shaded relief values */
    int16 *dem=NULL;         /* input DEM data (meters) */
    Input_t *toa_input=NULL; /* input structure for both the TOA reflectance
                                and brightness temperature products */
    FILE *dem_fptr=NULL;     /* input scene-based DEM file pointer */
    FILE *scm_fptr=NULL;     /* snow cover mask file pointer */
    FILE *sc_prob_fptr=NULL; /* snow cover probability file pointer */
    FILE *cm_fptr=NULL;      /* cloud mask file pointer */
    FILE *dsm_fptr=NULL;     /* deep shadow mask file pointer */
    FILE *relief_fptr=NULL;  /* shade relief file pointer */
    FILE *refl_qa_fptr=NULL; /* TOA reflectance QA */
    FILE *btemp_qa_fptr=NULL;/* brightness temp QA */

    printf ("Starting scene-based snow cover processing ...\n");

    /* Temporary -- open the mask output files for raw binary output */
    scm_fptr = fopen ("snow_cover_mask.bin", "wb");
    sc_prob_fptr = fopen ("snow_cover_probability_score.bin", "wb");
    cm_fptr = fopen ("cloud_mask.bin", "wb");
    dsm_fptr = fopen ("deep_shadow_mask.bin", "wb");
    relief_fptr = fopen ("shade_relief.bin", "wb");
    refl_qa_fptr = fopen ("toa_refl_qa.bin", "wb");
    btemp_qa_fptr = fopen ("btemp_qa.bin", "wb");

    /* Read the command-line arguments, including the name of the input
       Landsat TOA reflectance product and the DEM */
    retval = get_args (argc, argv, &toa_infile, &btemp_infile, &dem_infile,
        &sc_outfile, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Provide user information is verbose is turned on */
    if (verbose)
    {
        printf ("  TOA reflectance input file: %s\n", toa_infile);
        printf ("  Brightness temp input file: %s\n", btemp_infile);
        printf ("  DEM input file: %s\n", dem_infile);
        printf ("  Snow cover output file: %s\n", sc_outfile);
    }

    /* Open the TOA reflectance and brightness temperature products, set up
       the input data structure, allocate memory for the data buffers, and
       read the associated metadata and attributes. */
    toa_input = open_input (toa_infile, btemp_infile);
    if (toa_input == (Input_t *)NULL)
    {
        sprintf (errmsg, "Error opening/reading the TOA reflectance file: %s "
            "and the brightness temperature file: %s", toa_infile,
            btemp_infile);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Output some information from the input files if verbose */
    if (verbose)
    {
        printf ("  WRS path/row: %03d/%02d\n", toa_input->meta.path,
            toa_input->meta.row);
        printf ("  Number of lines/samples: %d/%d\n", toa_input->nlines,
            toa_input->nsamps);
        printf ("  Number of reflective bands: %d\n", toa_input->nrefl_band);
        printf ("  Number of thermal bands: %d\n", toa_input->nbtemp_band);
        printf ("  Pixel size: %f\n", toa_input->meta.pixsize);
        printf ("  Solar elevation angle: %f radians (%f degrees)\n",
            toa_input->meta.solar_elev, toa_input->meta.solar_elev * DEG);
        printf ("  Solar azimuth angle: %f radians (%f degrees)\n",
            toa_input->meta.solar_az, toa_input->meta.solar_az * DEG);
        printf ("  Fill value (refl, btemp): %d, %d\n", toa_input->refl_fill,
            toa_input->btemp_fill);
        printf ("  Scale factor (refl, btemp): %f, %f\n",
            toa_input->refl_scale_fact, toa_input->btemp_scale_fact);
        printf ("  Saturation value (refl, btemp): %d, %d\n",
            toa_input->refl_saturate_val, toa_input->btemp_saturate_val);
    }

    /* Allocate memory for the TOA reflectance QA mask */
    refl_qa_mask = (uint8 *) calloc (PROC_NLINES * toa_input->nsamps,
        sizeof (uint8));
    if (refl_qa_mask == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the TOA reflectance "
            "QA mask");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Allocate memory for the brightness temperature QA mask */
    btemp_qa_mask = (uint8 *) calloc (PROC_NLINES * toa_input->nsamps,
        sizeof (uint8));
    if (btemp_qa_mask == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the brightness temp "
            "QA mask");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Allocate memory for the cloud mask */
    cloud_mask = (uint8 *) calloc (PROC_NLINES * toa_input->nsamps,
        sizeof (uint8));
    if (cloud_mask == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the cloud mask");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Allocate memory for the snow cover mask (full scene to handle
       post-processing) */
    snow_mask = (uint8 *) calloc (toa_input->nlines * toa_input->nsamps,
        sizeof (uint8));
    if (snow_mask == NULL)
    {
        sprintf (errmsg, "Error allocating memory (full scene) for the snow "
            "mask");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Allocate memory for the snow cover probability (full scene to handle
       post-processing) */
    snow_prob = (uint8 *) calloc (toa_input->nlines * toa_input->nsamps,
        sizeof (uint8));
    if (snow_prob == NULL)
    {
        sprintf (errmsg, "Error allocating memory (full scene) for the snow "
            "cover probability");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
    {
        printf ("  Processing %d lines at a time\n", PROC_NLINES);
        printf ("  Cloud and snow cover -- %% complete: 0%%\r");
    }

    /* Loop through the lines and samples in the TOA reflectance and
       brightness temperature products, computing the cloud and snow cover */
    nlines_proc = PROC_NLINES;
    k = 0;
    for (line = 0; line < toa_input->nlines; line += PROC_NLINES)
    {
        /* Do we have nlines_proc left to process? */
        if (line + nlines_proc >= toa_input->nlines)
            nlines_proc = toa_input->nlines - line;

        /* Update processing status? */
        if (verbose && (100 * line / toa_input->nlines > k))
        {
            k = 100 * line / toa_input->nlines;
            if (k % 10 == 0)
            {
                printf ("  Cloud and snow cover -- %% complete: %d%%\r", k);
                fflush (stdout);
            }
        }

        /* Read the current lines from the TOA file for each of the TOA
           reflectance bands */
        for (band = 0; band < toa_input->nrefl_band; band++)
        {
            if (get_input_refl_lines (toa_input, band, line, nlines_proc) !=
                SUCCESS)
            {
                sprintf (errmsg, "Error reading %d lines from band %d of the "
                    "TOA reflectance file starting at line %d", nlines_proc,
                    band, line);
                error_handler (true, FUNC_NAME, errmsg);
                close_input (toa_input);
                free_input (toa_input);
                exit (ERROR);
            }
        }  /* end for band */

        /* Read the current lines from the brightness temp file */
        if (get_input_btemp_lines (toa_input, line, nlines_proc) != SUCCESS)
        {
            sprintf (errmsg, "Error reading %d lines from the brightness "
                "temperature file starting at line %d", nlines_proc, line);
            error_handler (true, FUNC_NAME, errmsg);
            close_input (toa_input);
            free_input (toa_input);
            exit (ERROR);
        }

        /* Set up mask for the TOA reflectance values */
        refl_mask (toa_input->refl_buf[0] /*b1*/, toa_input->refl_buf[1] /*b2*/,
            toa_input->refl_buf[2] /*b3*/, toa_input->refl_buf[3] /*b4*/,
            toa_input->refl_buf[4] /*b5*/, toa_input->refl_buf[5] /*b7*/,
            nlines_proc, toa_input->nsamps, toa_input->refl_fill,
            refl_qa_mask);

        /* Set up mask for the brightness temperature values */
        btemp_mask (toa_input->btemp_buf, nlines_proc, toa_input->nsamps,
            toa_input->btemp_fill, btemp_qa_mask);

        /* Compute the cloud mask */
        cloud_cover_class (toa_input->refl_buf[0] /*b1*/,
            toa_input->refl_buf[3] /*b4*/, toa_input->btemp_buf /*b6*/,
            toa_input->refl_buf[5] /*b7*/, nlines_proc, toa_input->nsamps,
            toa_input->refl_scale_fact, toa_input->btemp_scale_fact,
            refl_qa_mask, btemp_qa_mask, cloud_mask);

        /* Find the location of the current line in the snow cover and snow
           probability masks, since they are full scene array buffers */
        curr_snow_pix = line * toa_input->nsamps;

        /* Compute the snow cover mask */
        snow_cover_class (toa_input->refl_buf[0] /*b1*/,
            toa_input->refl_buf[1] /*b2*/, toa_input->refl_buf[2] /*b3*/,
            toa_input->refl_buf[3] /*b4*/, toa_input->refl_buf[4] /*b5*/,
            toa_input->btemp_buf /*b6*/, toa_input->refl_buf[5] /*b7*/,
            nlines_proc, toa_input->nsamps, toa_input->refl_scale_fact,
            toa_input->btemp_scale_fact, toa_input->refl_saturate_val,
            refl_qa_mask, &snow_mask[curr_snow_pix], &snow_prob[curr_snow_pix]);

        /* Temporary - write the non snow-related masks to raw binary output */
        fwrite (refl_qa_mask, 1, nlines_proc*toa_input->nsamps * sizeof(uint8),
            refl_qa_fptr);
        fwrite (btemp_qa_mask, 1, nlines_proc*toa_input->nsamps * sizeof(uint8),
            btemp_qa_fptr);
        fwrite (cloud_mask, 1, nlines_proc*toa_input->nsamps * sizeof(uint8),
            cm_fptr);
    }  /* end for line */

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Cloud and snow cover -- %% complete: 100%%\n");

    /* Full scene post-processing of the snow cover pixels to deal with
       false positives in the 90% confidence branch */
    if (verbose)
        printf ("  Post-processing snow cover mask.\n");
    post_process_snow_cover_class (toa_input->nlines, toa_input->nsamps,
        snow_mask, snow_prob);

    /* Temporary - write the snow-related masks to raw binary output */
    fwrite (snow_mask, 1, toa_input->nlines * toa_input->nsamps * sizeof(uint8),
        scm_fptr);
    fwrite (snow_prob, 1, toa_input->nlines * toa_input->nsamps * sizeof(uint8),
        sc_prob_fptr);

    /* Temporary -- close the mask output files */
    fclose (scm_fptr);
    fclose (sc_prob_fptr);
    fclose (cm_fptr);

    /* Free the mask pointers */
    if (refl_qa_mask != NULL)
    {
        free (refl_qa_mask);
        refl_qa_mask = NULL;
    }
    if (btemp_qa_mask != NULL)
    {
        free (btemp_qa_mask);
        btemp_qa_mask = NULL;
    }
    if (cloud_mask != NULL)
    {
        free (cloud_mask);
        cloud_mask = NULL;
    }
    if (snow_mask != NULL)
    {
        free (snow_mask);
        snow_mask = NULL;
    }
    if (snow_prob != NULL)
    {
        free (snow_prob);
        snow_prob = NULL;
    }

    /* Open the DEM for reading raw binary */
    dem_fptr = fopen (dem_infile, "rb");
    if (dem_fptr == NULL)
    {
        sprintf (errmsg, "Error opening the DEM file: %s", dem_infile);
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Allocate memory for the DEM, which will hold PROC_NLINES of data.  The
       DEM should be the same size as the input scene, since the scene was
       used to resample the DEM.  To process the shaded relieve we need to
       read an extra two lines to process a 3x3 window. */
    dem = (int16 *) calloc ((PROC_NLINES+2) * toa_input->nsamps, sizeof(int16));
    if (dem == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the DEM data");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Allocate memory for the shaded relief and deep shadow mask */
    deep_shad_mask = (uint8 *) calloc (PROC_NLINES * toa_input->nsamps,
        sizeof (uint8));
    if (deep_shad_mask == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the deep shadow mask");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    shaded_relief = (uint8 *) calloc (PROC_NLINES * toa_input->nsamps,
        sizeof (uint8));
    if (shaded_relief == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the shaded relief");
        error_handler (true, FUNC_NAME, errmsg);
        close_input (toa_input);
        free_input (toa_input);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
    {
        printf ("  Processing %d lines at a time\n", PROC_NLINES);
        printf ("  Shaded relief -- %% complete: 0%%\r");
    }

    /* Loop through the lines and samples in the DEM and compute the
       relief shading */
    nlines_proc = PROC_NLINES;
    k = 0;
    for (line = 0; line < toa_input->nlines; line += nlines_proc)
    {
        /* Do we have nlines_proc left to process? */
        if (line + nlines_proc >= toa_input->nlines)
            nlines_proc = toa_input->nlines - line;

        /* Update processing status? */
        if (verbose && (100 * line / toa_input->nlines > k))
        {
            k = 100 * line / toa_input->nlines;
            if (k % 10 == 0)
            {
                printf ("  Shaded relief -- %% complete: %d%%\r", k);
                fflush (stdout);
            }
        }

        /* Prepare to read the current lines from the DEM.  We need an extra
           line at the start and end for the shaded relief.  If we are just
           starting at line 0, then read an extra line from the end of the
           image window.  If we aren't at line 0, then read an extra line at
           the start and end of the image window.  If we are at the end of the
           image, then read an extra line from the start of the image window. */
        start_line = line - 1;
        extra_lines_start = 1;
        dem_top = false;
        dem_bottom = false;
        if (start_line < 0)
        {   /* adjust for the number of lines at the start which we
               don't have then read at line 0 */
            extra_lines_start = 0;
            start_line = 0;
            dem_top = true;
        }
        if (line + nlines_proc < toa_input->nlines)
            extra_lines_end = 1;
        else
        {
            extra_lines_end = 0;
            dem_bottom = true;
        }

        /* Start reading from the start_line and read
           nlines_proc + extra_lines_start + extra_lines_end lines */
        nvals_to_read = (nlines_proc+extra_lines_start+extra_lines_end) *
            toa_input->nsamps;
        offset = sizeof (int16) * start_line * toa_input->nsamps;
        fseek (dem_fptr, offset, SEEK_SET);
        if (fread (dem, sizeof (int16), nvals_to_read, dem_fptr)
            != nvals_to_read)
        {
            sprintf (errmsg, "Error reading %d values from the DEM file "
                "starting at line %d.", nvals_to_read, start_line);
            error_handler (true, FUNC_NAME, errmsg);
            close_input (toa_input);
            free_input (toa_input);
            fclose (dem_fptr);
            exit (ERROR);
        }

        /* Initialize the shaded relief and deep shadow masks to 0s for the
           current window.  The first and last pixel will not get processed
           in the deep_shadow mask.  The first and last line in the entire
           image also will not get processed in the deep_shadow mask. */
        memset ((void *) deep_shad_mask, 0, PROC_NLINES * toa_input->nsamps
            * sizeof (uint8));
        memset ((void *) shaded_relief, 0, PROC_NLINES * toa_input->nsamps
            * sizeof (uint8));

        /* Compute the shaded relief and associated terrain-derived deep
           shadow mask */
        deep_shadow (dem, dem_top, dem_bottom, nlines_proc, toa_input->nsamps,
            toa_input->meta.pixsize, toa_input->meta.pixsize,
            toa_input->meta.solar_elev, toa_input->meta.solar_az,
            shaded_relief, deep_shad_mask);

        /* Temporary - write the shaded relief and deep shadow mask to raw
           binary output */
        fwrite (deep_shad_mask, 1, nlines_proc*toa_input->nsamps *
            sizeof(uint8), dsm_fptr);
        fwrite (shaded_relief, 1, nlines_proc*toa_input->nsamps * sizeof(uint8),
            relief_fptr);
    }  /* end for line */

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Shaded relief -- %% complete: 100%%\n");

    /* Temporary -- close the mask output files for raw binary output */
    fclose (dem_fptr);
    fclose (dsm_fptr);
    fclose (relief_fptr);

    /* Free the TOA reflectance and brightness temperature pointers and
       close the TOA reflectance and brightness temperature products */
    close_input (toa_input);
    free_input (toa_input);
    toa_input = NULL;

    /* Free the filename pointers */
    if (toa_infile != NULL)
        free (toa_infile);
    if (btemp_infile != NULL)
        free (btemp_infile);
    if (dem_infile != NULL)
        free (dem_infile);
    if (sc_outfile != NULL)
        free (sc_outfile);

    /* Free the mask pointers */
    if (dem != NULL)
        free (dem);
    if (deep_shad_mask != NULL)
        free (deep_shad_mask);
    if (shaded_relief != NULL)
        free (shaded_relief);

    /* Indicate successful completion of processing */
    printf ("Scene-based snow cover processing complete!\n");
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/2/2013    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("scene_based_snow_cover handles the snow cover classification for "
            "the input Landsat scene (top of atmosphere reflection), "
            "including a cloud cover assessment and determining the relief "
            "shading based on the terrain.\n\n");
    printf ("usage: scene_based_snow_cover "
            "--toa=input_TOA_reflectance_Landsat_filename "
            "--btemp=input_brightness_temperature_Landsat_filename "
            "--dem=input_DEM_filename "
            "--snow_cover=output_snow_cover_filename "
            "[--verbose]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -toa: name of the input Landsat TOA reflectance file to be "
            "classified (HDF)\n");
    printf ("    -btemp: name of the input Landsat brightness temperature file "
            "to be classified (HDF)\n");
    printf ("    -dem: name of the DEM associated with the Landsat TOA file "
            "(raw binary 16-bit integers)\n");
    printf ("    -snow_cover: name of the output snow cover file (HDF)\n");
    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");
    printf ("\nscene_based_snow_cover --help will print the usage statement\n");
    printf ("\nExample: snow_cover "
            "--toa=lndcal.LT50400331995173AAA02.hdf "
            "--btemp=lndth.LT50400331995173AAA02.hdf "
            "--dem=lsrd_scene_based_dem.bin "
            "--snow_cover=snow_cover.hdf "
            "--verbose\n");
}
