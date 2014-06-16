#include <cv.h>
#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  revised_cloud_mask

PURPOSE:  Revised the cfmask cloud mask which is included as part of the ESPA
surface reflectance product.  The cfmask flags many snow pixels as cloud pixels
and this revision will clean up the cloud mask to correctly identify the snow
pixels.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the cloud mask
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
5/19/2014     Gail Schmidt     Original Development, based on an algorithm
                               provided by David Selkowitz

NOTES:
  1. The rule-based models expect that the variances for the TOA reflectance
     bands were calculated using scaled values straight from the TOA bands.
     Thus the scale-factor was not applied to unscale these values.  The
     models also expect the variances for the NDVI and NDSI values to be 
     computed on the original values between -1.0 and 1.0.  Thus these values
     need to be unscaled after being written to the output file as scaled.
******************************************************************************/
int main (int argc, char *argv[])
{
    bool verbose;              /* verbose flag for printing messages */
    bool toa_refl=true;        /* process TOA reflectance by default, but leave
                                  it open to use surface reflectance in the
                                  future */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];     /* error message */
    char envi_file[STR_SIZE];  /* name of the output ENVI header file */
    char short_cm_names[MAX_OUT_BANDS][STR_SIZE]; /* output short names for new
                                                     cloud mask bands */
    char long_cm_names[MAX_OUT_BANDS][STR_SIZE];  /* output long names for new
                                                     cloud mask bands */
    char cm_data_units[MAX_OUT_BANDS][STR_SIZE];  /* output data units for new
                                                     cloud mask bands */
    char *xml_infile=NULL;     /* input XML filename */
    char *cptr=NULL;           /* pointer to the file extension */
    int retval;                /* return status */
    int i;                     /* looping variable */
    int ib;                    /* looping variable for bands */
    int line, samp;            /* current line,samp to be processed */
    int nlines_proc;           /* number of lines to process at one time */
    int num_cm;                /* number of cloud mask products to be output */
    int16 *ndvi=NULL;          /* NDVI values */
    int16 *ndsi=NULL;          /* NDSI values */
    int16 *band_arr=NULL;      /* values for current reflectance band */
    int16 *index_arr=NULL;     /* values for current index */
    int32 *variance_arr=NULL;  /* variance values for current band/index */
    int32 *b1_var=NULL;        /* band1 variance values */
    int32 *b2_var=NULL;        /* band2 variance values */
    int32 *b3_var=NULL;        /* band3 variance values */
    int32 *b4_var=NULL;        /* band4 variance values */
    int32 *b5_var=NULL;        /* band5 variance values */
    int32 *b7_var=NULL;        /* band7 variance values */
    int32 *ndvi_var=NULL;      /* NDVI variance values */
    int32 *ndsi_var=NULL;      /* NDSI variance values */
    uint8 *rev_cm=NULL;        /* revised cloud mask */
    uint8 *rev_lim_cm=NULL;    /* revised cloud mask without variances */
    uint8 *buff_cm=NULL;       /* revised cloud mask with buffering */
    uint8 *opencv_img=NULL;    /* pointer to the opencv image data */
    void *ptr=NULL;            /* generic pointer for reading data */
    Input_t *refl_input=NULL;  /* input structure for the TOA product */
    Output_t *cm_output=NULL;  /* output structure and metadata for the new
                                  cloud mask products */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Envi_header_t envi_hdr;   /* output ENVI header information */
    IplImage* cv_img=NULL;    /* image pointer to hold the cloud mask imagery
                                 for OpenCV */
    IplConvKernel* element=NULL; /* 5x5 kernel for the erosion and dilation */

    printf ("Starting revised cloud mask processing ...\n");

    /* Read the command-line arguments */
    retval = get_args (argc, argv, &xml_infile, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Provide user information if verbose is turned on */
    if (verbose)
        printf ("  XML input file: %s\n", xml_infile);

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile, ESPA_SCHEMA) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_infile, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Open the reflectance product, set up the input data structure, and
       allocate memory for the data buffers */
    refl_input = open_input (&xml_metadata, toa_refl);
    if (refl_input == (Input_t *) NULL)
    {
        sprintf (errmsg, "Error opening/reading the reflectance data: %s",
            xml_infile);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Output some information from the input files if verbose */
    if (verbose)
    {
        printf ("  Number of lines/samples: %d/%d\n", refl_input->nlines,
            refl_input->nsamps);
        printf ("  Number of reflective bands: %d\n", refl_input->nrefl_band);
        printf ("  Fill value: %d\n", refl_input->refl_fill);
        printf ("  Scale factor: %f\n", refl_input->refl_scale_fact);
        printf ("  Saturation value: %d\n", refl_input->refl_saturate_val);
    }

    /* Allocate memory for the NDVI and NDSI, holds PROC_NLINES */
    ndvi = calloc (PROC_NLINES*refl_input->nsamps, sizeof (int16));
    if (ndvi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the NDVI");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    ndsi = calloc (PROC_NLINES*refl_input->nsamps, sizeof (int16));
    if (ndsi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the NDSI");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Set up the output information for the NDVI and NDSI */
    num_cm = NUM_CM;
    strcpy (short_cm_names[CM_NDVI], "ndvi");
    strcpy (long_cm_names[CM_NDVI], "normalized difference vegetation index");
    strcpy (cm_data_units[CM_NDVI], "band ratio index value");

    strcpy (short_cm_names[CM_NDSI], "ndsi");
    strcpy (long_cm_names[CM_NDSI], "normalized difference snow index");
    strcpy (cm_data_units[CM_NDSI], "band ratio index value");

    for (i = VARIANCE_B1; i <= VARIANCE_B5; i++)
    {
        sprintf (short_cm_names[i], "varb%d", i-VARIANCE_B1+1);
        sprintf (long_cm_names[i], "variance for band %d", i-VARIANCE_B1+1);
        strcpy (cm_data_units[i], "squared deviation from mean");
    }

    strcpy (short_cm_names[VARIANCE_B7], "varb7");
    strcpy (long_cm_names[VARIANCE_B7], "variance for band 7");
    strcpy (cm_data_units[VARIANCE_B7], "squared deviation from mean");

    strcpy (short_cm_names[VARIANCE_NDVI], "varndvi");
    strcpy (long_cm_names[VARIANCE_NDVI], "variance for NDVI");
    strcpy (cm_data_units[VARIANCE_NDVI], "squared deviation from mean");

    strcpy (short_cm_names[VARIANCE_NDSI], "varndsi");
    strcpy (long_cm_names[VARIANCE_NDSI], "variance for NDSI");
    strcpy (cm_data_units[VARIANCE_NDSI], "squared deviation from mean");

    strcpy (short_cm_names[REVISED_CM], "revcm");
    strcpy (long_cm_names[REVISED_CM], "revised cloud mask");
    strcpy (cm_data_units[REVISED_CM], "quality/feature classification");

    strcpy (short_cm_names[REVISED_LIM_CM], "revlimcm");
    strcpy (long_cm_names[REVISED_LIM_CM], "revised limited cloud mask");
    strcpy (cm_data_units[REVISED_LIM_CM], "quality/feature classification");

    /* Open the specified output files and create the metadata structure */
    cm_output = open_output (&xml_metadata, refl_input, num_cm,
        short_cm_names, long_cm_names, cm_data_units, toa_refl);
    if (cm_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
    {
        printf ("  Processing spectral indices %d lines at a time\n",
            PROC_NLINES);
    }

    /* Loop through the lines and samples in the reflectance product,
       computing the NDVI and NDSI */
    nlines_proc = PROC_NLINES;
    for (line = 0; line < refl_input->nlines; line += PROC_NLINES)
    {
        /* Do we have nlines_proc left to process? */
        if (line + nlines_proc >= refl_input->nlines)
            nlines_proc = refl_input->nlines - line;

        /* Read the current lines from the reflectance file for each of the
           reflectance bands */
        for (ib = 0; ib < refl_input->nrefl_band; ib++)
        {
            if (get_input_refl_lines (refl_input, ib, line, nlines_proc, NULL)
                != SUCCESS)
            {
                sprintf (errmsg, "Error reading %d lines from band %d of the "
                    "reflectance file starting at line %d", nlines_proc, ib,
                    line);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }  /* end for ib */

        /* Compute the NDVI and write to output file
           NDVI = (nir - red) / (nir + red) */
        make_index (refl_input->refl_buf[3] /*b4*/,
            refl_input->refl_buf[2] /*b3*/, refl_input->refl_fill,
            refl_input->refl_saturate_val, nlines_proc, refl_input->nsamps,
            ndvi);

        if (put_output_lines (cm_output, ndvi, CM_NDVI, line, nlines_proc,
            sizeof (int16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output NDVI data for line %d", line);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Compute the NDSI and write to output file
           NDSI = (green - mir) / (green + mir) */
        make_index (refl_input->refl_buf[1] /*b2*/,
            refl_input->refl_buf[4] /*b5*/, refl_input->refl_fill,
            refl_input->refl_saturate_val, nlines_proc, refl_input->nsamps,
            ndsi);

        if (put_output_lines (cm_output, ndsi, CM_NDSI, line, nlines_proc,
            sizeof (int16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output NDSI data for line %d", line);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end for line */

    /* Free the index pointers */
    free (ndvi);
    free (ndsi);

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Spectral indices -- complete\n");

    /* Allocate memory for the variance array, whole band */
    variance_arr = calloc (refl_input->nlines * refl_input->nsamps,
        sizeof (int32));
    if (variance_arr == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the variance array");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for the band array to read the entire reflectance band
       vs just a subset of lines, whole band */
    band_arr = calloc (refl_input->nlines * refl_input->nsamps,
        sizeof (int16));
    if (band_arr == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band array");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Processing variance for reflectance bands and indices\n");

    /* Compute the variance for each reflectance band, NDVI, and NDSI.  Write
       these to the output.  Process the whole band or index at a time, one
       at a time. */ 
    /* reflectance bands */
    for (ib = 0; ib < refl_input->nrefl_band; ib++)
    {
        /* Read the band */
        if (verbose)
            printf ("    - reflectance band %d\n", ib);
        if (get_input_refl_lines (refl_input, ib, 0, refl_input->nlines,
            band_arr) != SUCCESS)
        {
            sprintf (errmsg, "Error reading band %d of reflectance file", ib);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Compute the variance.  Don't unscale the reflectance values or
           rescale them, as the rule-based models are expecting these
           varaiances to be computed on the scaled reflectances. */
        variance (band_arr, 1.0, refl_input->refl_fill, refl_input->nlines,
            refl_input->nsamps, variance_arr);

        /* Write variance to the output file */
        if (put_output_lines (cm_output, variance_arr, ib+VARIANCE_B1, 0,
            refl_input->nlines, sizeof (int32)) != SUCCESS)
        {
            sprintf (errmsg, "Error writing variance for band %d of "
                "reflectance file", ib);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end for ib */

    /* Free the band array */
    free (band_arr);

    /* Allocate memory for the index array, whole band */
    index_arr = calloc (refl_input->nlines * refl_input->nsamps,
        sizeof (int16));
    if (index_arr == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the index array");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* index bands */
    for (ib = CM_NDVI; ib <= CM_NDSI; ib++)
    {
        /* Read the index band from the output file */
        if (verbose && ib == CM_NDVI)
            printf ("    - NDVI band\n");
        else if (verbose && ib == CM_NDSI)
            printf ("    - NDSI band\n");
        if (get_output_lines (cm_output, ib, 0, refl_input->nlines,
            sizeof(int16), index_arr) != SUCCESS)
        {
            sprintf (errmsg, "Error reading band %d of indices", ib);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Compute the variance.  The indices will need to be unscaled to
           compute the variance and then rescaled to store as an int. */
        variance (index_arr, SCALE_FACTOR, refl_input->refl_fill,
            refl_input->nlines, refl_input->nsamps, variance_arr);

        /* Write variance to the output file */
        if (put_output_lines (cm_output, variance_arr, ib+VARIANCE_NDVI, 0,
            refl_input->nlines, sizeof (int32)) != SUCCESS)
        {
            sprintf (errmsg, "Error writing variance for band %d of "
                "reflectance file", ib);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end for ib */

    /* Free the index and variance arrays */
    free (variance_arr);
    free (index_arr);

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Computing variance -- complete\n");

    /* Allocate memory for the NDVI and NDSI to hold one line of data */
    ndvi = calloc (refl_input->nsamps, sizeof (int16));
    if (ndvi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the NDVI");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    ndsi = calloc (refl_input->nsamps, sizeof (int16));
    if (ndsi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the NDSI");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for the variances to hold one line of data */
    ndvi_var = calloc (refl_input->nsamps, sizeof (int32));
    if (ndvi_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the NDVI variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    ndsi_var = calloc (refl_input->nsamps, sizeof (int32));
    if (ndsi_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the NDSI variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    b1_var = calloc (refl_input->nsamps, sizeof (int32));
    if (b1_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band1 variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    b2_var = calloc (refl_input->nsamps, sizeof (int32));
    if (b2_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band2 variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    b3_var = calloc (refl_input->nsamps, sizeof (int32));
    if (b3_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band3 variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    b4_var = calloc (refl_input->nsamps, sizeof (int32));
    if (b4_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band4 variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    b5_var = calloc (refl_input->nsamps, sizeof (int32));
    if (b5_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band5 variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    b7_var = calloc (refl_input->nsamps, sizeof (int32));
    if (b7_var == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the band7 variance");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for the revised cfmask arrays, for one line */
    rev_cm = calloc (refl_input->nsamps, sizeof (uint8));
    if (rev_cm == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the revised cloud mask");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    rev_lim_cm = calloc (refl_input->nsamps, sizeof (uint8));
    if (rev_lim_cm == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the revised cloud mask");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Running the rule-based models\n");

    /* Run the rule-based models one line at a time to compute the revised
       cloud masks */
    for (line = 0; line < refl_input->nlines; line++)
    {
        /* Read the reflectance bands */
        for (ib = 0; ib < refl_input->nrefl_band; ib++)
        {
            if (get_input_refl_lines (refl_input, ib, line, 1, NULL) != SUCCESS)
            {
                sprintf (errmsg, "Error reading reflectance band %d", ib);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Read the cfmask */
        if (get_input_cfmask_lines (refl_input, line, 1, NULL) != SUCCESS)
        {
            sprintf (errmsg, "Error reading reflectance band %d", ib);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Read the index bands */
        for (ib = CM_NDVI; ib <= CM_NDSI; ib++)
        {
            /* Read the index band from the output file */
            switch (ib)
            {
                case CM_NDVI:
                    ptr = ndvi; break;
                case CM_NDSI:
                    ptr = ndsi; break;
            }
            if (get_output_lines (cm_output, ib, line, 1, sizeof(int16), ptr)
                != SUCCESS)
            {
                sprintf (errmsg, "Error reading band %d of indices", ib);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Read the variance bands */
        for (ib = VARIANCE_B1; ib <= VARIANCE_NDSI; ib++)
        {
            /* Read the index or variance band from the output file */
            switch (ib)
            {
                case VARIANCE_B1:
                    ptr = b1_var; break;
                case VARIANCE_B2:
                    ptr = b2_var; break;
                case VARIANCE_B3:
                    ptr = b3_var; break;
                case VARIANCE_B4:
                    ptr = b4_var; break;
                case VARIANCE_B5:
                    ptr = b5_var; break;
                case VARIANCE_B7:
                    ptr = b7_var; break;
                case VARIANCE_NDVI:
                    ptr = ndvi_var; break;
                case VARIANCE_NDSI:
                    ptr = ndsi_var; break;
            }
            if (get_output_lines (cm_output, ib, line, 1, sizeof(int32), ptr)
                != SUCCESS)
            {
                sprintf (errmsg, "Error reading band %d of variances", ib);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Run the rule-based models on all the input data */
        rule_based_model (refl_input, ndsi, ndvi, b1_var, b2_var, b4_var,
            b5_var, b7_var, ndvi_var, ndsi_var, refl_input->nsamps, rev_cm,
            rev_lim_cm);

        /* Write the revised cloud masks */
        if (put_output_lines (cm_output, rev_cm, REVISED_CM, line, 1,
            sizeof (uint8)) != SUCCESS)
        {
            sprintf (errmsg, "Writing revised cloud mask for line %d", line);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        if (put_output_lines (cm_output, rev_lim_cm, REVISED_LIM_CM, line, 1,
            sizeof (uint8)) != SUCCESS)
        {
            sprintf (errmsg, "Writing revised limited cloud mask for line %d",
                line);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end for line */

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Rule-based models -- complete\n");

    /* Free the index and variance pointers */
    free (ndvi);
    free (ndsi);
    free (b1_var);
    free (b2_var);
    free (b3_var);
    free (b4_var);
    free (b5_var);
    free (b7_var);
    free (ndvi_var);
    free (ndsi_var);
    free (rev_cm);
    free (rev_lim_cm);

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Running the erosion and dilation filters on the revised "
            "cloud mask\n");

    /* Set up a 5x5 kernel for applying the erosion and dilation */
    element = cvCreateStructuringElementEx (5, 5, 1, 1, CV_SHAPE_RECT, NULL);
    if (element == NULL)
    {
        sprintf (errmsg, "Unable to establish a 5x5 kernel for erosion and "
            "dilation.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Create the OpenCV image to be used for the erosion and dilation */
    cv_img = cvCreateImage (cvSize (cm_output->nsamps, cm_output->nlines),
        IPL_DEPTH_8U, 1);
    if (cv_img == NULL)
    {
        sprintf (errmsg, "Unable to create the OpenCV image");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    opencv_img = (uint8 *)cv_img->imageData;

    /* Allocate memory for the revised cloud mask */
    rev_cm = calloc (refl_input->nlines * refl_input->nsamps, sizeof (uint8));
    if (rev_cm == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the revised cloud mask");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate memory for the limited revised cloud mask */
    rev_lim_cm = calloc (refl_input->nlines * refl_input->nsamps,
        sizeof (uint8));
    if (rev_lim_cm == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the limited revised "
            "cloud mask");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the revised cloud mask */
    if (get_output_lines (cm_output, REVISED_CM, 0, cm_output->nlines,
        sizeof(uint8), rev_cm) != SUCCESS)
    {
        sprintf (errmsg, "Reading the revised cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Loop through the revised cloud mask and put it in the OpenCV image */
    for (line = 0; line < cm_output->nlines; line++)
        for (samp = 0; samp < cm_output->nsamps; samp++)
        {
            opencv_img[line*cv_img->widthStep + samp] =
                rev_cm[line*cm_output->nsamps + samp];
        }

    /* Apply the erosion and dilation filters to the revised cloud mask */
    cvErode (cv_img, cv_img, element, 1);
    cvDilate (cv_img, cv_img, element, 1);

    /* Loop through OpenCV image and put it back in the revised cloud mask */
    opencv_img = (uint8 *)cv_img->imageData;
    for (line = 0; line < cm_output->nlines; line++)
        for (samp = 0; samp < cm_output->nsamps; samp++)
        {
            rev_cm[line*cm_output->nsamps + samp] =
                opencv_img[line*cv_img->widthStep + samp];
        }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Erosion and dilation -- complete\n");

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Running the erosion and dilation filters on the revised "
            "limited cloud mask\n");

    /* Read the limited cloud mask */
    if (get_output_lines (cm_output, REVISED_LIM_CM, 0, cm_output->nlines,
        sizeof(uint8), rev_lim_cm) != SUCCESS)
    {
        sprintf (errmsg, "Reading the revised limited cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Loop through the revised limited cloud mask and put it in the OpenCV
       image */
    for (line = 0; line < cm_output->nlines; line++)
        for (samp = 0; samp < cm_output->nsamps; samp++)
        {
            opencv_img[line*cv_img->widthStep + samp] =
                rev_lim_cm[line*cm_output->nsamps + samp];
        }

    /* Apply the erosion and dilation filters to the revised limited cloud
       mask.  Given there are more false clouds (speckles), let's use a
       2-pass erosion followed by dilation. */
    cvErode (cv_img, cv_img, NULL, 2);
    cvDilate (cv_img, cv_img, NULL, 2);

    /* Loop through OpenCV image and put it back in the revised cloud mask
       for buffering */
    for (line = 0; line < cm_output->nlines; line++)
        for (samp = 0; samp < cm_output->nsamps; samp++)
        {
            rev_lim_cm[line*cm_output->nsamps + samp] =
                opencv_img[line*cv_img->widthStep + samp];
        }

    /* Release the OpenCV image and structuring element */
    cvReleaseImage (&cv_img);
    cvReleaseStructuringElement (&element);
    
    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Erosion and dilation -- complete\n");

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Buffering the revised cloud mask\n");

    /* Allocate memory for the buffered cloud mask */
    buff_cm = calloc (refl_input->nlines * refl_input->nsamps, sizeof (uint8));
    if (buff_cm == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the buffered cloud mask");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Apply a 7 pixel buffer to all cloudy pixels in the revised cloud mask */
    if (buffer (rev_cm, 6, cm_output->nlines, cm_output->nsamps, buff_cm) !=
        SUCCESS)
    {
        sprintf (errmsg, "Buffering revised cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the revised buffered cloud mask */
    if (put_output_lines (cm_output, buff_cm, REVISED_CM, 0, cm_output->nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing revised cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Buffering the limited revised cloud mask\n");

    /* Apply a 7 pixel buffer to all cloudy pixels in the limited revised
       cloud mask */
    if (buffer (rev_lim_cm, 6, cm_output->nlines, cm_output->nsamps, buff_cm) !=
        SUCCESS)
    {
        sprintf (errmsg, "Buffering limited revised cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the revised buffered cloud mask */
    if (put_output_lines (cm_output, buff_cm, REVISED_LIM_CM, 0,
        cm_output->nlines, sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing limited revised cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Cloud buffering -- complete\n");

    /* Free the revised and buffered cloud masks */
    free (rev_cm);
    free (rev_lim_cm);
    free (buff_cm);

    /* Close the reflectance product */
    close_input (refl_input);
    free_input (refl_input);

    /* Write the ENVI header for spectral indices files */
    for (ib = 0; ib < cm_output->nband; ib++)
    {
        /* Create the ENVI header file this band */
        if (create_envi_struct (&cm_output->metadata.band[ib],
            &xml_metadata.global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
  
        /* Write the ENVI header */
        strcpy (envi_file, cm_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }
  
    /* Append the spectral index bands to the XML file */
    if (append_metadata (cm_output->nband, cm_output->metadata.band,
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending revised cloud mask bands to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
  
    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the output revised cloud mask product */
    close_output (cm_output);
    free_output (cm_output);

    /* Free the filename pointers */
    free (xml_infile);

    /* Indicate successful completion of processing */
    printf ("Revised cloud mask processing complete!\n");
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/20/2014    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("revised_cloud_mask revises the cfmask cloud mask in the input "
            "surface reflectance product, using the TOA reflectance bands. The "
            "cfmask flags many snow pixels as cloud pixels and this revision "
            "will clean up the cloud mask to correctly identify the snow "
            "pixels.\n\n");
    printf ("usage: revised_cloud_mask "
            "--xml=input_xml_filename [--verbose]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file to be processed\n");

    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");
    printf ("\nrevised_cloud_mask --help will print the usage statement\n");
    printf ("\nExample: revised_cloud_mask "
            "--xml=LT50400331995173AAA02.xml --verbose\n");
}

