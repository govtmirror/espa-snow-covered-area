#include <cv.h>
#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  revised_cloud_mask

PURPOSE:  Revised the cfmask cloud mask which is included as part of the ESPA
surface reflectance product.  The cfmask flags many snow pixels as cloud pixels
and this revision will clean up the cloud mask to correctly identify the snow
pixels.  This routine also reads the DSWE PSCCSS band and flags water pixels
in the output mask.

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
12/2/2015     Gail Schmidt     Mark any clear pixels that are flagged as cloud
                               in cfmask as possibly cloudy
12/2/2015     Gail Schmidt     Removed the mask computed using the variances
                               and stick to the better mask without variances
12/2/2015     Gail Schmidt     Commented out the buffer around the
                               erosion/dilation results because it seems to
                               better represent the clouds without it
12/4/2015     Gail Schmidt     Erosion/dilation fills in fill pixels.  Remask
                               all fill pixels as fill.

NOTES:
  1. The DSWE PSCCSS band is used and water-high confidence (value of 1),
     water-moderate confidence (value of 2), and partial surface water (value
     of 3) pixels are marked as water pixels in the output mask.
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
    int ib;                    /* looping variable for bands */
    int pix;                   /* current pixel */
    int line, samp;            /* current line,samp to be processed */
    int num_cm;                /* number of cloud mask products to be output */
    uint8 *rev_cm=NULL;        /* revised cloud mask */
#ifdef BUFFER
    uint8 *buff_cm=NULL;       /* revised cloud mask with buffering */
#endif
    uint8 *opencv_img=NULL;    /* pointer to the opencv image data */
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
    if (validate_xml_file (xml_infile) != SUCCESS)
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

    /* Set up the output information for the output cloud mask */
    num_cm = NUM_CM;
    strcpy (short_cm_names[REVISED_CM], "revcm");
    strcpy (long_cm_names[REVISED_CM], "revised cloud mask");
    strcpy (cm_data_units[REVISED_CM], "quality/feature classification");

    /* Open the specified output files and create the metadata structure */
    cm_output = open_output (&xml_metadata, refl_input, num_cm,
        short_cm_names, long_cm_names, cm_data_units, toa_refl);
    if (cm_output == NULL)
    {   /* error message already printed */
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
            sprintf (errmsg, "Error reading cfmask band");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Run the rule-based model for the revised cloud mask */
        rule_based_model (refl_input, refl_input->nsamps, rev_cm);

        /* Write the revised cloud mask */
        if (put_output_lines (cm_output, rev_cm, REVISED_CM, line, 1,
            sizeof (uint8)) != SUCCESS)
        {
            sprintf (errmsg, "Writing revised cloud mask for line %d", line);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end for line */

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Rule-based model -- complete\n");

    /* Close the reflectance product */
    close_input (refl_input);
    free_input (refl_input);

    /* Free the cloud mask pointer */
    free (rev_cm);

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
    {
        for (samp = 0; samp < cm_output->nsamps; samp++)
        {
            opencv_img[line*cv_img->widthStep + samp] =
                rev_cm[line*cm_output->nsamps + samp];
        }
    }

    /* Apply the erosion and dilation filters to the revised cloud mask */
    cvErode (cv_img, cv_img, element, 1);
    cvDilate (cv_img, cv_img, element, 1);

    /* Loop through OpenCV image and put it back in the revised cloud mask */
    opencv_img = (uint8 *)cv_img->imageData;
    for (line = 0; line < cm_output->nlines; line++)
    {
        for (samp = 0; samp < cm_output->nsamps; samp++)
        {
            rev_cm[line*cm_output->nsamps + samp] =
                opencv_img[line*cv_img->widthStep + samp];
        }
    }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Erosion and dilation -- complete\n");

    /* Release the OpenCV image and structuring element */
    cvReleaseImage (&cv_img);
    cvReleaseStructuringElement (&element);
    
#ifdef BUFFER
    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Buffering the revised cloud mask\n");
#endif

    /* Reopen the cfmask and open DSWE.  Allocate memory for both buffers to
       contain data for the entire image. */
    if (!open_cfmask_dswe (refl_input))
    {
        sprintf (errmsg, "Error opening the cfmask and dswe bands");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the cfmask for the entire image */
    if (get_input_cfmask_lines (refl_input, 0, refl_input->nlines, NULL) !=
        SUCCESS)
    {
        sprintf (errmsg, "Error reading cfmask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the DSWE for the entire image, if it's available */
    if (refl_input->dswe_file_name)
    {
        if (get_input_dswe_lines (refl_input, 0, refl_input->nlines, NULL) !=
            SUCCESS)
        {
            sprintf (errmsg, "Error reading DSWE band");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

#ifdef BUFFER
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
#endif

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Marking cfmask and dswe in revised cloud mask\n");

    /* Loop through the pixels and mark the clear pixels that are flagged as
       cloud in cfmask as possible cloud pixels */
    for (pix = 0; pix < refl_input->nlines * refl_input->nsamps; pix++)
    {
        /* Make sure fill pixels stay as fill pixels and the erosion/dilation
           doesn't overwrite fill pixels */
        if (refl_input->cfmask_buf[pix] == CFMASK_FILL)
            rev_cm[pix] = 0;

        /* Mark cfmask pixels */
        if (rev_cm[pix] == OUT_NOCLOUD &&
            refl_input->cfmask_buf[pix] == CFMASK_CLOUD)
            rev_cm[pix] = OUT_POSS_CLOUD;
    }

    /* Loop through the pixels and mark the water pixels from the DSWE band as
       water pixels, if hte DSWE band is available */
    if (refl_input->dswe_file_name)
    {
        for (pix = 0; pix < refl_input->nlines * refl_input->nsamps; pix++)
        {
            /* Mark cfmask pixels */
            if (refl_input->dswe_buf[pix] > 0 &&
                refl_input->dswe_buf[pix] <= 3)
                rev_cm[pix] = OUT_WATER;
        }
    }

    /* Write the revised buffered cloud mask */
    if (put_output_lines (cm_output, rev_cm, REVISED_CM, 0, cm_output->nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing revised cloud mask band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Marking cfmask/dswe -- complete\n");

#ifdef BUFFER
    /* Print the processing status if verbose */
    if (verbose)
        printf ("  Cloud buffering -- complete\n");
#endif

    /* Free the revised and buffered cloud masks */
    free (rev_cm);
#ifdef BUFFER
    free (buff_cm);
#endif

    /* Close and free the cfmask and dswe pointers and buffers, then free the
       input product */
    close_cfmask_dswe (refl_input);
    free_cfmask_dswe (refl_input);
    free (refl_input);

    /* Write the ENVI header for the cloud mask */
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
  
    /* Append the revised cloud mask band to the XML file */
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

