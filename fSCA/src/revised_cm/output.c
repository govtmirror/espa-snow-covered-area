#include <time.h>
#include <ctype.h>
#include "output.h"


/******************************************************************************
MODULE:  open_output

PURPOSE:  Set up the output data structure.  Open the output file for write
and read access.

RETURN VALUE:
Type = Output_t
Value          Description
-----          -----------
NULL           Error occurred opening the file
not-NULL       Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
  1. Don't allocate space for buf, since pointers to existing buffers will
     be assigned in the output structure.
******************************************************************************/
Output_t *open_output
(
    Espa_internal_meta_t *in_meta,  /* I: input metadata structure */
    Input_t *input,                 /* I: input reflectance band data */
    int nband,                      /* I: number of bands to be created */
    char short_names[][STR_SIZE],   /* I: array of short names for new bands */
    char long_names[][STR_SIZE],    /* I: array of long names for new bands */
    char data_units[][STR_SIZE],    /* I: array of data units for new bands */
    bool toa         /* I: are we processing TOA reflectance data, otherwise
                           process surface reflectance data */
)
{
    Output_t *this = NULL;
    char FUNC_NAME[] = "open_output";   /* function name */
    char errmsg[STR_SIZE];       /* error message */
    char *upper_str = NULL;      /* upper case version of the SI short name */
    char *mychar = NULL;         /* pointer to '_' */
    char scene_name[STR_SIZE];   /* scene name for the current scene */
    char production_date[MAX_DATE_LEN+1]; /* current date/time for production */
    time_t tp;                   /* time structure */
    struct tm *tm = NULL;        /* time structure for UTC time */
    int ib;    /* looping variable for bands */
    int refl_indx = -1;          /* band index in XML file for the reflectance
                                    band */
    Espa_band_meta_t *bmeta = NULL;  /* pointer to the band metadata array
                                        within the output structure */

    /* Check parameters */
    if (nband < 1 || nband > MAX_OUT_BANDS)
    {
        sprintf (errmsg, "Invalid number of image bands");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    /* Create the Output data structure */
    this = (Output_t *) malloc (sizeof (Output_t));
    if (this == NULL) 
    {
        sprintf (errmsg, "Error allocating Output data structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    /* Find the representative band for metadata information.  SR vs. TOA
       shouldn't matter in this case.  They should be the same size and
       resolution. */
    for (ib = 0; ib < in_meta->nbands; ib++)
    {
        if (!strcmp (in_meta->band[ib].name, "toa_band1") &&
            !strcmp (in_meta->band[ib].product, "toa_refl"))
        {
            /* this is the index we'll use for reflectance band info */
            refl_indx = ib;
            break;
        }
    }

    /* Make sure we found the TOA band 1 */
    if (refl_indx == -1)
    {
        sprintf (errmsg, "Unable to find the TOA reflectance bands in the "
            "XML file for initializing the output metadata.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Initialize the internal metadata for the output product. The global
       metadata won't be updated, however the band metadata will be updated
       and used later for appending to the original XML file. */
    init_metadata_struct (&this->metadata);

    /* Allocate memory for the total bands */
    if (allocate_band_metadata (&this->metadata, nband) != SUCCESS)
    {
        sprintf (errmsg, "Allocating band metadata.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
    bmeta = this->metadata.band;

    /* Determine the scene name */
    strcpy (scene_name, in_meta->band[refl_indx].file_name);
    mychar = strchr (scene_name, '_');
    if (mychar != NULL)
      *mychar = '\0';
  
    /* Get the current date/time (UTC) for the production date of each band */
    if (time (&tp) == -1)
    {
        sprintf (errmsg, "Unable to obtain the current time.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    tm = gmtime (&tp);
    if (tm == NULL)
    {
        sprintf (errmsg, "Converting time to UTC.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (strftime (production_date, MAX_DATE_LEN, "%Y-%m-%dT%H:%M:%SZ", tm) == 0)
    {
        sprintf (errmsg, "Formatting the production date/time.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Populate the data structure */
    this->open = false;
    this->nband = nband;
    this->nlines = input->nlines;
    this->nsamps = input->nsamps;
    for (ib = 0; ib < this->nband; ib++)
        this->fp_bin[ib] = NULL;
 
    for (ib = 0; ib < nband; ib++)
    {
        strncpy (bmeta[ib].short_name, in_meta->band[refl_indx].short_name, 3);
        bmeta[ib].short_name[3] = '\0';
        upper_str = upper_case_str (short_names[ib]);
        strcat (bmeta[ib].short_name, upper_str);
        strcpy (bmeta[ib].product, "fsca");
        if (toa)
            strcpy (bmeta[ib].source, "toa_refl");
        else
            strcpy (bmeta[ib].source, "sr_refl");
        strcpy (bmeta[ib].category, "index");
        bmeta[ib].nlines = this->nlines;
        bmeta[ib].nsamps = this->nsamps;
        bmeta[ib].pixel_size[0] = input->pixsize[0];
        bmeta[ib].pixel_size[1] = input->pixsize[1];
        strcpy (bmeta[ib].pixel_units, "meters");
        sprintf (bmeta[ib].app_version, "revised_cloud_mask_%s",
            CLOUD_MASK_VERSION);
        strcpy (bmeta[ib].production_date, production_date);
        strcpy (bmeta[ib].name, short_names[ib]);
        strcpy (bmeta[ib].long_name, long_names[ib]);
        strcpy (bmeta[ib].data_units, data_units[ib]);

        /* Handle the cloud mask bands differently */
        if (ib == REVISED_CM)
        {
            bmeta[ib].data_type = ESPA_UINT8;
            bmeta[ib].fill_value = CFMASK_FILL_VALUE;
            bmeta[ib].valid_range[0] = 0;
            bmeta[ib].valid_range[1] = 3;

            /* Set up class values information */
            if (allocate_class_metadata (&bmeta[ib], 4) != SUCCESS)
            {
                sprintf (errmsg, "Allocating cfmask classes.");
                error_handler (true, FUNC_NAME, errmsg);
                return (NULL);
            }
          
            /* Identify the class values for the mask */
            bmeta[ib].class_values[0].class = 0;
            bmeta[ib].class_values[1].class = 1;
            bmeta[ib].class_values[2].class = 2;
            bmeta[ib].class_values[3].class = 3;
            strcpy (bmeta[ib].class_values[0].description, "clear");
            strcpy (bmeta[ib].class_values[1].description, "cloud in cfmask");
            strcpy (bmeta[ib].class_values[2].description, "cloud");
            strcpy (bmeta[ib].class_values[3].description, "water");
        }

        /* Set up the filename with the scene name and band name and open the
           file for read/write access */
        sprintf (bmeta[ib].file_name, "%s_%s.img", scene_name, bmeta[ib].name);
        this->fp_bin[ib] = open_raw_binary (bmeta[ib].file_name, "w+");
        if (this->fp_bin[ib] == NULL)
        {
            sprintf (errmsg, "Unable to open output band %d file: %s", ib,
                bmeta[ib].file_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        /* Free the memory for the upper-case string */
        free (upper_str);
    }  /* for ib */
    this->open = true;

    /* Successful completion */
    return this;
}


/******************************************************************************
MODULE:  close_output

PURPOSE:  Closes the output file.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred closing the HDF file
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/12/2012    Gail Schmidt     Original Development (based on input routines
                              from the LEDAPS lndsr application)
2/14/2014    Gail Schmidt     Modified to work with ESPA internal raw binary
                              file format

NOTES:
******************************************************************************/
int close_output
(
    Output_t *this    /* I/O: Output data structure to close */
)
{
    char FUNC_NAME[] = "close_output";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int ib;                   /* looping variable */

    if (!this->open)
    {
        sprintf (errmsg, "File is not open, so it cannot be closed.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close raw binary products */
    for (ib = 0; ib < this->nband; ib++)
        close_raw_binary (this->fp_bin[ib]);
    this->open = false;

    return (SUCCESS);
}


/******************************************************************************
MODULE:  free_output

PURPOSE:  Frees the memory for the output data structure

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred freeing the data structure
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/12/2012    Gail Schmidt     Original Development (based on input routines
                              from the LEDAPS lndsr application)
2/14/2014    Gail Schmidt     Modified to work with ESPA internal raw binary
                              file format

NOTES:
******************************************************************************/
int free_output
(
    Output_t *this    /* I/O: Output data structure to free */
)
{
    char FUNC_NAME[] = "free_output";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
  
    if (this->open) 
    {
        sprintf (errmsg, "Spectral index file is still open, so cannot free "
            "memory.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (this != NULL)
    {
        /* Free the band data */
        free (this->metadata.band[REVISED_CM].class_values);
        free (this->metadata.band);

        /* Free the data structure */
        free (this);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_output_lines

PURPOSE:  Writes a line or lines of data to the output file.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred writing the output data
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
******************************************************************************/
int put_output_lines
(
    Output_t *this,    /* I: Output data structure; buf contains the line to
                             be written */
    void *buf,         /* I: buffer to be written */
    int iband,         /* I: current band to be written (0-based) */
    int iline,         /* I: current line to be written (0-based) */
    int nlines,        /* I: number of lines to be written */
    int nbytes         /* I: number of bytes per pixel in this band */
)
{
    char FUNC_NAME[] = "put_output_lines";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the output file */
  
    /* Check the parameters */
    if (this == (Output_t *)NULL) 
    {
        sprintf (errmsg, "Invalid input structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open)
    {
        sprintf (errmsg, "File is not open.  Cannot write data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband)
    {
        sprintf (errmsg, "Invalid band number.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->nlines)
    {
        sprintf (errmsg, "Invalid line number.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (nlines < 0 || iline+nlines > this->nlines)
    {
        sprintf (errmsg, "Line plus number of lines to be written exceeds "
            "the predefined size of the image.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Write the data, but first seek to the correct line */
    loc = (long) iline * this->nsamps * nbytes;
    if (fseek (this->fp_bin[iband], loc, SEEK_SET))
    {
        sprintf (errmsg, "Seeking to the current line in the output file for "
            "band %d", iband);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (write_raw_binary (this->fp_bin[iband], nlines, this->nsamps, nbytes,
        buf) != SUCCESS)
    {
        sprintf (errmsg, "Error writing the output line(s) for band %d.",
            iband);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_output_lines

PURPOSE:  Reads the cfmask data for the current lines, and populates the buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading cfmask data for this line
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
  1. The Output_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_output to do that.
******************************************************************************/
int get_output_lines
(
    Output_t *this,  /* I: pointer to output data structure */
    int iband,       /* I: current band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int nbytes,      /* I: number of bytes per pixel in this band */
    void *buf        /* I: pointer to the buffer to be returned */
)
{
    char FUNC_NAME[] = "get_output_lines";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open)
    {
        strcpy (errmsg, "Output file has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->nlines)
    {
        strcpy (errmsg, "Invalid line number for cfmask band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->nsamps * nbytes;
    if (fseek (this->fp_bin[iband], loc, SEEK_SET))
    {
        sprintf (errmsg, "Seeking to the current line in the output file for "
            "band %d", iband);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin[iband], nlines, this->nsamps, nbytes,
        buf) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from band %d starting at line "
            "%d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  upper_case_str

PURPOSE:  Returns the upper case version of the input string.

RETURN VALUE:
Type = char *
Value      Description
-----      -----------
up_str     Upper case version of the input string

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/14/2012    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
char *upper_case_str
(
    char *str    /* I: string to be converted to upper case */
)
{
    char *up_str = NULL;    /* upper case version of the input string */
    char *ptr = NULL;       /* pointer to the upper case string */

    up_str = strdup (str);
    ptr = up_str;
    while (*ptr != '\0')
    {
        if (islower (*ptr))
            *ptr = toupper (*ptr);
        ptr++;
    }

    return up_str;
}

