#include "input.h"

/******************************************************************************
MODULE:  open_input

PURPOSE:  Sets up the input data structure, opens the input reflectance file
for read access, allocates space, and stores some of the metadata for later
reference.

RETURN VALUE:
Type = Input_t*
Value      Description
-----      -----------
NULL       Error occurred opening or reading the file
non-NULL   Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)
6/9/2015     Gail Schmidt     Updated the cfmask band name to 'cfmask' since
                              it had changed in ESPA a while ago

NOTES:
  1. This routine opens the input reflectance files.  It also allocates memory
     for pointers in the input structure.  It is up to the caller to use
     close_input and free_input to close the files and free up the memory when
     done using the input data structure.
  2. The read buffers are only set up to read NPROC_LINES at a time.
******************************************************************************/
Input_t *open_input
(
    Espa_internal_meta_t *metadata,     /* I: input metadata */
    bool toa         /* I: are we processing TOA reflectance data, otherwise
                           process surface reflectance data */
)
{
    char FUNC_NAME[] = "open_input";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    Input_t *this = NULL;     /* input data structure to be initialized,
                                 populated, and returned to the caller */
    int ib;                   /* loop counter for bands */
    int refl_indx = -1;       /* band index in XML file for the reflectance
                                 band */
    int16 *buf = NULL;        /* temporary buffer to allocate memory for
                                 the reflectance bands */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */
  
    /* Create the Input data structure */
    this = malloc (sizeof (Input_t));
    if (this == NULL) 
    {
        strcpy (errmsg, "Error allocating memory for Input data structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Initialize the input pointers */
    this->refl_open = false;
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
    {
        this->file_name[ib] = NULL;
        this->fp_bin[ib] = NULL;
        this->refl_buf[ib] = NULL;
    }
    this->cfmask_file_name = NULL;
    this->fp_cfmask = NULL;
    this->cfmask_buf = NULL;

    /* Initialize the input fields using information from the metadata
       structure */
    if (!strcmp (gmeta->instrument, "TM") ||
        !strncmp (gmeta->instrument, "ETM", 3))
    {
        this->nrefl_band = 6;     /* number of reflectance bands */
        this->refl_band[0] = 1;
        this->refl_band[1] = 2;
        this->refl_band[2] = 3;
        this->refl_band[3] = 4;
        this->refl_band[4] = 5;
        this->refl_band[5] = 7;
    }
    else
    {
        free_input (this);
        sprintf (errmsg, "Unsupported instrument type.  Currently only "
            "TM and ETM+ are supported");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* If processing TOA bands, then grab the toa_refl product bands.
       Otherwise grab the surface reflectance bands. */
    if (toa)
    {
        for (ib = 0; ib < metadata->nbands; ib++)
        {
            if (!strcmp (metadata->band[ib].name, "toa_band1") &&
                !strcmp (metadata->band[ib].product, "toa_refl"))
            {
                /* this is the index we'll use for reflectance band info */
                refl_indx = ib;
                this->file_name[0] = strdup (metadata->band[ib].file_name);
            }
            else if (!strcmp (metadata->band[ib].name, "toa_band2") &&
                !strcmp (metadata->band[ib].product, "toa_refl"))
                this->file_name[1] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "toa_band3") &&
                !strcmp (metadata->band[ib].product, "toa_refl"))
                this->file_name[2] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "toa_band4") &&
                !strcmp (metadata->band[ib].product, "toa_refl"))
                this->file_name[3] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "toa_band5") &&
                !strcmp (metadata->band[ib].product, "toa_refl"))
                this->file_name[4] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "toa_band7") &&
                !strcmp (metadata->band[ib].product, "toa_refl"))
                this->file_name[5] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "cfmask") &&
                !strcmp (metadata->band[ib].product, "cfmask"))
                this->cfmask_file_name = strdup (metadata->band[ib].file_name);
        }  /* for ib */
    }
    else
    {
        for (ib = 0; ib < metadata->nbands; ib++)
        {
            if (!strcmp (metadata->band[ib].name, "sr_band1") &&
                !strcmp (metadata->band[ib].product, "sr_refl"))
            {
                /* this is the index we'll use for reflectance band info */
                refl_indx = ib;
                this->file_name[0] = strdup (metadata->band[ib].file_name);
            }
            else if (!strcmp (metadata->band[ib].name, "sr_band2") &&
                !strcmp (metadata->band[ib].product, "sr_refl"))
                this->file_name[1] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "sr_band3") &&
                !strcmp (metadata->band[ib].product, "sr_refl"))
                this->file_name[2] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "sr_band4") &&
                !strcmp (metadata->band[ib].product, "sr_refl"))
                this->file_name[3] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "sr_band5") &&
                !strcmp (metadata->band[ib].product, "sr_refl"))
                this->file_name[4] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "sr_band7") &&
                !strcmp (metadata->band[ib].product, "sr_refl"))
                this->file_name[5] = strdup (metadata->band[ib].file_name);
            else if (!strcmp (metadata->band[ib].name, "cfmask") &&
                !strcmp (metadata->band[ib].product, "cfmask"))
                this->cfmask_file_name = strdup (metadata->band[ib].file_name);
        }  /* for ib */
    }

    /* Make sure we found the bands */
    if (refl_indx == -1)
    {
        free_input (this);
        if (toa)
            sprintf (errmsg, "Unable to find the TOA reflectance bands in the "
                "XML file.");
        else
            sprintf (errmsg, "Unable to find the surface reflectance bands in "
                "the XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Pull the reflectance info from representative band1 in the XML file */
    this->nsamps = metadata->band[refl_indx].nsamps;
    this->nlines = metadata->band[refl_indx].nlines;
    this->pixsize[0] = metadata->band[refl_indx].pixel_size[0];
    this->pixsize[1] = metadata->band[refl_indx].pixel_size[1];
    this->refl_fill = metadata->band[refl_indx].fill_value;
    this->refl_scale_fact = metadata->band[refl_indx].scale_factor;
    this->refl_saturate_val = metadata->band[refl_indx].saturate_value;

    /* Open each of the reflectance files then the cfmask file */
    for (ib = 0; ib < this->nrefl_band; ib++)
    {
        this->fp_bin[ib] = open_raw_binary (this->file_name[ib], "rb");
        if (this->fp_bin[ib] == NULL)
        {
            free_input (this);
            sprintf (errmsg, "Opening raw binary file: %s",
                this->file_name[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
    }

    this->fp_cfmask = open_raw_binary (this->cfmask_file_name, "rb");
    if (this->fp_cfmask == NULL)
    {
        free_input (this);
        sprintf (errmsg, "Opening raw binary file: %s", this->cfmask_file_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
    this->refl_open = true;

    /* Validate the input data type */
    if (metadata->band[refl_indx].data_type != ESPA_INT16)
    {
        close_input (this);
        free_input (this);
        sprintf (errmsg, "Input data type is assumed to be int16, but "
            "the data in the XML file for the reflectance bands doesn't "
            "match this data type.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Allocate input buffer.  Reflectance buffer has multiple bands.
       Allocate PROC_NLINES of data for each band. */
    buf = calloc (PROC_NLINES * this->nsamps * this->nrefl_band,
        sizeof (int16));
    if (buf == NULL)
    {
        close_input (this);
        free_input (this);
        sprintf (errmsg, "Allocating memory for input reflectance buffer "
            "containing %d lines.", PROC_NLINES);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
    else
    {
        /* Set up the memory buffers for each band */
        this->refl_buf[0] = buf;
        for (ib = 1; ib < this->nrefl_band; ib++)
            this->refl_buf[ib] = this->refl_buf[ib-1] +
                PROC_NLINES * this->nsamps;
    }

    this->cfmask_buf = calloc (PROC_NLINES * this->nsamps, sizeof (uint8));
    if (this->cfmask_buf == NULL)
    {
        close_input (this);
        free_input (this);
        sprintf (errmsg, "Allocating memory for input cfmask buffer "
            "containing %d lines.", PROC_NLINES);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Success */
    return (this);
}


/******************************************************************************
MODULE:  close_input

PURPOSE:  Ends SDS access and closes the input file.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
******************************************************************************/
void close_input
(
    Input_t *this    /* I: pointer to input data structure */
)
{
    int ib;      /* loop counter for bands */
  
    /* Close the raw binary files */
    if (this->refl_open)
    {
        /* Close reflectance and cfmask files */
        for (ib = 0; ib < this->nrefl_band; ib++)
            close_raw_binary (this->fp_bin[ib]);
        close_raw_binary (this->fp_cfmask);
        this->refl_open = false;
    }
}


/******************************************************************************
MODULE:  free_input

PURPOSE:  Frees memory in the input data structure.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
******************************************************************************/
void free_input
(
    Input_t *this    /* I: pointer to input data structure */
)
{
    char FUNC_NAME[] = "free_input";   /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int ib;                            /* loop counter for bands */
   
    if (this != NULL)
    {
        if (this->refl_open) 
        {
            strcpy (errmsg, "Freeing input data structure, but reflectance "
                "file is still open. Use close_input to close the file");
            error_handler (false, FUNC_NAME, errmsg);
        }
  
        /* Free image band pointers */
        for (ib = 0; ib < this->nrefl_band; ib++)
            free (this->file_name[ib]);
        free (this->cfmask_file_name);
  
        /* Free the data buffers */
        free (this->refl_buf[0]);
        free (this->cfmask_buf);

        /* Free the data structure */
        free (this);
    } /* end if */
}


/******************************************************************************
MODULE:  get_input_refl_lines

PURPOSE:  Reads the reflectance data for the current band and lines, and
populates the refl_buf buffer in the Input_t data structure.  If the output_arr
is not NULL, read the data into the output_arr array instead of the refl_buf
buffer.  This allows the caller to read more data than the basic NPROC_LINES
that has been allocated by default.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_refl_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int16 *out_arr   /* O: output array to populate, if not NULL */
)
{
    char FUNC_NAME[] = "get_input_refl_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
    void *buf = NULL;         /* pointer to the buffer for the current band */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->refl_open)
    {
        strcpy (errmsg, "Reflectance file has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nrefl_band)
    {
        strcpy (errmsg, "Invalid band number for the reflectance file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->nlines)
    {
        strcpy (errmsg, "Invalid line number for reflectance band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* If the output array is not NULL then use it, otherwise use the refl_buf
       in the Input_t data structure */
    if (out_arr == NULL)
        buf = (void *) this->refl_buf[iband];
    else
        buf = (void *) out_arr;

    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->nsamps * sizeof (int16);
    if (fseek (this->fp_bin[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin[iband], nlines, this->nsamps,
        sizeof (int16), buf) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from reflectance band %d starting "
            "at line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_cfmask_lines

PURPOSE:  Reads the cfmask data for the current lines, and populates the
cfmask_buf buffer in the Input_t data structure.  If the output_arr is not
NULL, read the data into the output_arr array instead of the refl_buf buffer.
This allows the caller to read more data than the basic NPROC_LINES that has
been allocated by default.

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
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_cfmask_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int16 *out_arr   /* O: output array to populate, if not NULL */
)
{
    char FUNC_NAME[] = "get_input_cfmask_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
    void *buf = NULL;         /* pointer to the buffer for the current band */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->refl_open)
    {
        strcpy (errmsg, "Reflectance file has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->nlines)
    {
        strcpy (errmsg, "Invalid line number for cfmask band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* If the output array is not NULL then use it, otherwise use the refl_buf
       in the Input_t data structure */
    if (out_arr == NULL)
        buf = (void *) this->cfmask_buf;
    else
        buf = (void *) out_arr;

    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->nsamps * sizeof (uint8);
    if (fseek (this->fp_cfmask, loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input cfmask file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_cfmask, nlines, this->nsamps, sizeof (uint8),
        buf) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from cfmask band starting at line "
            "%d", nlines, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}

