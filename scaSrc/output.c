#include <time.h>
#include "output.h"

#define SDS_PREFIX ("band")

#define OUTPUT_PROVIDER ("DataProvider")
#define OUTPUT_SAT ("Satellite")
#define OUTPUT_INST ("Instrument")
#define OUTPUT_ACQ_DATE ("AcquisitionDate")
#define OUTPUT_L1_PROD_DATE ("Level1ProductionDate")
#define OUTPUT_SUN_ZEN ("SolarZenith")
#define OUTPUT_SUN_AZ ("SolarAzimuth")
#define OUTPUT_WRS_SYS ("WRS_System")
#define OUTPUT_WRS_PATH ("WRS_Path")
#define OUTPUT_WRS_ROW ("WRS_Row")
#define OUTPUT_NBAND ("NumberOfBands")
#define OUTPUT_SHORT_NAME ("ShortName")
#define OUTPUT_LOCAL_GRAN_ID ("LocalGranuleID")
#define OUTPUT_PROD_DATE ("SCAProductionDate")
#define OUTPUT_SCAVERSION ("SCAVersion")

#define OUTPUT_WEST_BOUND  ("WestBoundingCoordinate")
#define OUTPUT_EAST_BOUND  ("EastBoundingCoordinate")
#define OUTPUT_NORTH_BOUND ("NorthBoundingCoordinate")
#define OUTPUT_SOUTH_BOUND ("SouthBoundingCoordinate")

#define OUTPUT_LONG_NAME        ("long_name")
#define OUTPUT_UNITS            ("units")
#define OUTPUT_VALID_RANGE      ("valid_range")
#define OUTPUT_FILL_VALUE       ("_FillValue")
#define OUTPUT_SATU_VALUE       ("_SaturateValue")
#define OUTPUT_SCALE_FACTOR     ("scale_factor")
#define OUTPUT_ADD_OFFSET       ("add_offset")
#define OUTPUT_SCALE_FACTOR_ERR ("scale_factor_err")
#define OUTPUT_ADD_OFFSET_ERR   ("add_offset_err")
#define OUTPUT_CALIBRATED_NT    ("calibrated_nt")
#define OUTPUT_QAMAP_INDEX      ("qa_bitmap_index")

/******************************************************************************
MODULE:  create_output

PURPOSE:  Create a new HDF output file

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred creating the HDF file
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/12/2012    Gail Schmidt     Original Development (based on input routines
                              from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int create_output
(
    char *file_name    /* I: name of HDF file to be created */
)
{
    int32 hdf_file_id;        /* HDF file ID */
    char FUNC_NAME[] = "create_output";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
  
    /* Create the file with HDF open */
    hdf_file_id = Hopen (file_name, DFACC_CREATE, DEF_NDDS); 
    if (hdf_file_id == HDF_ERROR)
    {
        sprintf (errmsg, "Error creating the HDF file: %s", file_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Close the file */
    Hclose (hdf_file_id);
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  open_output

PURPOSE:  Set up the output data structure.  Open the output file for write
    access, and create the output Science Data Set (SDS).

RETURN VALUE:
Type = Output_t
Value          Description
-----          -----------
NULL           Error occurred opening the HDF file
Valid output   Successful completion
data structure

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/12/2012    Gail Schmidt     Original Development (based on input routines
                              from the LEDAPS lndsr application)

NOTES:
  1. Don't allocate space for buf, since pointers to existing buffers will
     be assigned in the output structure.
******************************************************************************/
Output_t *open_output
(
    char *file_name,                /* I: name of output HDF file */
    int nband,                      /* I: number of image bands (SDSs) to be
                                          created */
    char *sds_names[NUM_OUT_SDS],   /* I: array of SDS names for each band */
    int nlines,                     /* I: number of lines in image */
    int nsamps                      /* I: number of samples in image */
)
{
    Output_t *this = NULL;
    char FUNC_NAME[] = "create_output";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    Myhdf_dim_t *dim[MYHDF_MAX_RANK];     /* dimension information */
    Myhdf_sds_t *sds = NULL;  /* SDS information */
    int ir;    /* looping variable for rank/dimension */
    int ib;    /* looping variable for bands */

    /* Check parameters */
    if (nlines < 1)
    {
        sprintf (errmsg, "Invalid number of input lines (< 1)");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (nsamps < 1)
    {
        sprintf (errmsg, "Invalid number of input samples (< 1)");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (nband < 1 || nband > NUM_OUT_SDS)
    {
        sprintf (errmsg, "Invalid number of image bands");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    /* Create the Output data structure */
    this = (Output_t *) malloc (sizeof(Output_t));
    if (this == NULL) 
    {
        sprintf (errmsg, "Error allocating Output data structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    /* Populate the data structure */
    this->file_name = dup_string (file_name);
    if (this->file_name == (char *)NULL)
    {
        free (this);
        sprintf (errmsg, "Error duplicating file name");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    this->open = false;
    this->nband = nband;
    this->size.l = nlines;
    this->size.s = nsamps;
    for (ib = 0; ib < this->nband; ib++)
    {
        this->sds[ib].name = NULL;
        this->sds[ib].dim[0].name = NULL;
        this->sds[ib].dim[1].name = NULL;
        this->buf[ib] = NULL;
    }
  
    /* Open file for SD access */
    this->sds_file_id = SDstart ((char *)file_name, DFACC_RDWR);
    if (this->sds_file_id == HDF_ERROR)
    {
        free (this->file_name);
        free (this);  
        sprintf (errmsg, "Error opening output file for SD access");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
    this->open = true;
  
    /* Set up the image SDSs */
    for (ib = 0; ib < this->nband; ib++)
    {
        sds = &this->sds[ib];
        sds->rank = 2;
        sds->type = DFNT_UINT8;
        sds->name = dup_string (sds_names[ib]);
        if (sds->name == NULL)
        {
            free_output (this);
            close_output (this);
            sprintf (errmsg, "Error duplicating SDS name");
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
    
        dim[0] = &sds->dim[0];
        dim[1] = &sds->dim[1];
    
        dim[0]->nval = this->size.l;
        dim[1]->nval = this->size.s;
    
        dim[0]->type = dim[1]->type = sds->type;
    
        dim[0]->name = dup_string("YDim_Grid");
        if (dim[0]->name == NULL)
        {
            free_output (this);
            close_output (this);
            sprintf (errmsg, "Error duplicating y dim name");
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        dim[1]->name = dup_string("XDim_Grid");
        if (dim[1]->name == NULL)
        {
            free_output (this);
            close_output (this);
            sprintf (errmsg, "Error duplicating x dim name");
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
    
        if (put_sds_info (this->sds_file_id, sds) != SUCCESS)
        {
            free_output (this);
            close_output (this);
            sprintf (errmsg, "Error setting up the SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
    
        for (ir = 0; ir < sds->rank; ir++)
        {
            if (put_sds_dim_info (sds->id, ir, dim[ir]) != SUCCESS)
            {
                free_output (this);
                close_output (this);
                sprintf (errmsg, "Error setting up the dimension");
                error_handler (true, FUNC_NAME, errmsg);
                return (NULL);
            }
        }
    }  /* end for image bands */
  
    return this;
}


/******************************************************************************
MODULE:  close_output

PURPOSE:  Ends SDS access and closes the output file.

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

    /* Close image SDSs */
    for (ib = 0; ib < this->nband; ib++)
    {
        if (SDendaccess(this->sds[ib].id) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error ending SDS access for band %d.", ib);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF file itself */
    SDend (this->sds_file_id);
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

NOTES:
******************************************************************************/
int free_output
(
    Output_t *this    /* I/O: Output data structure to free */
)
{
    char FUNC_NAME[] = "free_output";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int ir;    /* looping variable for rank/dimension */
    int ib;    /* looping variable for bands */
  
    if (this->open) 
    {
        sprintf (errmsg, "File is still open, so cannot free memory.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (this != NULL)
    {
        /* Free image band SDSs */
        for (ib = 0; ib < this->nband; ib++)
        {
            for (ir = 0; ir < this->sds[ib].rank; ir++)
            {
                if (this->sds[ib].dim[ir].name != NULL) 
                    free (this->sds[ib].dim[ir].name);
            }
            if (this->sds[ib].name != NULL) 
                free (this->sds[ib].name);
        }
    
        if (this->file_name != NULL)
            free (this->file_name);

        free (this);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_output_line

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
2/12/2012    Gail Schmidt     Original Development (based on input routines
                              from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int put_output_line
(
    Output_t *this,    /* I: Output data structure; buf contains the line to
                             be written */
    int iband,         /* I: current band to be written (0-based) */
    int iline,         /* I: current line to be written (0-based) */
    int nlines         /* I: number of lines to be written */
)
{
    char FUNC_NAME[] = "put_output_line";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    int32 start[MYHDF_MAX_RANK];  /* starting location for writing to HDF */
    int32 nval[MYHDF_MAX_RANK];   /* number of values to write to HDF */
    void *buf = NULL;             /* buffer to write to HDF */
  
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
    if (iline < 0 || iline >= this->size.l)
    {
        sprintf (errmsg, "Invalid line number.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (nlines < 0 || iline+nlines > this->size.l)
    {
        sprintf (errmsg, "Line plus number of lines to be written exceeds "
            "the predefined size of the image.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Write the data */
    start[0] = iline;
    start[1] = 0;
    nval[0] = nlines;
    nval[1] = this->size.s;
    buf = (void *)this->buf[iband];
    if (SDwritedata (this->sds[iband].id, start, NULL, nval, buf) == HDF_ERROR)
    {
        sprintf (errmsg, "Error writing the output line(s) to HDF.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    
    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_metadata

PURPOSE:  Writes metadata to the output file.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred writing the metadata
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/12/2012    Gail Schmidt     Original Development (based on input routines
                              from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int put_metadata
(
    Output_t *this,      /* I: Output data structure */
    int nband,           /* I: number of bands to write */
    char *band_names[NUM_OUT_SDS],  /* I: band names to write */
    char *QA_on[NUM_OUT_SDS],  /* I: metadata info for the current band "on" */
    char *QA_off[NUM_OUT_SDS], /* I: metadata info for the current band "off" */
    Input_meta_t *meta         /* I: metadata to be written */
)
{
    char FUNC_NAME[] = "put_metadata";   /* function name */
    char errmsg[STR_SIZE];        /* error message */
    Myhdf_attr_t attr;            /* HDF attributes */
    char date[MAX_DATE_LEN + 1];  /* date string */
    char prod_date[MAX_DATE_LEN + 1];  /* production date for SCA products */
    double dval[NUM_OUT_SDS];     /* data value to write */
    char string[250];             /* string to be written */
    char long_name[250];          /* long name for the attribute */
    char process_ver[100];        /* SCA processing version */

    int ib;                       /* looping variable for bands */
    char* units_b=NULL;           /* units string for current band */
    char* qa_string=NULL;         /* QA string */
    char sca_QAMAP[1000];         /* QA values for the SCA product */
    time_t tp;                    /* structure for obtaining current time */
    struct tm *tm = NULL;         /* structure for obtaining current time
                                     in UTC format */
  
    /* Check the parameters */
    if (!this->open)
    {
        sprintf (errmsg, "File not open therefore cannot write metadata.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (nband < 1 || nband > NUM_OUT_SDS)
    {
        sprintf (errmsg, "Invalid number of bands.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Write the metadata */
    attr.id = -1;
  
    strcpy (string, "USGS/EROS");
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_PROVIDER;
    if (put_attr_string (this->sds_file_id, &attr, string) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (data provider)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    sprintf (string, "%s", meta->sat);
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_SAT;
    if (put_attr_string (this->sds_file_id, &attr, string) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (satellite)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
 
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(meta->inst);
    attr.name = OUTPUT_INST;
    if (put_attr_string (this->sds_file_id, &attr, meta->inst) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (instrument)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (format_date (&meta->acq_date, DATE_FORMAT_DATEA_TIME, date) != SUCCESS)
    {
        sprintf (errmsg, "Error formatting acquisition date");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    attr.type = DFNT_CHAR8;
    attr.nval = strlen(date);
    attr.name = OUTPUT_ACQ_DATE;
    if (put_attr_string (this->sds_file_id, &attr, date) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (acquisition date)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (format_date (&meta->prod_date, DATE_FORMAT_DATEA_TIME, date) != SUCCESS)
    {
        sprintf (errmsg, "Error formatting production date");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    attr.type = DFNT_CHAR8;
    attr.nval = strlen(date);
    attr.name = OUTPUT_L1_PROD_DATE;
    if (put_attr_string (this->sds_file_id, &attr, date) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (production date)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    attr.type = DFNT_FLOAT32;
    attr.nval = 1;
    attr.name = OUTPUT_SUN_ZEN;
    dval[0] = (double)meta->solar_elev * DEG;
    if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (solar zenith)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    attr.type = DFNT_FLOAT32;
    attr.nval = 1;
    attr.name = OUTPUT_SUN_AZ;
    dval[0] = (double)meta->solar_az * DEG;
    if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (solar azimuth)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    sprintf (string, "%s", meta->wrs_sys);
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_WRS_SYS;
    if (put_attr_string (this->sds_file_id, &attr, string) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (WRS system)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    attr.type = DFNT_INT16;
    attr.nval = 1;
    attr.name = OUTPUT_WRS_PATH;
    dval[0] = (double)meta->path;
    if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (WRS path)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    attr.type = DFNT_INT16;
    attr.nval = 1;
    attr.name = OUTPUT_WRS_ROW;
    dval[0] = (double)meta->row;
    if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (WRS row)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get the current time for the production date, convert it to UTC, and
       format it based on other date formats */
    if (time (&tp) == -1)
    {
        sprintf (errmsg, "Error obtaining the current time");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    tm = gmtime (&tp);
    if (tm == (struct tm *) NULL)
    {
        sprintf (errmsg, "Error converting current time to UTC");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (strftime (prod_date, (MAX_DATE_LEN+1), "%Y-%m-%dT%H:%M:%SZ", tm) == 0)
    {
        sprintf (errmsg, "Error formatting the production date and time");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    attr.type = DFNT_CHAR8;
    attr.nval = strlen (prod_date);
    attr.name = OUTPUT_PROD_DATE;
    if (put_attr_string (this->sds_file_id, &attr, prod_date) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (SCA production date)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    sprintf (process_ver, "%s", SCA_VERSION);
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(process_ver);
    attr.name = OUTPUT_SCAVERSION;
    if (put_attr_string (this->sds_file_id, &attr, process_ver) != SUCCESS)
    {
        sprintf (errmsg, "Error writing attribute (SCA Version)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* output the geographic bounding coordinates if they are available */
    if (!meta->bounds.is_fill)
    {
        attr.type = DFNT_FLOAT64;
        attr.nval = 1;
        attr.name = OUTPUT_WEST_BOUND;
        dval[0] = meta->bounds.min_lon;
        if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (West Bounding Coord)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
  
        attr.type = DFNT_FLOAT64;
        attr.nval = 1;
        attr.name = OUTPUT_EAST_BOUND;
        dval[0] = meta->bounds.max_lon;
        if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (East Bounding Coord)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
  
        attr.type = DFNT_FLOAT64;
        attr.nval = 1;
        attr.name = OUTPUT_NORTH_BOUND;
        dval[0] = meta->bounds.max_lat;
        if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (North Bounding Coord)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
  
        attr.type = DFNT_FLOAT64;
        attr.nval = 1;
        attr.name = OUTPUT_SOUTH_BOUND;
        dval[0] = meta->bounds.min_lat;
        if (put_attr_double (this->sds_file_id, &attr, dval) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (South Bounding Coord)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }  /* if geographic bounds are not fill */

    /* now write out the per sds attributes */
    for (ib = 0; ib < nband; ib++)
    {
        sprintf (long_name,"%s", band_names[ib]);
        attr.type = DFNT_CHAR8;
        attr.nval = strlen(long_name);
        attr.name = OUTPUT_LONG_NAME;
        if (put_attr_string (this->sds[ib].id, &attr, long_name) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (long name)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
      
        attr.type = DFNT_CHAR8;
        units_b=dup_string("quality/feature classification");
        attr.nval = strlen(units_b);
        attr.name = OUTPUT_UNITS;
        if (put_attr_string (this->sds[ib].id, &attr, units_b) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (units ref)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
      
        attr.type = DFNT_CHAR8;
        sprintf (sca_QAMAP,
            "\n\tQA pixel values are either off or on:\n"
            "\tValue  Description\n"
            "\t0\t%s\n"
            "\t255\t%s", QA_off[ib], QA_on[ib]);
        qa_string=dup_string(sca_QAMAP);
        attr.nval = strlen(qa_string);
        attr.name = "QA index";
        if (put_attr_string (this->sds[ib].id, &attr, qa_string) != SUCCESS)
        {
            sprintf (errmsg, "Error writing attribute (QA index)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }  /* end for ib */
  
    return (SUCCESS);
}
