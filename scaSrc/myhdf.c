#include <stdlib.h>
#include "myhdf.h"
#include "error_handler.h"
#include "hdf.h"
#include "mfhdf.h"
#include "mystring.h"

/* Constants */
#define DIM_MAX_NCHAR (80)  /* Maximum size of a dimension name */

/* Possible ranges for data types */
#define MYHDF_CHAR8H     (        255  )
#define MYHDF_CHAR8L     (          0  )
#define MYHDF_INT8H      (        127  )
#define MYHDF_INT8L      (       -128  )
#define MYHDF_UINT8H     (        255  )
#define MYHDF_UINT8L     (          0  )
#define MYHDF_INT16H     (      32767  )
#define MYHDF_INT16L     (     -32768  )
#define MYHDF_UINT16H    (      65535u )
#define MYHDF_UINT16L    (          0u )
#define MYHDF_INT32H     ( 2147483647l )
#define MYHDF_INT32L     (-2147483647l )
#define MYHDF_UINT32H    ( 4294967295ul)
#define MYHDF_UINT32L    (          0ul)
#define MYHDF_FLOAT32H   (3.4028234e+38f)
#define MYHDF_FLOAT32L   (1.1754943e-38f)
#define MYHDF_FLOAT64H   (1.797693134862316e+308)
#define MYHDF_FLOAT64L   (2.225073858507201e-308)

/******************************************************************************
MODULE:  get_sds_info

PURPOSE:  Reads the HDF information for the given SDS.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading information for this SDS
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int get_sds_info
(
    int32 sds_file_id,   /* I: HDF file ID */
    Myhdf_sds_t *sds     /* I/O: attributes for the SDS; SDS name needs to
                                 be populated */
)
{
    char FUNC_NAME[] = "get_sds_info";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int32 dims[MYHDF_MAX_RANK];
  
    sds->index = SDnametoindex (sds_file_id, sds->name);
    if (sds->index == HDF_ERROR)
    {
        sprintf (errmsg, "Error getting index of SDS %s", sds->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    sds->id = SDselect (sds_file_id, sds->index);
    if (sds->id == HDF_ERROR)
    {
        sprintf (errmsg, "Error selecting ID for SDS %s", sds->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (SDgetinfo(sds->id, sds->name, &sds->rank, dims, 
                  &sds->type, &sds->nattr) == HDF_ERROR)
    {
        SDendaccess(sds->id);
        sprintf (errmsg, "Error getting SDS information for %s", sds->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (sds->rank > MYHDF_MAX_RANK)
    {
        SDendaccess(sds->id);
        sprintf (errmsg, "SDS rank too large for %s", sds->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_sds_dim_info

PURPOSE:  Reads the dimension information for the given SDS.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading information for this SDS
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int get_sds_dim_info
(
    int32 sds_id,        /* I: SDS ID */
    int irank,           /* I: rank to get dimension information for */
    Myhdf_dim_t *dim     /* O: dimension structure to populate */
)
{
    char FUNC_NAME[] = "get_sds_dim_info";   /* function name */
    char errmsg[STR_SIZE];          /* error message */
    char dim_name[DIM_MAX_NCHAR];   /* dimension name */
  
    dim->id = SDgetdimid(sds_id, irank);
    if (dim->id == HDF_ERROR) 
    {
        sprintf (errmsg, "Error getting dimension ID for dimension %d", irank);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (SDdiminfo (dim->id, dim_name, &dim->nval, &dim->type, &dim->nattr)
        == HDF_ERROR)
    {
        sprintf (errmsg, "Error getting dimension information");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    dim->name = dup_string (dim_name);
    if (dim->name == (char *)NULL)
    {
        sprintf (errmsg, "Error copying dimension name %s", dim_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_sds_info

PURPOSE:  Creates an SDS and writes the SDS information.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred creating/writing data for this SDS
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int put_sds_info
(
    int32 sds_file_id,   /* I: HDF file ID */
    Myhdf_sds_t *sds     /* I/O: attributes for the SDS to be written; and
                                 attributes to be obtained */
)
{
    char FUNC_NAME[] = "put_sds_info";   /* function name */
    char errmsg[STR_SIZE];          /* error message */
    int irank;                      /* looping index for dimensions */
    int32 dims[MYHDF_MAX_RANK];     /* array of dimensions */
  
    /* Copy the dimensions for the SDS */
    for (irank = 0; irank < sds->rank; irank++)
        dims[irank] = sds->dim[irank].nval;
  
    /* Create the SDS */
    sds->id = SDcreate (sds_file_id, sds->name, sds->type, sds->rank, dims);
    if (sds->id == HDF_ERROR)
    {
        sprintf (errmsg, "Error creating the SDS %s", sds->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    sds->index = SDnametoindex (sds_file_id, sds->name);
    if (sds->index == HDF_ERROR)
    {
        sprintf (errmsg, "Error getting the index of SDS %s", sds->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_sds_dim_info

PURPOSE:  Writes information for the SDS dimension.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred writing data for this SDS
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int put_sds_dim_info
(
    int32 sds_id,        /* I: SDS ID */
    int irank,           /* I: rank to get dimension information for */
    Myhdf_dim_t *dim     /* I/O: dimension structure to populate */
)
{
    char FUNC_NAME[] = "put_sds_dim_info";   /* function name */
    char errmsg[STR_SIZE];          /* error message */

    dim->id = SDgetdimid (sds_id, irank);
    if (dim->id == HDF_ERROR) 
    {
        sprintf (errmsg, "Error getting the dimension ID for dimension %d",
            irank);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (SDsetdimname (dim->id, dim->name) == HDF_ERROR)
    {
        sprintf (errmsg, "Error setting the dimension name for dimension %d",
            irank);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Set dimension type */
      /* !! do it !! */
 
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_attr_double

PURPOSE:  Reads an attribute of the defined data type into a parameter of type
double.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading this attribute
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int get_attr_double
(
    int32 sds_id,         /* I: SDS ID to read attribute from */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    double *val           /* O: array of values returned as doubles */
)
{
    char FUNC_NAME[] = "get_attr_double";   /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variable for the number of values */
    char attr_name[80];      /* attribute name */
    char8 val_char8[MYHDF_MAX_NATTR_VAL];
    uint8 val_int8[MYHDF_MAX_NATTR_VAL];
    uint8 val_uint8[MYHDF_MAX_NATTR_VAL];
    int16 val_int16[MYHDF_MAX_NATTR_VAL];
    uint16 val_uint16[MYHDF_MAX_NATTR_VAL];
    int32 val_int32[MYHDF_MAX_NATTR_VAL];
    uint32 val_uint32[MYHDF_MAX_NATTR_VAL];
    float32 val_float32[MYHDF_MAX_NATTR_VAL];
    float64 val_float64[MYHDF_MAX_NATTR_VAL];
    
    /* Get the attribute ID and attribute information */
    if ((attr->id = SDfindattr(sds_id, attr->name)) == HDF_ERROR)
    {
        sprintf (errmsg, "Error getting attribute ID for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (SDattrinfo(sds_id, attr->id, attr_name, &attr->type, &attr->nval) == 
        HDF_ERROR)
    {
        sprintf (errmsg, "Error getting attribute info for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Validate the number of attribute values */
    if (attr->nval < 1)
    {
        sprintf (errmsg, "No attribute values for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (attr->nval > MYHDF_MAX_NATTR_VAL) 
    {
        sprintf (errmsg, "Too many attribute values for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the attribute based on its data type */
    switch (attr->type)
    {
    case DFNT_CHAR8:
        if (SDreadattr (sds_id, attr->id, val_char8) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (char8)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_char8[i];
        break;

    case DFNT_INT8:
        if (SDreadattr(sds_id, attr->id, val_int8) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (int8)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_int8[i];
        break;

    case DFNT_UINT8:
        if (SDreadattr(sds_id, attr->id, val_uint8) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (uint8)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_uint8[i];
        break;

    case DFNT_INT16:
        if (SDreadattr(sds_id, attr->id, val_int16) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (int16)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_int16[i];
        break;

    case DFNT_UINT16:
        if (SDreadattr(sds_id, attr->id, val_uint16) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (uint16)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_uint16[i];
        break;

    case DFNT_INT32:
        if (SDreadattr(sds_id, attr->id, val_int32) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (int32)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
          val[i] = (double)val_int32[i];
        break;

    case DFNT_UINT32:
        if (SDreadattr(sds_id, attr->id, val_uint32) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (uint32)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_uint32[i];
        break;

    case DFNT_FLOAT32:
        if (SDreadattr(sds_id, attr->id, val_float32) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (float32)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_float32[i];
        break;

    case DFNT_FLOAT64:
        if (SDreadattr(sds_id, attr->id, val_float64) == HDF_ERROR) 
        {
            sprintf (errmsg, "Error reading attribute (float64)");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        for (i = 0; i < attr->nval; i++) 
            val[i] = (double)val_float64[i];
        break;

    default:
        sprintf (errmsg, "Unknown attribute data type");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_attr_double

PURPOSE:  Writes an attribute from a parameter type of double to the HDF file.
The double value is converted to the native data type before writing.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred writing this attribute
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int put_attr_double
(
    int32 sds_id,         /* I: SDS ID to write attribute to */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    double *val           /* I: array of values to be written as native type */
)
{
    char FUNC_NAME[] = "put_attr_double";   /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variable for the number of values */
    void *buf = NULL;        /* void pointer to actual data array */
    char8 val_char8[MYHDF_MAX_NATTR_VAL];
    int8 val_int8[MYHDF_MAX_NATTR_VAL];
    uint8 val_uint8[MYHDF_MAX_NATTR_VAL];
    int16 val_int16[MYHDF_MAX_NATTR_VAL];
    uint16 val_uint16[MYHDF_MAX_NATTR_VAL];
    int32 val_int32[MYHDF_MAX_NATTR_VAL];
    uint32 val_uint32[MYHDF_MAX_NATTR_VAL];
    float32 val_float32[MYHDF_MAX_NATTR_VAL];
    float64 val_float64[MYHDF_MAX_NATTR_VAL];

    if (attr->nval <= 0 || attr->nval > MYHDF_MAX_NATTR_VAL) 
    {
        sprintf (errmsg, "Invalid number of attribute values");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Write the attribute based on its data type */
    switch (attr->type)
    {
    case DFNT_CHAR8:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_CHAR8H))
                val_char8[i] = MYHDF_CHAR8H;
            else if (val[i] <= ((double)MYHDF_CHAR8L))
                val_char8[i] = MYHDF_CHAR8L;
            else if (val[i] >= 0.0)
                val_char8[i] = (char8)(val[i] + 0.5);
            else
                val_char8[i] = -((char8)(-val[i] + 0.5));
        }
        buf = (void *)val_char8;
        break;
  
    case DFNT_INT8:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_INT8H))
                val_int8[i] = MYHDF_INT8H;
            else if (val[i] <= ((double)MYHDF_INT8L))
                val_int8[i] = MYHDF_INT8L;
            else if (val[i] >= 0.0)
                val_int8[i] = (int8)(val[i] + 0.5);
            else
                val_int8[i] = -((int8)(-val[i] + 0.5));
        }
        buf = (void *)val_int8;
        break;
  
    case DFNT_UINT8:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_UINT8H))
                val_uint8[i] = MYHDF_UINT8H;
            else if (val[i] <= ((double)MYHDF_UINT8L))
                val_uint8[i] = MYHDF_UINT8L;
            else if (val[i] >= 0.0)
                val_uint8[i] = (uint8)(val[i] + 0.5);
            else
                val_uint8[i] = -((uint8)(-val[i] + 0.5));
        }
        buf = (void *)val_uint8;
        break;
  
    case DFNT_INT16:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_INT16H))
                val_int16[i] = MYHDF_INT16H;
            else if (val[i] <= ((double)MYHDF_INT16L))
                val_int16[i] = MYHDF_INT16L;
            else if (val[i] >= 0.0)
                val_int16[i] = (int16)( val[i] + 0.5);
            else
                val_int16[i] = -((int16)(-val[i] + 0.5));
        }
        buf = (void *)val_int16;
        break;
  
    case DFNT_UINT16:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_UINT16H))
                val_uint16[i] = MYHDF_UINT16H;
            else if (val[i] <= ((double)MYHDF_UINT16L))
                val_uint16[i] = MYHDF_UINT16L;
            else if (val[i] >= 0.0)
                val_uint16[i] = (uint16)( val[i] + 0.5);
            else
                val_uint16[i] = -((uint16)(-val[i] + 0.5));
        }
        buf = (void *)val_uint16;
        break;
  
    case DFNT_INT32:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_INT32H))
                val_int32[i] = MYHDF_INT32H;
            else if (val[i] <= ((double)MYHDF_INT32L))
                val_int32[i] = MYHDF_INT32L;
            else if (val[i] >= 0.0)
                val_int32[i] = (int32)( val[i] + 0.5);
            else
                val_int32[i] = -((int32)(-val[i] + 0.5));
        }
        buf = (void *)val_int32;
        break;
  
    case DFNT_UINT32:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_UINT32H))
                val_uint32[i] = MYHDF_UINT32H;
            else if (val[i] <= ((double)MYHDF_UINT32L))
                val_uint32[i] = MYHDF_UINT32L;
            else if (val[i] >= 0.0)
                val_uint32[i] = (uint32)( val[i] + 0.5);
            else
                val_uint32[i] = -((uint32)(-val[i] + 0.5));
        }
        buf = (void *)val_uint32;
        break;
  
    case DFNT_FLOAT32:
        for (i = 0; i < attr->nval; i++)
        {
            if (val[i] >= ((double)MYHDF_FLOAT32H))
                val_float32[i] = MYHDF_FLOAT32H;
            else if (val[i] <= ((double)MYHDF_FLOAT32L))
                val_float32[i] = MYHDF_FLOAT32L;
            else if (val[i] >= 0.0)
                val_float32[i] = (float32)( val[i] + 0.5);
            else
                val_float32[i] = -((float32)(-val[i] + 0.5));
        }
        buf = (void *)val_float32;
        break;
  
    case DFNT_FLOAT64:
        if (sizeof (float64) == sizeof (double))
            buf = (void *)val;
        else
        {
            for (i = 0; i < attr->nval; i++)
                val_float64[i] = val[i];
            buf = (void *)val_float64;
        }
        break;
  
    default: 
        sprintf (errmsg, "Unsupported attribute data type");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (SDsetattr(sds_id, attr->name, attr->type, attr->nval, buf) == HDF_ERROR)
    {
        sprintf (errmsg, "Error writing attribute");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_attr_string

PURPOSE:  Reads a string attribute.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading this attribute
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int get_attr_string
(
    int32 sds_id,         /* I: SDS ID to read attribute from */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    char *string          /* O: attribute string read from SDS */
)
{
    char FUNC_NAME[] = "get_attr_string";   /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variable for the number of values */
    char attr_name[80];      /* attribute name */
    int i_length;            /* length of string of this attribute */
    void *buf;               /* void pointer to buffer for reading */
    char8 val_char8[MYHDF_MAX_NATTR_VAL];
  
    /* Get the attribute ID and attribute information */
    if ((attr->id = SDfindattr(sds_id, attr->name)) == HDF_ERROR)
    {
        sprintf (errmsg, "Error getting attribute ID for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (SDattrinfo(sds_id, attr->id, attr_name, &attr->type, &attr->nval) == 
        HDF_ERROR)
    {
        sprintf (errmsg, "Error getting attribute info for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Validate the number of attribute values */
    if (attr->nval < 1)
    {
        sprintf (errmsg, "No attribute values for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (attr->nval > MYHDF_MAX_NATTR_VAL) 
    {
        sprintf (errmsg, "Too many attribute values for %s", attr->name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Validate the data type is character */
    if (attr->type != DFNT_CHAR8) 
    {
        sprintf (errmsg, "Invalid data type - should be string (char8)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Set up the void buffer and read the attribute */
    if (sizeof(char8) == sizeof(char))
        buf = (void *)string;
    else
        buf = (void *)val_char8;

    if (SDreadattr (sds_id, attr->id, buf) == HDF_ERROR) 
    {
        sprintf (errmsg, "Error reading attribute");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (sizeof(char8) != sizeof(char))
    {
        for (i = 0; i < attr->nval; i++) 
            string[i] = (char)val_char8[i];
    }

    /* Terminate the string since HDF doesn't do that correctly */
    i_length= (int)attr->nval;
    string[i_length]= '\0';

    return (SUCCESS);
}


/******************************************************************************
MODULE:  put_attr_string

PURPOSE:  Writes a string attribute to the HDF file.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred writing this attribute
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int put_attr_string
(
    int32 sds_id,         /* I: SDS ID to write attribute to */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    char * string         /* I: string value to be written */
)
{
    char FUNC_NAME[] = "put_attr_string";   /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variable for the number of values */
    void *buf = NULL;        /* void pointer to actual data array */
    char8 val_char8[MYHDF_MAX_NATTR_VAL];

    if (attr->nval <= 0 || attr->nval > MYHDF_MAX_NATTR_VAL) 
    {
        sprintf (errmsg, "Invalid number of attribute values");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Validate the data type is character */
    if (attr->type != DFNT_CHAR8) 
    {
        sprintf (errmsg, "Invalid data type - should be string (char8)");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Set up the void buffer and write the attribute */
    if (sizeof(char8) == sizeof(char))
        buf = (void *)string;
    else
    {
        for (i = 0; i < attr->nval; i++) 
            val_char8[i] = (char8)string[i];
        buf = (void *)val_char8;
    }

    if (SDsetattr(sds_id, attr->name, attr->type, attr->nval, buf) == HDF_ERROR)
    {
        sprintf (errmsg, "Error writing attribute");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}
