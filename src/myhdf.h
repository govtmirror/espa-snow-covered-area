#ifndef MYHDF_H
#define MYHDF_H

#include "hdf.h"
#include "mfhdf.h"
#include "bool.h"

#define MYHDF_MAX_RANK (4)     /* maximum rank of an SDS expected */
#define MYHDF_MAX_NATTR_VAL (3000)
                               /* maximum number of attribute values expected */
#define HDF_ERROR (-1)

/* Structure to store information about the HDF SDS */
typedef struct
{
  int32 nval, id;        /* number of values and id */ 
  int32 type, nattr;     /* HDF data type and number of attributes */
  char *name;            /* dimension name */
} Myhdf_dim_t;

typedef struct
{
  int32 index, id, rank;           /* index, id and rank */
  int32 type, nattr;               /* HDF data type and number of attributes */
  char *name;                      /* SDS name */
  Myhdf_dim_t dim[MYHDF_MAX_RANK]; /* dimension data structure */
} Myhdf_sds_t;

/* Structure to store information about the HDF attribute */
typedef struct
{
  int32 id, type, nval;	 /* id, data type and number of values */
  char *name;            /* attribute name */
} Myhdf_attr_t;

/* Prototypes */
int get_sds_info
(
    int32 sds_file_id,   /* I: HDF file ID */
    Myhdf_sds_t *sds     /* I/O: attributes for the SDS; SDS name needs to
                                 be populated */
);

int get_sds_dim_info
(
    int32 sds_id,        /* I: SDS ID */
    int irank,           /* I: rank to get dimension information for */
    Myhdf_dim_t *dim     /* O: dimension structure to populate */
);

int put_sds_info
(
    int32 sds_file_id,   /* I: HDF file ID */
    Myhdf_sds_t *sds     /* I/O: attributes for the SDS to be written; and
                                 attributes to be obtained */
);

int put_sds_dim_info
(
    int32 sds_id,        /* I: SDS ID */
    int irank,           /* I: rank to get dimension information for */
    Myhdf_dim_t *dim     /* I/O: dimension structure to populate */
);

int get_attr_double
(
    int32 sds_id,         /* I: SDS ID to read attribute from */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    double *val           /* O: array of values returned as doubles */
);

int put_attr_double
(
    int32 sds_id,         /* I: SDS ID to write attribute to */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    double *val           /* I: array of values to be written as native type */
);

int get_attr_string
(
    int32 sds_id,         /* I: SDS ID to read attribute from */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    char *string          /* O: attribute string read from SDS */
);

int put_attr_string
(
    int32 sds_id,         /* I: SDS ID to write attribute to */
    Myhdf_attr_t *attr,   /* I: attribute data structure */
    char * string         /* I: string value to be written */
);

#endif
