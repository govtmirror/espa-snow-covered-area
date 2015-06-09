#ifndef _COMMON_H
#define _COMMON_H

typedef signed short int16;
typedef unsigned char uint8;
typedef signed int int32;

/* Define the output products to be processed */
typedef enum {CM_NDVI=0, CM_NDSI, VARIANCE_B1, VARIANCE_B2, VARIANCE_B3,
    VARIANCE_B4, VARIANCE_B5, VARIANCE_B7, VARIANCE_NDVI, VARIANCE_NDSI,
    REVISED_CM, REVISED_LIM_CM, NUM_CM} Mycm_list_t;

/* Application version */
#define CLOUD_MASK_VERSION "1.2.0"

/* How many lines of data should be processed at one time */
#define PROC_NLINES 1000

#endif
