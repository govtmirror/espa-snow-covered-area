#ifndef _COMMON_H
#define _COMMON_H

typedef signed short int16;
typedef unsigned char uint8;
typedef signed int int32;

/* Define the output products to be processed */
typedef enum {REVISED_CM=0, NUM_CM} Mycm_list_t;

/* Application version */
#define CLOUD_MASK_VERSION "1.3.0"

/* How many lines of data should be processed at one time */
#define PROC_NLINES 1000

#endif
