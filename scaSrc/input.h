#ifndef _INPUT_H_
#define _INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bool.h"
#include "mystring.h"
#include "myhdf.h"
#include "error_handler.h"
#include "date.h"

/* Conversion factor for degrees to radians */
#ifndef PI
#ifndef M_PI
#define PI (3.141592653589793238)
#else
#define PI (M_PI)
#endif
#endif

#define TWO_PI (2.0 * PI)
#define HALF_PI (PI / 2.0)

#define DEG (180.0 / PI)
#define RAD (PI / 180.0)

/* There are currently a maximum of 6 reflective bands in the output surface
   reflectance product */
#define NBAND_REFL_MAX 6

/* How many lines of TOA reflectance and brightness temperature data should be
   processed at one time */
#define PROC_NLINES 10

/* How many lines of DEM data should be processed at one time, multiple of 3
   since the shade relief input needs to be a factor of 3. */
#define DEM_PROC_NLINES 300

/* Structure for the global metadata */
typedef struct {
    char provider[STR_SIZE]; /* data provider type */
    char sat[STR_SIZE];      /* satellite */
    char inst[STR_SIZE];     /* instrument */
    Date_t acq_date;         /* acqsition date/time (scene center) */
    Date_t prod_date;        /* production date */
    float solar_elev;        /* solar elevation angle (radians; scene center) */
    float solar_az;          /* solar azimuth angle (radians; scene center) */
    char wrs_sys[STR_SIZE];  /* WRS system */
    int path;                /* WRS path number */
    int row;                 /* WRS row number */
    float pixsize;           /* pixel size */
    int refl_band[NBAND_REFL_MAX]; /* band numbers for TOA reflectance data */
    int btemp_band;          /* band number for brightness temp data */
} Input_meta_t;

/* Structure for the 'input' data type, particularly to handle the file/SDS
   IDs and the band-specific information */
typedef struct {
    Input_meta_t meta;       /* input metadata */
    char *refl_file_name;    /* input TOA reflectance file name */
    bool refl_open;          /* open reflectance file flag; open = true */
    char *btemp_file_name;   /* input brightness temp file name */
    bool btemp_open;         /* open brightness temp file flag; open = true */
    int nrefl_band;          /* number of input TOA reflectance bands */
    int nbtemp_band;         /* number of input brightness temp bands */
    int nlines;              /* number of input lines */
    int nsamps;              /* number of input samples */
    int32 refl_sds_file_id;  /* SDS file id for TOA reflectance */
    int32 btemp_sds_file_id; /* SDS file id for brightness temp */
    Myhdf_sds_t refl_sds[NBAND_REFL_MAX]; /* SDS data structures for TOA
                                reflectance data */
    int16 *refl_buf[NBAND_REFL_MAX]; /* input data buffer for unscaled TOA
                                reflectance data (PROC_NLINES lines of data) */
    Myhdf_sds_t btemp_sds;   /* SDS data structure for brightness temp data */
    int16 *btemp_buf;        /* input data buffer for unscaled brightness temp
                                data (PROC_NLINES lines of thermal data) */
    int refl_fill;           /* fill value for TOA reflectance bands */
    int btemp_fill;          /* fill value for brightness temperature band */
    float refl_scale_fact;   /* scale factor for TOA reflectance bands */
    float btemp_scale_fact;  /* scale factor for brightness temp bands */
    int refl_saturate_val;   /* saturation value for TOA reflectance bands */
    int btemp_saturate_val;  /* saturation value for brightness temp bands */
} Input_t;

/* Prototypes */
Input_t *open_input
(
    char *refl_file_name,     /* I: input TOA reflectance filename */
    char *btemp_file_name     /* I: input brightness temp filename */
);

void close_input
(
    Input_t *this    /* I: pointer to input data structure */
);

void free_input
(
    Input_t *this    /* I: pointer to input data structure */
);

int get_input_refl_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current TOA band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines       /* I: number of lines to read */
);

int get_input_btemp_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iline,       /* I: current line to read (0-based) */
    int nlines       /* I: number of lines to read */
);

int get_input_meta
(
    Input_t *this    /* I: pointer to input data structure */
);

#endif
