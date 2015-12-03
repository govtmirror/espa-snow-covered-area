#ifndef _INPUT_H_
#define _INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "error_handler.h"
#include "raw_binary_io.h"
#include "espa_metadata.h"

/* There are currently a maximum of 6 reflective bands in the output surface
   reflectance product.  The reflective bands plus the cfmask will be the bands
   which are input for this application. */
#define NBAND_REFL_MAX 6

/* Structure for the 'input' data type, particularly to handle the file/SDS
   IDs and the band-specific information */
typedef struct {
    bool refl_open;          /* open reflectance file flag; open = true */
    int nrefl_band;          /* number of input reflectance bands */
    int nlines;              /* number of input lines */
    int nsamps;              /* number of input samples */
    float pixsize[2];        /* pixel size x, y */
    int refl_band[NBAND_REFL_MAX];   /* band numbers for reflectance data */
    char *file_name[NBAND_REFL_MAX]; /* name of the input image files */
    int16 *refl_buf[NBAND_REFL_MAX]; /* input data buffer for unscaled
                                        reflectance and cfmask data
                                        (PROC_NLINES lines of data) */
    FILE *fp_bin[NBAND_REFL_MAX];    /* file pointer for binary files */
    char *cfmask_file_name;  /* name of the input cfmask files */
    uint8 *cfmask_buf;       /* input data buffer for cfmask data
                                (PROC_NLINES lines of data) */
    FILE *fp_cfmask;         /* file pointer for cfmask file */
    char *dswe_file_name;    /* name of the input dswe files */
    uint8 *dswe_buf;         /* input data buffer for dswe data */
    FILE *fp_dswe;           /* file pointer for dswe file */
    int16 refl_fill;         /* fill value for reflectance bands */
    float refl_scale_fact;   /* scale factor for reflectance bands */
    int refl_saturate_val;   /* saturation value for reflectance bands */
} Input_t;

/* Prototypes */
Input_t *open_input
(
    Espa_internal_meta_t *metadata,     /* I: input metadata */
    bool toa         /* I: are we processing TOA reflectance data, otherwise
                           process surface reflectance data */
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
    int iband,       /* I: current band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int16 *out_arr   /* O: output array to populate, if not NULL */
);

int get_input_cfmask_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int16 *out_arr   /* O: output array to populate, if not NULL */
);

int get_input_dswe_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int16 *out_arr   /* O: output array to populate, if not NULL */
);

bool open_cfmask_dswe
(
    Input_t *input_struct   /* I: existing input structure */
);

void close_cfmask_dswe
(
    Input_t *this    /* I: pointer to input data structure */
);

void free_cfmask_dswe
(
    Input_t *this    /* I: pointer to input data structure */
);

#endif
