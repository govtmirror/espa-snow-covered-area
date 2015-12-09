#ifndef _REVISED_CLOUD_MASK_H_
#define _REVISED_CLOUD_MASK_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "common.h"
#include "input.h"
#include "output.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "error_handler.h"

/* Defines for cloud mask and revised cloud mask */
#define CFMASK_FILL 255       /* fill value in input cfmask */
#define CFMASK_CLOUD 4        /* cloud value in input cfmask */
#define OUT_NOCLOUD 0         /* pixel is not cloud (clear) */
#define OUT_POSS_CLOUD 1      /* pixel was cloud in input cfmask (possibly
                                 cloud) */
#define OUT_CLOUD 2           /* pixel is cloud in revised cfmask (cloud) */
#define OUT_WATER 3           /* pixel is water in DSWE */

/* Prototypes */
void usage ();
void version ();

short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **xml_infile,    /* O: address of input XML file */
    bool *verbose         /* O: verbose flag */
);

float make_index
(
    int16 band1,          /* I: input scaled reflectance value for the spectral
                                index */
    int16 band2,          /* I: input scaled reflectance value for the spectral
                                index */
    int fill_value,       /* I: fill value for the reflectance values */
    int satu_value        /* I: saturation value for the reflectance values */
);

void variance
(
    float *array,       /* I: input array of data for which to compute the
                              covariances */
    int fill_value,     /* I: fill value for the band */
    int nlines,         /* I: number of lines in the data arrays */
    int nsamps,         /* I: number of samples in the data arrays */
    float *variance     /* O: output variance array */
);

void rule_based_model
(
    Input_t *input_img,     /* I: pointer to input data structure containing
                                  the scaled reflectance and cloud mask
                                  buffers (reflectance bands are scaled) */
    int nsamps,             /* I: number of samples in the input arrays */
    uint8 *rev_cloud_mask   /* O: revised cloud mask */
);

short buffer
(
    uint8 *array,       /* I: input array of data for which to buffer by the
                              distance value */
    int distance,       /* I: distance to buffer */
    int nlines,         /* I: number of lines in the data arrays */
    int nsamps,         /* I: number of samples in the data arrays */
    uint8 *buff_array   /* O: output array with buffer applied */
);

#endif
