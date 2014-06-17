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

/* Prototypes */
void usage ();

short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **xml_infile,    /* O: address of input XML file */
    bool *verbose         /* O: verbose flag */
);

void make_index
(
    int16 *band1,         /* I: input array of scaled reflectance data for
                                the spectral index */
    int16 *band2,         /* I: input array of scaled reflectance data for
                                the spectral index */
    int fill_value,       /* I: fill value for the reflectance values */
    int satu_value,       /* I: saturation value for the reflectance values */
    int nlines,           /* I: number of lines in the data arrays */
    int nsamps,           /* I: number of samples in the data arrays */
    float *spec_indx      /* O: output spectral index */
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
    float *ndsi_arr,        /* I: NDSI scaled values */
    float *ndvi_arr,        /* I: NDVI scaled values */
    float *b1_var_arr,      /* I: band1 variance values */
    float *b2_var_arr,      /* I: band2 variance values */
    float *b4_var_arr,      /* I: band4 variance values */
    float *b5_var_arr,      /* I: band5 variance values */
    float *b7_var_arr,      /* I: band7 variance values */
    float *ndvi_var_arr,    /* I: NDVI variance values */
    float *ndsi_var_arr,    /* I: NDSI variance values */
    int nsamps,             /* I: number of samples in the input arrays */
    uint8 *rev_cloud_mask,      /* O: revised cloud mask */
    uint8 *rev_lim_cloud_mask   /* O: revised cloud mask without variances */
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
