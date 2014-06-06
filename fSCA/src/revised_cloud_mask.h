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
    int16 *spec_indx      /* O: output spectral index */
);

void variance
(
    int16 *array,       /* I: input array of data for which to compute the
                              covariances */
    float scale_factor, /* I: scale factor for the input data */
    int fill_value,     /* I: fill value for the band */
    int nlines,         /* I: number of lines in the data arrays */
    int nsamps,         /* I: number of samples in the data arrays */
    int32 *variance     /* O: output variance array; if the scale_factor is 1.0
                              then the outputs are not scaled, otherwise these
                              are scaled values */
);

void boosted_model
(
    Input_t *input_img,     /* I: pointer to input data structure containing
                                  the scaled reflectance and cloud mask
                                  buffers (reflectance bands are scaled) */
    int16 *ndsi_arr,        /* I: NDSI scaled values (scaled) */
    int16 *ndvi_arr,        /* I: NDVI scaled values (scaled) */
    int32 *b1_var_arr,      /* I: band1 variance values */
    int32 *b2_var_arr,      /* I: band2 variance values */
    int32 *b4_var_arr,      /* I: band4 variance values */
    int32 *b5_var_arr,      /* I: band5 variance values */
    int32 *b7_var_arr,      /* I: band7 variance values */
    int32 *ndvi_var_arr,    /* I: NDVI variance values (scaled) */
    int32 *ndsi_var_arr,    /* I: NDSI variance values (scaled) */
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
