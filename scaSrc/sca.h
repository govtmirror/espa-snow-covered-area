#ifndef _SCA_H_
#define _SCA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bool.h"
#include "mystring.h"
#include "error_handler.h"
#include "input.h"

/* Set up the fill, cloud, snow, and deep shadow mask values */
#define NO_DATA 255
#define VALID_DATA 0
#define CLOUD_COVER 255
#define NO_CLOUD 0
#define SNOW_COVER 255
#define NO_SNOW 0
#define DEEP_SHADOW 255
#define NO_DEEP_SHADOW 0

/* Define the terrain-derived deep shadow threshold */
#define TERRAIN_DEEP_SHADOW_THRESH 0.03

/* Prototypes */
void usage ();

short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **toa_infile,    /* O: address of input TOA filename */
    char **btemp_infile,  /* O: address of input TOA filename */
    char **dem_infile,    /* O: address of input DEM filename */
    char **sc_outfile,    /* O: address of output snow cover filename */
    bool *verbose         /* O: verbose flag */
);

void cloud_cover_class
(
    int16 *b1,     /* I: array of unscaled band 1 TOA reflectance values */
    int16 *b4,     /* I: array of unscaled band 4 TOA reflectance values */
    int16 *b6,     /* I: array of unscaled band 6 brightness temp values */
    int16 *b7,     /* I: array of unscaled band 7 TOA reflectance values */
    int nlines,    /* I: number of lines in the data arrays */
    int nsamps,    /* I: number of samples in the data arrays */
    float refl_scale_fact,  /* I: scale factor for the TOA reflectance values */
    float btemp_scale_fact, /* I: scale factor for the broghtness temp values */
    uint8 *refl_qa_mask,  /* I: array of masked values for processing (non-zero
                                values are not to be processed) reflectance
                                bands */
    uint8 *therm_qa_mask, /* I: array of masked values for processing (non-zero
                                values are not to be processed) thermal bands */
    uint8 *cloud_mask     /* O: array of cloud cover masked values (non-zero
                                values represent clouds) */
);

void snow_cover_class
(
    int16 *b1,     /* I: array of unscaled band 1 TOA reflectance values */
    int16 *b2,     /* I: array of unscaled band 2 TOA reflectance values */
    int16 *b3,     /* I: array of unscaled band 3 TOA reflectance values */
    int16 *b4,     /* I: array of unscaled band 4 TOA reflectance values */
    int16 *b5,     /* I: array of unscaled band 5 TOA reflectance values */
    int16 *b6,     /* I: array of unscaled band 6 brightness temp values */
    int16 *b7,     /* I: array of unscaled band 7 TOA reflectance values */
    int nlines,    /* I: number of lines in the data arrays */
    int nsamps,    /* I: number of samples in the data arrays */
    float refl_scale_fact,  /* I: scale factor for the TOA reflectance values */
    float btemp_scale_fact, /* I: scale factor for the broghtness temp values */
    int refl_sat_value,     /* I: saturation value for TOA reflectance values */
    uint8 *refl_qa_mask, /* I: array of masked values for processing (non-zero
                               values are not to be processed) reflectance
                               bands */
    uint8 *snow_mask,    /* O: array of snow cover masked values (non-zero
                               values represent snow) */
    uint8 *probability_score /* O: probability pixel was classified correctly;
                               this is stored as a percentage between 0-100% */
);

void post_process_snow_cover_class
(
    int nlines,              /* I: number of lines in the data arrays */
    int nsamps,              /* I: number of samples in the data arrays */
    uint8 *snow_mask,        /* I/O: array of snow cover masked values
                                (non-zero values represent snow) */
    uint8 *probability_score /* I/O: probability pixel was classified correctly;
                                stored as a percentage between 0-100% */
);

float hillshade
(
    int16 *elev_window,   /* I: 3x3 array of elevation values in meters */
    float ew_res,         /* I: east/west resolution of the elevation data in
                                meters */
    float ns_res,         /* I: north/south resolution of the elevation data in
                                meters */
    float sun_elev,       /* I: sun elevation angle in radians */
    float solar_azimuth   /* I: solar azimuth angle in radians */
);

void deep_shadow
(
    int16 *dem,          /* I: array of DEM values in meters (nlines+[1or2] x
                               nsamps values - see NOTES);  if processing
                               at the top of the image, then an extra line
                               before will not be available;  if processing
                               at the bottom of the image, then an extra line
                               at the end will not be available */
    bool top,            /* I: are we at the top of the dem and therefore no
                               extra lines at the start of the dem? */
    bool bottom,         /* I: are we at the bottom of the dem and therefore no
                               extra lines at the end of the dem? */
    int nlines,          /* I: number of lines of data to be processed in the
                               mask array; dem array will have one or two lines
                               more depending on top, middle, bottom */
    int nsamps,          /* I: number of samples of data to be processed in the
                               mask array; dem array will have the same number
                               of samples therefore the first and last sample
                               will not be processed as part of the mask since
                               a 3x3 window won't be available */
    float ew_res,        /* I: east/west resolution of the elevation data in
                               meters */
    float ns_res,        /* I: north/south resolution of the elevation data in
                               meters */
    float sun_elev,      /* I: sun elevation angle in radians */
    float solar_azimuth, /* I: solar azimuth angle in radians */
    uint8 *shaded_relief,    /* O: array of shaded relief values (multiplied
                                   by 255 to take advantage of the 8-bit int)
                                   of size nlines * nsamps */
    uint8 *deep_shadow_mask  /* O: array of deep shadow masked values (non-zero
                                   values represent terrain-derived deep
                                   shadow areas) of size nlines * nsamps */
);
void refl_mask
(
    int16 *b1,     /* I: array of unscaled band 1 TOA reflectance values */
    int16 *b2,     /* I: array of unscaled band 2 TOA reflectance values */
    int16 *b3,     /* I: array of unscaled band 3 TOA reflectance values */
    int16 *b4,     /* I: array of unscaled band 4 TOA reflectance values */
    int16 *b5,     /* I: array of unscaled band 5 TOA reflectance values */
    int16 *b7,     /* I: array of unscaled band 7 TOA reflectance values */
    int nlines,    /* I: number of lines in the data arrays */
    int nsamps,    /* I: number of samples in the data arrays */
    int fill_value,   /* I: fill value for the TOA reflectance values */
    uint8 *refl_qa_mask  /* O: array of masked values for processing (non-zero
                               values are not to be processed) reflectance
                               bands */
);

void btemp_mask
(
    int16 *b6,     /* I: array of unscaled brightness temperature values */
    int nlines,    /* I: number of lines in the data arrays */
    int nsamps,    /* I: number of samples in the data arrays */
    int fill_value,   /* I: fill value for the brightness temp values */
    uint8 *btemp_qa_mask  /* O: array of masked values for processing (non-zero
                                values are not to be processed) brightness
                                temp bands */
);

#endif
