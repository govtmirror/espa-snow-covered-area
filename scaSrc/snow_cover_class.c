#include "sca.h"

/******************************************************************************
MODULE:  snow_cover_class

PURPOSE:  Performs snow cover classification on TOA reflectance products

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/28/2012  Gail Schmidt     Original Development
1/7/2013    Gail Schmidt     Updated to add the water test and thermal test
                             as a post-processing step
1/9/2013    Gail Schmidt     If pixel is saturated then set it to the maximum
                             TOA reflectance value.
                             Added an additional thermal test for more general
                             false-positives.

NOTES:
  1. Algorithm is based on the snow cover classification algorithm provided by
     Dave Selkowitz, Research Geographer, USGS Alaska Science Center.
  2. Input and output arrays are 1D arrays of size nlines * nsamps.
******************************************************************************/
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
    float btemp_scale_fact, /* I: scale factor for the brightness temp values */
    int refl_sat_value,     /* I: saturation value for TOA reflectance values */
    uint8 *refl_qa_mask, /* I: array of masked values for processing (non-zero
                               values are not to be processed) reflectance
                               bands */
    uint8 *snow_mask,    /* O: array of snow cover masked values (non-zero
                               values represent snow) */
    uint8 *probability_score /* O: probability pixel was classified correctly;
                               this is stored as a percentage between 0-100% */
)
{
    uint8 sc_mask;    /* snow cover mask for the current pixel */
    uint8 prob_score; /* probability score (percentage) for current pixel */
    float ndsi;       /* normalized differenced snow index */
    float ndvi;       /* normalized differenced vegetation index */
    int pix;          /* current pixel being processed */
    float b1_pix;     /* scaled band 1 value for current pixel */
    float b2_pix;     /* scaled band 2 value for current pixel */
    float b3_pix;     /* scaled band 3 value for current pixel */
    float b4_pix;     /* scaled band 4 value for current pixel */
    float b5_pix;     /* scaled band 5 value for current pixel */
    float b6_pix;     /* scaled band 6 value for current pixel */
    float b7_pix;     /* scaled band 7 value for current pixel */

    /* Loop through the pixels in the array to determine the snow cover
       classification */
    for (pix = 0; pix < nlines*nsamps; pix++)
    {
        /* If the reflective QA mask is turned on, then skip the snow cover
           processing for this pixel */
        if (refl_qa_mask[pix] != 0)
        {
            snow_mask[pix] = NO_SNOW;
            probability_score[pix] = 0;
            continue;
        }

        /* Scale the current pixel for each band.  If the pixel is saturated,
           then set it to it's maximum instead of the saturated value.  The
           maximum TOA reflectance value is 1.0.  Given that we are focused
           on cold pixels (clouds, snow), we won't worry about saturated
           thermal pixels and will just use the thermal values as-is. */
        b6_pix = b6[pix] * btemp_scale_fact;

        if (b1[pix] == refl_sat_value)
            b1_pix = 1.0;
        else
            b1_pix = b1[pix] * refl_scale_fact;

        if (b1[pix] == refl_sat_value)
            b2_pix = 1.0;
        else
            b2_pix = b2[pix] * refl_scale_fact;

        if (b1[pix] == refl_sat_value)
            b3_pix = 1.0;
        else
            b3_pix = b3[pix] * refl_scale_fact;

        if (b1[pix] == refl_sat_value)
            b4_pix = 1.0;
        else
            b4_pix = b4[pix] * refl_scale_fact;

        if (b1[pix] == refl_sat_value)
            b5_pix = 1.0;
        else
            b5_pix = b5[pix] * refl_scale_fact;

        if (b1[pix] == refl_sat_value)
            b7_pix = 1.0;
        else
            b7_pix = b7[pix] * refl_scale_fact;

        /* Run the short-circuit tests for water and false-positives first
           to save time from running the other tests */
        /* Water test - if this is water then it's not snow */
        if (b4_pix < 0.11)
        {
            snow_mask[pix] = NO_SNOW;
            probability_score[pix] = 3;
            continue;
        }

        /* Thermal test - to catch barren and other non-snow areas that
           appear as snow */
        if (b6_pix > 24.85) /* 298 K */
        {
            snow_mask[pix] = NO_SNOW;
            probability_score[pix] = 3;
            continue;
        }

        /* If band 2 and band 5 are zero, then the NDSI cannot be computed */
        if (b2_pix == 0.0 && b5_pix == 0.0)
        {
            snow_mask[pix] = NO_SNOW;
            probability_score[pix] = 3;
            continue;
        }

        /* Compute the NDSI */
        ndsi = (b2_pix - b5_pix) / (b2_pix + b5_pix);

        /* If band 3 and band 4 are zero, then the NDVI cannot be computed */
        if (b3_pix == 0.0 && b4_pix == 0.0)
        {
            snow_mask[pix] = NO_SNOW;
            probability_score[pix] = 3;
            continue;
        }

        /* Compute the NDVI */
        ndvi = (b4_pix - b3_pix) / (b4_pix + b3_pix);

        /* Initialize the pixel to no snow cover */
        sc_mask = NO_SNOW;
        prob_score = 3;

        /* Determine snow cover via the binary snow cover tree */
        if (ndsi < 0.25)
        {
            if (b5_pix >= 0.072)
            {
                if (b3_pix < 0.35)
                {
                    sc_mask = NO_SNOW;
                    prob_score = 98;
                }
                else
                {
                    sc_mask = SNOW_COVER;
                    prob_score = 100;
                }
            }
            else
            {
                sc_mask = SNOW_COVER;
                prob_score = 90;
            }
        }
        else
        {
            if (b7_pix >= 0.14)
            {
                if (b4_pix < 0.32)
                {
                    if (ndvi < 0.19)
                    {
                        if (b1_pix < 0.3)
                        {
                            sc_mask = NO_SNOW;
                            prob_score = 96;
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 100;
                        }
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 91;
                    }
                }
                else
                {
                    if (ndsi < 0.57)
                    {
                        if (ndvi < 0.18)
                        {
                            if (b1_pix < 0.35)
                            {
                                sc_mask = NO_SNOW;
                                prob_score = 79;
                            }
                            else
                            {
                                if (b7_pix >= 0.22)
                                {
                                    sc_mask = NO_SNOW;
                                    prob_score = 84;
                                }
                                else
                                {
                                    sc_mask = SNOW_COVER;
                                    prob_score = 88;
                                }
                            }
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 95;
                        }
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 93;
                    }
                }
            }
            else
            {
                if (b5_pix >= 0.1)
                {
                    if (b7_pix >= 0.11)
                    {
                        if (b1_pix < 0.34)
                        {
                            sc_mask = NO_SNOW;
                            prob_score = 90;
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 100;
                        }
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 81;
                    }
                }
                else
                {
                    sc_mask = SNOW_COVER;
                    prob_score = 99;
                }
            }
        }

        /* Post-processing thermal test - to catch forested areas which are
           false positives */
        if (b6_pix > 18.85 /* 292 K */ && prob_score < 98)
        {
            sc_mask = NO_SNOW;
            prob_score = 2;
        }

        snow_mask[pix] = sc_mask;
        probability_score[pix] = prob_score;
    }  /* end for pix */
}
