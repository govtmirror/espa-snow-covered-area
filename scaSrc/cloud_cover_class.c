#include "sca.h"

/******************************************************************************

MODULE:  cloud_cover_class

PURPOSE:  Performs cloud cover classification on TOA reflectance products

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/28/2012  Gail Schmidt     Original Development
1/8/2013    Gail Schmidt     Converted the brightness temp constants from
                             degrees Kelvin to degrees Celsius since the
                             LEDAPS brightness temps are in Celsius

NOTES:
  1. Algorithm is based on the cloud cover classification tree provided by
     Dave Selkowitz, Research Geographer, USGS Alaska Science Center.
  2. Input and output arrays are 1D arrays of size nlines * nsamps.
******************************************************************************/
void cloud_cover_class
(
    int16 *b1,     /* I: array of unscaled band 1 TOA reflectance values */
    int16 *b4,     /* I: array of unscaled band 4 TOA reflectance values */
    int16 *b6,     /* I: array of unscaled band 6 brightness temp values,
                         in degrees Celsius */
    int16 *b7,     /* I: array of unscaled band 7 TOA reflectance values */
    int nlines,    /* I: number of lines in the data arrays */
    int nsamps,    /* I: number of samples in the data arrays */
    float refl_scale_fact,  /* I: scale factor for the TOA reflectance values */
    float btemp_scale_fact, /* I: scale factor for the brightness temp values */
    uint8 *refl_qa_mask,  /* I: array of masked values for processing (non-zero
                                values are not to be processed) reflectance
                                bands */
    uint8 *therm_qa_mask, /* I: array of masked values for processing (non-zero
                                values are not to be processed) thermal bands */
    uint8 *cloud_mask     /* O: array of cloud cover masked values (non-zero
                                values represent clouds) */
)
{
    uint8 cc_mask;    /* cloud cover mask for the current pixel */
    int pix;          /* current pixel being processed */
    float b1_pix;     /* scaled band 1 value for current pixel */
    float b4_pix;     /* scaled band 4 value for current pixel */
    float b6_pix;     /* scaled band 6 value for current pixel */
    float b7_pix;     /* scaled band 7 value for current pixel */

    /* Loop through the pixels in the array to determine the cloud cover
       classification */
    for (pix = 0; pix < nlines*nsamps; pix++)
    {
        /* Scale the current pixel for each band */
        b1_pix = b1[pix] * refl_scale_fact;
        b4_pix = b4[pix] * refl_scale_fact;
        b6_pix = b6[pix] * btemp_scale_fact;
        b7_pix = b7[pix] * refl_scale_fact;

        /* Initialize the pixel to no cloud cover */
        cc_mask = NO_CLOUD;

        /* If the thermal or reflective QA mask is turned on, then skip the
           cloud cover processing for this pixel */
        if (refl_qa_mask[pix] != 0 || therm_qa_mask[pix] != 0)
        {
            cloud_mask[pix] = cc_mask;
            continue;
        }

        /* Determine cloud cover */
        if (b1_pix < 0.30095)
        {
            if (b1_pix < 0.20055)
                cc_mask = NO_CLOUD;
            else
            {
                if (b7_pix < 0.08255)
                    cc_mask = NO_CLOUD;
                else
                {
                    if (b6_pix < -7.052)   /* 266.098 K */
                        cc_mask = CLOUD_COVER;
                    else
                        cc_mask = NO_CLOUD;
                }
            }
        }
        else
        {
            if (b7_pix < 0.1166)
            {
                if (b6_pix < -19.316)   /* 253.834 K */
                    cc_mask = CLOUD_COVER;
                else
                    cc_mask = NO_CLOUD;
            }
            else
            {
                if (b7_pix < 0.15305)
                {
                    if (b6_pix < -20.036)    /* 253.114 K */
                        cc_mask = CLOUD_COVER;
                    else
                        cc_mask = NO_CLOUD;
                }
                else
                {
                    if (b6_pix < 8.788)    /* 281.938 K */
                    {
                        if (b4_pix < 1.04525)
                            cc_mask = CLOUD_COVER;
                        else
                            cc_mask = NO_CLOUD;
                    }
                    else
                        cc_mask = NO_CLOUD;
                }
            }
        }

        cloud_mask[pix] = cc_mask;
    }  /* end for pix */
}
