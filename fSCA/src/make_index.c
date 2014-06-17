#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  make_index

PURPOSE:  Computes the spectral index using the specified input bands.
index = (band1 - band2) / (band1 + band2)

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
---------     ---------------  -------------------------------------
5/19/2014     Gail Schmidt     Original Development
6/17/2014     Gail Schmidt     Don't handle the saturated values as special
                               cases

NOTES:
  1. Input and output arrays are 1D arrays of size nlines * nsamps.
  2. The index products will be created using the scaled reflectance
     values as it doesn't matter if they are scaled or unscaled for these
     simple band ratios.  Both bands are scaled by the same amount.
  3. If the current pixel is saturated in either band, then the output pixel
     value for the index will also be saturated.  The same applies for fill.
******************************************************************************/
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
)
{
    int pix;                /* current pixel being processed */
    float ratio;            /* band ratio */

    /* Loop through the pixels in the array and compute the spectral index */
    for (pix = 0; pix < nlines * nsamps; pix++)
    {
        /* If the current pixel is saturated in either band then the output
           is saturated.  Ditto for fill. */
        if (band1[pix] == fill_value || band2[pix] == fill_value)
            spec_indx[pix] = (float) FILL_VALUE;
        else
        {
            /* Compute the band ratio */
            ratio = (float) (band1[pix] - band2[pix]) /
                    (float) (band1[pix] + band2[pix]);

            /* Keep the ratio between -1.0, 1.0 */
            if (ratio > 1.0)
                ratio = 1.0;
            else if (ratio < -1.0)
                ratio = -1.0;
            spec_indx[pix] = ratio;
        }
    }
}

