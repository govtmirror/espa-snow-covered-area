#include "sca.h"

/******************************************************************************
MODULE:  refl_mask

PURPOSE:  Generates a QA mask for no data pixels (fill) for the reflectance
bands.  If a pixel is "no data" in any of the bands, then the pixel will be
marked as no data in the QA mask.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/28/2012  Gail Schmidt     Original Development

NOTES:
******************************************************************************/
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
)
{
    int pix;          /* current pixel being processed */

    /* Loop through the pixels in the array to determine the quality mask */
    for (pix = 0; pix < nlines*nsamps; pix++)
    {
        /* If the current pixel in any band is fill then set the QA value to
           not be processed. Otherwise it's valid data. */
        if (b1[pix] == fill_value || b2[pix] == fill_value ||
            b3[pix] == fill_value || b4[pix] == fill_value ||
            b5[pix] == fill_value || b7[pix] == fill_value)
        {
            refl_qa_mask[pix] = NO_DATA;
        }
        else
            refl_qa_mask[pix] = VALID_DATA;
    }
}


/******************************************************************************
MODULE:  btemp_mask

PURPOSE:  Generates a QA mask for no data pixels (fill) for the brightness
temperature band.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/28/2012  Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void btemp_mask
(
    int16 *b6,     /* I: array of unscaled brightness temperature values */
    int nlines,    /* I: number of lines in the data arrays */
    int nsamps,    /* I: number of samples in the data arrays */
    int fill_value,   /* I: fill value for the brightness temp values */
    uint8 *btemp_qa_mask  /* O: array of masked values for processing (non-zero
                                values are not to be processed) brightness
                                temp bands */
)
{
    int pix;          /* current pixel being processed */

    /* Loop through the pixels in the array to determine the quality mask */
    for (pix = 0; pix < nlines*nsamps; pix++)
    {
        /* If the current pixel is fill then set the QA value to not be
           processed.  Otherwise it's valid data. */
        if (b6[pix] == fill_value)
            btemp_qa_mask[pix] = NO_DATA;
        else
            btemp_qa_mask[pix] = VALID_DATA;
    }
}

