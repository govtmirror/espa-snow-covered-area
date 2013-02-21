#include "sca.h"

/******************************************************************************
MODULE:  combine_qa_mask

PURPOSE:  Combines the QA masks from cloud, deep shadow, and fill into one
overall mask.  If the current pixel is flagged as any of these, then the
combined QA mask is turned on.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/21/2013    Gail Schmidt     Original Development

NOTES:
  1. Input and output arrays are 1D arrays of size nlines * nsamps.
  2. Non-zero values represent snow, cloud, fill in the input masks.
******************************************************************************/
void combine_qa_mask
(
    int nlines,           /* I: number of lines in the data arrays */
    int nsamps,           /* I: number of samples in the data arrays */
    uint8 *cloud_mask,    /* I: array of cloud masked values */
    uint8 *shadow_mask,   /* I: array of deep shadow masked values */
    uint8 *refl_qa_mask,  /* I: quality mask for the TOA reflectance products,
                                where fill and saturated values are flagged */
    uint8 *btemp_qa_mask, /* I: quality mask for the brightness temp products,
                                where fill and saturated values are flagged */
    uint8 *combined_qa    /* O: combined mask representing cloud, deep shadow,
                                and fill for the current pixel */
)
{
    int pix;                /* current pixel being processed */

    /* Loop through the pixels in the array to count the adjacent snow cover
       pixels */
    for (pix = 0; pix < nlines * nsamps; pix++)
    {
        /* If the current pixel is cloud, deep shadow, or fill, then flag it
           in the combined mask. */
        if (cloud_mask[pix] == CLOUD_COVER || refl_qa_mask[pix] == NO_DATA ||
            btemp_qa_mask[pix] == NO_DATA || shadow_mask[pix] == DEEP_SHADOW)
        {
            combined_qa[pix] = COMBINED_MASK;
        }
    }
}
