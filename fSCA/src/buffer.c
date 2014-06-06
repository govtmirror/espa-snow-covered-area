#include <stdbool.h>
#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  buffer

PURPOSE:  Buffers the non-zero pixels by the specified distance.  Thus any
pixel which is distance pixels away from the nearest non-zero pixel is marked
with the same value as the non-zero pixel.  The distance is calculated at right
angles and not diagonally.

RETURN VALUE:
Type = bool
Value          Description
-----          -----------
ERROR          Error occurred buffering the array
SUCCESS        Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
---------     ---------------  -------------------------------------
6/4/2014      Gail Schmidt     Original Development

NOTES:
  1. Input and output arrays are 1D arrays of size nlines * nsamps.
******************************************************************************/
short buffer
(
    uint8 *array,       /* I: input array of data for which to buffer by the
                              distance value */
    int distance,       /* I: distance to buffer */
    int nlines,         /* I: number of lines in the data arrays */
    int nsamps,         /* I: number of samples in the data arrays */
    uint8 *buff_array   /* O: output array with buffer applied */
)
{
    char FUNC_NAME[] = "buffer";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int line;           /* current line being processed */
    int samp;           /* current sample being processed */
    int line_filter;    /* current line in the filter being processed */
    int samp_filter;    /* current sample in the filter being processed */
    int filter_size;    /* size of the filter based on the distance specified
                           for buffering */
    int fill_loc;       /* sample location for the current line to be filled
                           in the filter */
    int nfill;          /* number of samples to the left/right of the fill
                           location to be filled in the filter */
    long pix;           /* current pixel being processed */
    long filter_pix;    /* current pixel in the filter being processed */
    long buff_pix;      /* current pixel in the output buffered array */
    int HALF_WINDOW;    /* half the filter window size */
    uint8 *filter=NULL; /* filter to be applied; this will be a 2D filter of
                           size (distance * 2 + 1) */

    /* Determine the filter size.  Allocate memory for the filter and set the
       filter to zeros. */
    filter_size = 2 * distance + 1;
    filter = calloc (filter_size * filter_size, sizeof (uint8));
    if (filter == NULL)
    {
        strcpy (errmsg, "Error allocating memory for buffering filter.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Set up the filter to set the pixels which are distance from the center
       to one.  Example filter for distance of 3:
       0 0 0 1 0 0 0
       0 0 1 1 1 0 0
       0 1 1 1 1 1 0
       1 1 1 x 1 1 1
       0 1 1 1 1 1 0
       0 0 1 1 1 0 0
       0 0 0 1 0 0 0 */
    for (line_filter = 0; line_filter < filter_size; line_filter++)
    {
        /* Determine how many pixels from the center sample will be set/filled
           as one in the filter */
        fill_loc = abs (distance - line_filter);
        nfill = distance-fill_loc;
        filter_pix = line_filter * filter_size + (distance - nfill);
        for (samp_filter = distance-nfill; samp_filter <= distance+nfill;
             samp_filter++, filter_pix++)
        {
            filter[filter_pix] = 1;
        }
    }

#ifdef DEBUG
    printf ("DEBUG: Filter -->\n");
    for (line_filter = 0; line_filter < filter_size; line_filter++)
    {
        filter_pix = line_filter * filter_size;
        for (samp_filter = 0; samp_filter < filter_size; samp_filter++,
             filter_pix++)
            printf ("%d ", filter[filter_pix]);
        printf ("\n");
    }
#endif

    /* Set up the filter window parameters for processing the image */
    HALF_WINDOW = (int) (filter_size / 2);

    /* Initialize the buffered window to all zeros */
    memset (buff_array, 0, nlines * nsamps * sizeof (uint8));

    /* Loop through the lines and samples in the array and compute the
       variance.  Start at the HALF_WINDOW location, as no variance values
       will be calculated for windows in the fill areas. */
    for (line = 0; line < nlines; line++)
    {
        pix = line * nsamps;
        for (samp = 0; samp < nsamps; samp++, pix++)
        {
            /* If this pixel is zero then skip to the next pixel */
            if (array[pix] == 0)
                continue;

            /* Loop through the lines and samps in the window for this pixel
               and apply the buffer/filter */
            for (line_filter = -HALF_WINDOW; line_filter <= HALF_WINDOW;
                 line_filter++)
            {
                /* Skip this line in the filter window if it is outside the
                   image */
                if (line + line_filter < 0 || line + line_filter >= nlines)
                    continue;

                filter_pix = (line_filter + HALF_WINDOW) * filter_size;
                for (samp_filter = -HALF_WINDOW; samp_filter <= HALF_WINDOW;
                     samp_filter++, filter_pix++)
                {
                    /* Skip this sample in the filter window if it is outside
                       the image */
                    if (samp + samp_filter < 0 || samp + samp_filter >= nsamps)
                        continue;

                    /* If the pixel in the array is non-zero and the filter
                       for this pixel is set, then mark this pixel as a
                       buffered pixel */
                    if (filter[filter_pix] == 1)
                    {
                        buff_pix = (line + line_filter) * nsamps +
                            (samp + samp_filter);
                        buff_array[buff_pix] = array[pix];
                    }
                }  /* for samp_filter */
            }  /* for line_filter */
        }  /* for samp */
    }  /* for line */

    /* Free the memory */
    free (filter);

    /* Successful completion of the buffer */
    return (SUCCESS);
}
