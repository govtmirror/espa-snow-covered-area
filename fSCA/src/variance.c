#include <math.h>
#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  variance

PURPOSE:  Computes the variance of the 9x9 window of pixels in the input array.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
---------     ---------------  -------------------------------------
5/19/2014     Gail Schmidt     Original Development

NOTES:
  1. Input and output arrays are 1D arrays of size nlines * nsamps.
  2. Any window containing a fill value will not have a variance calculated.
******************************************************************************/
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
)
{
    bool skip_window;   /* should the current window be skipped due to fill */
    int i;              /* looping variable */
    int line;           /* current line being processed */
    int samp;           /* current sample being processed */
    int line_window;    /* current line in the window being processed */
    int samp_window;    /* current sample in the window being processed */
    long pix;           /* current pixel being processed */
    long win_pix;       /* current pixel in the window being processed */
    int count;          /* count of the pixels in the window */
    int WINDOW = 9;     /* variance window 9x9 around the center pixel */
    int HALF_WINDOW = (int) (WINDOW / 2); /* half window size */
    int SQR_WINDOW = (WINDOW * WINDOW); /* squared window size */
    float out_scale;    /* output scale factor to be applied to convert to
                           an integer for writing */
    double diff;        /* difference between the current pixel and the
                           average */
    double avg;         /* average/mean value of the window */
    double sum;         /* sum of values in the window */
    double var;         /* variance of values in the window */
    double unscaled_arr[SQR_WINDOW];  /* array of unscaled values */
    double epsilon = 0.00001;  /* value for comparing the floating points */

    /* Fill the variance array will fill values by default, since any window
       with a fill pixel will not have a variance calculated */
    for (pix = 0; pix < nlines * nsamps; pix++)
        variance[pix] = fill_value;

    /* If the input scale factor is 1.0 then don't apply an output scale
       to convert to an integer */
    if (fabs (scale_factor - 1.0) < epsilon)
        out_scale = 1.0;
    else
        out_scale = FLOAT_TO_INT;

    /* Loop through the lines and samples in the array and compute the
       variance.  Start at the HALF_WINDOW location, as no variance values
       will be calculated for windows in the fill areas. */
    for (line = HALF_WINDOW; line < nlines-HALF_WINDOW; line++)
    {
        pix = line * nsamps + HALF_WINDOW;
        for (samp = HALF_WINDOW; samp < nsamps-HALF_WINDOW; samp++, pix++)
        {
            /* Loop through the lines and samps in the window for this pixel
               and compute the sum */
            skip_window = false;
            sum = 0.0;
            count = 0;
            for (line_window = -HALF_WINDOW; line_window <= HALF_WINDOW;
                 line_window++)
            {
                win_pix = (line + line_window) * nsamps + (samp - HALF_WINDOW);
                for (samp_window = -HALF_WINDOW; samp_window <= HALF_WINDOW;
                     samp_window++, win_pix++, count++)
                {
                    /* If this pixel is fill break out, otherwise add it to
                       the window average */
                    if (array[win_pix] == fill_value)
                    {
                        skip_window = true;
                        break;
                    }
                    else
                    {
                        unscaled_arr[count] = array[win_pix] * scale_factor;
                        sum += unscaled_arr[count];
                    }
                }  /* for samp_window */

                /* Should the window be skipped */
                if (skip_window)
                    break;
            }  /* for line_window */

            /* Should the pixel be skipped */
            if (skip_window)
                continue;

            /* Compute the average and then loop back through computing the
               variance.  No need to check for fill values, since if there
               are any the code already jumped to the next pixel.  Use the
               unscaled values stored in the unscaled array. */
            avg = sum / SQR_WINDOW;
            sum = 0.0;
            for (i = 0; i < SQR_WINDOW; i++)
            {
                diff = unscaled_arr[i] - avg;
                sum += diff * diff;
            }

            /* Assign the variance to the current pixel */
            var = sum / (SQR_WINDOW - 1);

            /* Scale to an int32 */
            if (var >= 0.0)
                variance[pix] = (int32) (var * out_scale + 0.5);
            else
                variance[pix] = (int32) (var * out_scale - 0.5);
        }  /* for samp */
    }  /* for line */
}

