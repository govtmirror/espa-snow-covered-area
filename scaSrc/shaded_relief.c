#include "sca.h"

/******************************************************************************
MODULE:  hillshade

PURPOSE:  Performs the hillshade algorithm (from GDALDEM) to compute the
shaded relief for the current 3x3 window.

RETURN VALUE:
Type = float 
Value      Description
-----      -----------
0.0 - 1.0  Represents the shaded relief for the current pixel area

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/31/2012  Gail Schmidt     Original Development (based on GDALHillshade
                             algorithm in GDALDEM v1.9.2)

NOTES:
  1. Algorithm is based on Lambert's cosine law using Horn's algorithm for
     calculating the slope of the current point.  The other option is to use
     Zevenbergen and Thorn's algorithm.  The litterature suggests Zevenbergen
     and Thorne to be more suited to smooth landscapes, whereas Horn's formula
     to perform better on rougher terrain.
  2. Input and output arrays are 1D arrays of size 3 lines x 3 samples where
     the 3x3 indices are and the current pixel being processed is 4.
       0 1 2
       3 4 5
       6 7 8
******************************************************************************/
float hillshade
(
    int16 *elev_window,   /* I: 3x3 array of elevation values in meters */
    float ew_res,         /* I: east/west resolution of the elevation data in
                                meters */
    float ns_res,         /* I: north/south resolution of the elevation data in
                                meters */
    float sun_elev,       /* I: sun elevation angle in radians; 0 deg at the
                                horizon and 90 deg if directly above the DEM */
    float solar_azimuth   /* I: solar azimuth angle in radians; 0 deg=North,
                                90 deg=East, 180 deg=South, 270 deg=West */
)
{
    float x_slope;        /* slope at this point in east/west direction */
    float y_slope;        /* slope at this point in north/south direction */
    float xx_plus_yy;     /* value of x * x + y * y */
    float aspect;         /* aspect at this point in radians */
    float relief;         /* shaded relief value at this point */

    /* Compute the slope */
    x_slope = ((elev_window[0] + 2.0 * elev_window[3] + elev_window[6]) -
               (elev_window[2] + 2.0 * elev_window[5] + elev_window[8])) /
               ew_res;
    y_slope = ((elev_window[6] + 2.0 * elev_window[7] + elev_window[8]) -
               (elev_window[0] + 2.0 * elev_window[1] + elev_window[2])) /
               ew_res;
    xx_plus_yy = x_slope * x_slope + y_slope * y_slope;

    /* Compute the aspect */
    aspect = atan2 (y_slope, x_slope);

    /* Compute the shade value */
    relief = (sin (sun_elev) - cos (sun_elev) * sqrt (xx_plus_yy) *
        sin (aspect - solar_azimuth)) / sqrt (1.0 + xx_plus_yy);

    return relief;
}


/******************************************************************************
MODULE:  deep_shadow_mask

PURPOSE:  Computes the shaded relief based on the DEM, then masks terrain-based
deep shadow pixels.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/31/2012  Gail Schmidt     Original Development

NOTES:
  1. Algorithm is based on the terrain-derived deep shadow algorithm provided
     by Dave Selkowitz, Research Geographer, USGS Alaska Science Center.
  2. Input DEM arrays are 1D arrays of size nlines+[1or2] * nsamps.  An extra
     line should be provided before and after (if possible) the actual subset
     of data to be processed.  Since we are processing an entire line of data
     at a time, it will not be expected that there will be an extra sample on
     either end of the line.  The first and last sample simply will not be
     processed for the shaded relief and the deep shadow mask, thus they should
     already be initialized to an appropriate value before calling this
     function.  The hillshade module requires a 3x3 window to be provided
     surrounding each pixel, thus the extra line(s) of data.
  3. Output mask arrays are 1D arrays of size nlines * nsamps.
******************************************************************************/
void deep_shadow
(
    int16 *dem,          /* I: array of DEM values in meters (nlines+[1or2] x
                               nsamps values - see NOTES);  if processing
                               at the top of the image, then an extra line
                               before will not be available;  if processing
                               at the bottom of the image, then an extra line
                               at the end will not be available */
    bool dem_top,        /* I: are we at the top of the dem and therefore no
                               extra lines at the start of the dem? */
    bool dem_bottom,     /* I: are we at the bottom of the dem and therefore no
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
)
{
    int line, samp;        /* line and sample being processed */
    int curr_pix;          /* current input pixel being processed */
    int out_pix;           /* current output pixel being processed */
    int start_line;        /* which line to start processing of the output
                              shaded relief and mask */
    int dem_line;          /* line in the DEM associated with the start_line
                              in the mask; should always be a value of 1 */
    int proc_nlines;       /* number of lines to process in the output mask */
    float shade;           /* shaded relief value for the current pixel */
    int16 elev_window[9];  /* 3x3 window of elevation values surrounding the
                              current pixel */

    /* Loop through the lines samples in the array to calculate the relief
       shading and determine the terrain-derived deep shadow mask.  The first
       line and column in the input array of DEM data is padding for the 3x3
       window, so they won't get processed.  However if we are at the top of
       the DEM then we don't have padding to do a 3x3 window, so just start at
       line 1.  At the end of the DEM, the same applies, so don't process the
       last line. */
    start_line = 0;
    proc_nlines = nlines;
    if (dem_top)
        start_line = 1;
    if (dem_bottom)
        proc_nlines = nlines-1;
    for (line = start_line, dem_line = 1; line < proc_nlines;
         line++, dem_line++)
    {
        for (samp = 1; samp < nsamps-1; samp++)
        {
            /* Determine the output pixel location for the current line/samp */
            out_pix = line * nsamps + samp;

            /* Fill in the 3x3 elevation window surrounding the current pixel */
            curr_pix = (dem_line-1) * nsamps + (samp-1);
            elev_window[0] = dem[curr_pix];
            elev_window[1] = dem[curr_pix+1];
            elev_window[2] = dem[curr_pix+2];
            curr_pix = dem_line * nsamps + (samp-1);
            elev_window[3] = dem[curr_pix];
            elev_window[4] = dem[curr_pix+1];
            elev_window[5] = dem[curr_pix+2];
            curr_pix = (dem_line+1) * nsamps + (samp-1);
            elev_window[6] = dem[curr_pix];
            elev_window[7] = dem[curr_pix+1];
            elev_window[8] = dem[curr_pix+2];

            /* Compute the shaded relief for the current pixel */
            shade = hillshade (elev_window, ew_res, ns_res, sun_elev,
                solar_azimuth);

            /* If the shaded relief value is below the shaded relief
               threshold, then mask this pixel as a terrain-derived deep
               shadow pixel */
            deep_shadow_mask[out_pix] = NO_DEEP_SHADOW;
            if (shade <= TERRAIN_DEEP_SHADOW_THRESH)
                deep_shadow_mask[out_pix] = DEEP_SHADOW;

            /* Scale the shaded relief values from 0.0 to 1.0 to 1-255, with
               0 representing DEM fill pixels */
            if (shade <= 0.0)
                shaded_relief[out_pix] = 1;
            else
                shaded_relief[out_pix] = 1 + (int) (254.0 * shade + 0.5);
        }
    }
}

