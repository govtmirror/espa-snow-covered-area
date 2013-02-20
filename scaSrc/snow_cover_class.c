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
2/6/2013    Gail Schmidt     Added subclassifications for nodes 3 and 15 to
                             work with false positives in heavy conifer areas

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
    uint8 *probability_score,/* O: probability pixel was classified correctly;
                               this is stored as a percentage between 0-100%
                               (used for debugging) */
    uint8 *tree_node,    /* O: node in tree used to classify each pixel */
    uint8 *ndsi_array,   /* O: NDSI outputs */
    uint8 *ndvi_array    /* O: NDVI outputs (used for debugging) */
)
{
    uint8 sc_mask;    /* snow cover mask for the current pixel */
    uint8 prob_score; /* probability score (percentage) for current pixel */
    uint8 node;       /* tree node for current pixel */
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
            tree_node[pix] = 0;
            ndsi_array[pix] = 0;
            ndvi_array[pix] = 0;
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
        snow_mask[pix] = NO_SNOW;
        probability_score[pix] = 3;
        tree_node[pix] = 0;
        ndsi_array[pix] = 0;
        ndvi_array[pix] = 0;

        /* Water test - if this is water then it's not snow */
        if (b4_pix < 0.11)
            continue;

        /* Thermal test - to catch barren and other non-snow areas that
           appear as snow */
        if (b6_pix > 24.85) /* 298 K */
            continue;

        /* If band 2 and band 5 are zero, then the NDSI cannot be computed */
        if (b2_pix == 0.0 && b5_pix == 0.0)
            continue;

        /* Compute the NDSI */
        ndsi = (b2_pix - b5_pix) / (b2_pix + b5_pix);

        /* If band 3 and band 4 are zero, then the NDVI cannot be computed */
        if (b3_pix == 0.0 && b4_pix == 0.0)
            continue;

        /* Compute the NDVI */
        ndvi = (b4_pix - b3_pix) / (b4_pix + b3_pix);

        /* Initialize the pixel to no snow cover */
        sc_mask = NO_SNOW;
        prob_score = 3;
        node = 0;
        if (ndsi < 0.0)
            ndsi_array[pix] = 0.0;
        else
            ndsi_array[pix] = ndsi * 100.0 + 0.5;
        if (ndvi < 0.0)
            ndvi_array[pix] = 0.0;
        else
            ndvi_array[pix] = ndvi * 100.0 + 0.5;

        /* Determine snow cover via the binary snow cover tree */
        if (ndsi < 0.25)
        {
            if (b5_pix >= 0.072)
            {
                if (b3_pix < 0.35)
                {
                    sc_mask = NO_SNOW;
                    prob_score = 98;
                    node = 1;
                }
                else
                {
                    sc_mask = SNOW_COVER;
                    prob_score = 100;
                    node = 2;
                }
            }
            else
            {
                if (b1_pix < 0.11)
                {
                    if (ndsi < 0.21)
                    {
                        sc_mask = NO_SNOW;
                        prob_score = 98;
                        node = 31;
                    }
                    else
                    {
                        if (b3_pix < 0.042)
                        {
                            sc_mask = NO_SNOW;
                            prob_score = 91;
                            node = 32;
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 85;
                            node = 33;
                        }
                    }
                }
                else
                {
                    if (ndsi < 0.15)
                    {
                        sc_mask = NO_SNOW;
                        prob_score = 83;
                        node = 34;
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 95;
                        node = 35;
                    }
                }
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
                            node = 4;
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 100;
                            node = 5;
                        }
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 91;
                        node = 6;
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
                                node = 7;
                            }
                            else
                            {
                                if (b7_pix >= 0.22)
                                {
                                    sc_mask = NO_SNOW;
                                    prob_score = 84;
                                    node = 8;
                                }
                                else
                                {
                                    sc_mask = SNOW_COVER;
                                    prob_score = 88;
                                    node = 9;
                                }
                            }
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 95;
                            node = 10;
                        }
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 93;
                        node = 11;
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
                            node = 12;
                        }
                        else
                        {
                            sc_mask = SNOW_COVER;
                            prob_score = 100;
                            node = 13;
                        }
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 81;
                        node = 14;
                    }
                }
                else
                {
                    if (b3_pix < 0.0432)
                    {
                        sc_mask = NO_SNOW;
                        prob_score = 93;
                        node = 151;
                    }
                    else
                    {
                        sc_mask = SNOW_COVER;
                        prob_score = 99;
                        node = 152;
                    }
                }
            }
        }

        /* Post-processing thermal test - to catch forested areas which are
           false positives */
        if (b6_pix > 18.85 /* 292 K */ && prob_score < 98)
        {
            sc_mask = NO_SNOW;
            prob_score = 2;
            node = 0;
        }

        /* Assign the snow cover and probability scores to the current pixel */
        snow_mask[pix] = sc_mask;
        probability_score[pix] = prob_score;
        tree_node[pix] = node;
    }  /* end for pix */
}


/******************************************************************************
MODULE:  post_process_snow_cover_class

PURPOSE:  Performs snow cover classification post-processing to clear out
    false positives from the snow cover algorithm, particularly in the dense
    conifer regions flagged by nodes 3 or 15 of the snow cover binary tree.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/4/2013    Gail Schmidt     Original Development

NOTES:
  1. Algorithm is based on the snow cover classification algorithm provided by
     Dave Selkowitz, Research Geographer, USGS Alaska Science Center.
  2. If a pixel is identified as snow covered from nodes 3 or 15 from the
     original snow cover binary tree (i.e. node pixels 3-1, 3-2 ... 15-1,
     15-2 from the extended binary tree), check the other pixels in a 9x9
     window to count how many of those pixels have been identified as snow
     covered.  If there are fewer than the snow cover threshold, then change
     the pixel to be snow free.
  3. Input and output arrays are 1D arrays of size nlines * nsamps.
******************************************************************************/
void post_process_snow_cover_class
(
    int nlines,         /* I: number of lines in the data arrays */
    int nsamps,         /* I: number of samples in the data arrays */
    uint8 *snow_mask,   /* I/O: array of snow cover masked values (non-zero
                                values represent snow) */
    uint8 *tree_node    /* I: node in binary tree used to classify each pixel */
)
{
    int count;              /* number of snow-covered pixels in NxN window */
    int line, samp;         /* current line and sample being processed */
    int pix;                /* current pixel being processed */
    int start_window_line;  /* starting line for the NxN window */
    int end_window_line;    /* ending line for the NxN window */
    int start_window_samp;  /* starting sample for the NxN window */
    int end_window_samp;    /* ending sample for the NxN window */
    int win_line, win_samp; /* current line and sample being processed in the
                               NxN window */
    int win_pix;            /* current window pixel being processed */
    int orig_node;          /* nodes 3 and 15 have been expanded and therefore
                               are of a value of 3-x (31, 33 .. 35) and 15-x
                               (151, 152).  this variable will represent the
                               original node number without the secondary
                               node numbers */
    static int HALF_WINDOW = 4;  /* half of 9x9 window around current pixel */
    static float SNOW_COUNT_THRESH = 7; /* threshold for count of pixels in
                                           the NxN window needing to be snow */

    /* Loop through the pixels in the array to determine the snow cover
       classification */
    for (line = 0; line < nlines; line++)
    {
        /* Find the valid NxN window for the current line */
        start_window_line = line - HALF_WINDOW;
        end_window_line = line + HALF_WINDOW;
        if (start_window_line < 0)
            start_window_line = 0;
        if (end_window_line >= nlines)
            end_window_line = nlines - 1;

        for (samp = 0; samp < nsamps; samp++)
        {
            /* Calculate the location of the current pixel in the 1D array */
            pix = line * nsamps + samp;

            /* If the current pixel is snow covered and the tree node is 15-x
               or 3-x then count the snow-covered pixels in the surrounding NxN
               window.  If that count exceeds the threshold, then leave the
               pixel as snow.  Otherwise change the mask to not snow covered. */
            orig_node = tree_node[pix] / 10;
            if (snow_mask[pix] == SNOW_COVER && (orig_node == 15 ||
                orig_node == 3))
            {
                /* Find the valid NxN window for the current line */
                start_window_samp = samp - HALF_WINDOW;
                end_window_samp = samp + HALF_WINDOW;
                if (start_window_samp < 0)
                    start_window_samp = 0;
                if (end_window_samp >= nsamps)
                    end_window_samp = nsamps - 1;

                /* Loop through the NxN window (or whatever smaller window is
                   available) */
                count = 0;
                for (win_line = start_window_line; win_line <= end_window_line;
                     win_line++)
                {
                    /* Calculate the starting location of the current window
                       pixel in the 1D array */
                    win_pix = win_line * nsamps + start_window_samp;
                    for (win_samp = start_window_samp;
                         win_samp <= end_window_samp; win_samp++, win_pix++)
                    {
                        if (snow_mask[win_pix] == SNOW_COVER)
                            count++;
                    }
                }

                /* If the snow cover count for this pixel is less than the
                   threshold, then reset it to no snow cover */
                if (count < SNOW_COUNT_THRESH)
                {
                    snow_mask[pix] = NO_SNOW;
                }
            }  /* end if snow cover and nodes 3 or 15 */
        }  /* end for samp */
    }  /* end for line */
}


/******************************************************************************
MODULE:  count_adjacent_snow_cover

PURPOSE:  Performs an assessment of the snow cover results by counting the
    3x3 window around the current pixel and to determine if 1) there are
    any adjacent cloud or fill values, and 2) count the number of adjacent
    snow cover pixels.  If the current pixel has adjacent cloud, deep shadow,
    or fill pixels, then flag that pixel as such.  Otherwise output the count
    of adjacent snow pixels.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/13/2013    Gail Schmidt     Original Development

NOTES:
  1. Algorithm is based on the snow cover classification algorithm provided by
     Dave Selkowitz, Research Geographer, USGS Alaska Science Center.
  2. Input and output arrays are 1D arrays of size nlines * nsamps.
  3. Non-zero values represent snow, cloud, fill in the input masks.
******************************************************************************/
void count_adjacent_snow_cover
(
    int nlines,           /* I: number of lines in the data arrays */
    int nsamps,           /* I: number of samples in the data arrays */
    uint8 *snow_mask,     /* I: array of snow cover masked values */
    uint8 *cloud_mask,    /* I: array of cloud masked values */
    uint8 *shadow_mask,   /* I: array of deep shadow masked values */
    uint8 *refl_qa_mask,  /* I: quality mask for the TOA reflectance products,
                                where fill and saturated values are flagged */
    uint8 *btemp_qa_mask, /* I: quality mask for the brightness temp products,
                                where fill and saturated values are flagged */
    uint8 *snow_count     /* O: count of the snow cover results in the adjacent
                                3x3 window, or high value if one or more of
                                the adjacent pixels are cloud/shadow/fill */

)
{
    int count;              /* number of snow-covered pixels in NxN window */
    int line, samp;         /* current line and sample being processed */
    int pix;                /* current pixel being processed */
    int start_window_line;  /* starting line for the NxN window */
    int end_window_line;    /* ending line for the NxN window */
    int start_window_samp;  /* starting sample for the NxN window */
    int end_window_samp;    /* ending sample for the NxN window */
    int win_line, win_samp; /* current line and sample being processed in the
                               NxN window */
    int win_pix;            /* current window pixel being processed */
    static int HALF_WINDOW = 1;  /* half of 3x3 window around current pixel */

    /* Loop through the pixels in the array to count the adjacent snow cover
       pixels */
    for (line = 0; line < nlines; line++)
    {
        /* Find the valid NxN window for the current line */
        start_window_line = line - HALF_WINDOW;
        end_window_line = line + HALF_WINDOW;
        if (start_window_line < 0)
            start_window_line = 0;
        if (end_window_line >= nlines)
            end_window_line = nlines - 1;

        for (samp = 0; samp < nsamps; samp++)
        {
            /* Calculate the location of the current pixel in the 1D array */
            pix = line * nsamps + samp;

            /* Find the valid NxN window for the current line */
            start_window_samp = samp - HALF_WINDOW;
            end_window_samp = samp + HALF_WINDOW;
            if (start_window_samp < 0)
                start_window_samp = 0;
            if (end_window_samp >= nsamps)
                end_window_samp = nsamps - 1;

            /* Loop through the NxN window (or whatever smaller window is
               available) */
            count = 0;
            for (win_line = start_window_line; win_line <= end_window_line;
                 win_line++)
            {
                /* Calculate the starting location of the current window
                   pixel in the 1D array */
                win_pix = win_line * nsamps + start_window_samp;
                for (win_samp = start_window_samp;
                     win_samp <= end_window_samp; win_samp++, win_pix++)
                {
                    /* If the current pixel is cloud, shadow, or fill, then
                       flag it and move to the next pixel in the image.
                       Otherwise keep track of the count of adjacent snow
                       pixels. */
                    if (cloud_mask[win_pix] == CLOUD_COVER ||
                        refl_qa_mask[win_pix] == NO_DATA ||
                        btemp_qa_mask[win_pix] == NO_DATA ||
                        shadow_mask[win_pix] == DEEP_SHADOW)
                    {
                        snow_count[pix] = ADJ_PIX_MASKED;
                        win_line = end_window_line+1;  /* make sure we break
                            the line loop for the window as well */
                        break;
                    }

                    /* Otherwise keep track of the count of adjacent snow
                       pixels */
                    if (snow_mask[win_pix] == SNOW_COVER)
                        count++;
                }
            }

            /* If the current pixel doesn't have adjacent cloud/shadow/fill
               pixels, then assign the count of adjacent snow pixels */
            if (snow_count[pix] != ADJ_PIX_MASKED)
                snow_count[pix] = count;
        }  /* end for samp */
    }  /* end for line */
}
