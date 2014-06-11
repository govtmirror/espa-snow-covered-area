#include <math.h>
#include <string.h>
#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  rule_based_model

PURPOSE:  Runs the rule-based models on the input line of data, including the
reflectance bands, NDVI, NDSI, and the variance for each of those products.
Five model runs are made per model, and there are two models.  The first model
uses the variances and the second model does not use the variances.

RETURN VALUE:
Type = none

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
---------     ---------------  -------------------------------------
5/21/2014     Gail Schmidt     Original Development

NOTES:
  1. Input and output arrays are 1D arrays of size NPROC_LINES * nsamps.
  2. This algorithm was provided by David Selkowitz, USGS Alaska Science Center.
  3. The algorithm uses the scaled reflectance values as-is.
  4. The algorithm will unscale the NDSI, NDVI, and index variance values, as
     they were scaled before writing to the output file as integers and they
     need to be unscaled before being used by the model.  The SCALE_FACTOR in
     output.h will be applied to unscale these values.
******************************************************************************/
void rule_based_model
(
    Input_t *input_img,     /* I: pointer to input data structure containing
                                  the scaled reflectance and cloud mask
                                  buffers (reflectance bands are scaled) */
    int16 *ndsi_arr,        /* I: NDSI scaled values (scaled) */
    int16 *ndvi_arr,        /* I: NDVI scaled values (scaled) */
    int32 *b1_var_arr,      /* I: band1 variance values */
    int32 *b2_var_arr,      /* I: band2 variance values */
    int32 *b4_var_arr,      /* I: band4 variance values */
    int32 *b5_var_arr,      /* I: band5 variance values */
    int32 *b7_var_arr,      /* I: band7 variance values */
    int32 *ndvi_var_arr,    /* I: NDVI variance values (scaled) */
    int32 *ndsi_var_arr,    /* I: NDSI variance values (scaled) */
    int nsamps,             /* I: number of samples in the input arrays */
    uint8 *rev_cloud_mask,      /* O: revised cloud mask */
    uint8 *rev_lim_cloud_mask   /* O: revised cloud mask without variances */
)
{
    int samp;             /* current sample being processed */
    int16 m1_cloud_code;  /* cloud mask for the 1st Rule-based model */
    int16 m2_cloud_code;  /* cloud mask for the 2nd Rule-based model */
    int16 m3_cloud_code;  /* cloud mask for the 3rd Rule-based model */
    int16 m4_cloud_code;  /* cloud mask for the 4th Rule-based model */
    int16 m5_cloud_code;  /* cloud mask for the 5th Rule-based model */
    int16 m1_lim_cloud_code;  /* limited cloud mask for 1st Rule-based model */
    int16 m2_lim_cloud_code;  /* limited cloud mask for 2nd Rule-based model */
    int16 m3_lim_cloud_code;  /* limited cloud mask for 3rd Rule-based model */
    int16 m4_lim_cloud_code;  /* limited cloud mask for 4th Rule-based model */
    int16 m5_lim_cloud_code;  /* limited cloud mask for 5th Rule-based model */
    int16 conserv_cloud_code; /* maximum of the model runs for the
                                 conservative cloud code */
    int16 limited_cloud_code; /* maximum of the model runs for the
                                 limited (non-variance) cloud code */
    int16 b1, b2, b3, b4, b5, b7;  /* pixel values for the reflectance bands */
    int32 b1_var, b2_var, b4_var, b5_var, b7_var; /* variances for the
                                      reflectance bands */
    double ndvi, ndsi;    /* pixel values for the NDVI and NDSI bands */
    double ndvi_var, ndsi_var;  /* pixel values for NDVI and NDSI variances */

    /* Initialize the revised cloud mask to all zeros */
    memset (rev_cloud_mask, 0, nsamps * sizeof (uint8));
    memset (rev_lim_cloud_mask, 0, nsamps * sizeof (uint8));

    /* Loop through the pixels in the array and run the Rule-based models.
       Any pixel which is not flagged as cloudy in the input cfmask will be
       skipped. */
    for (samp = 0; samp < nsamps; samp++)
    {
        /* If this isn't a cloudy pixel in the cfmask then skip to the next
           pixel */
        if (input_img->cfmask_buf[samp] != 4)
            continue;

        /* Initialize the cloud mask */
        m1_cloud_code = -1;
        m2_cloud_code = -1;
        m3_cloud_code = -1;
        m4_cloud_code = -1;
        m5_cloud_code = -1;

        /* Set up the individual pixel values. Don't unscale the reflectance
           values or the reflectance variance values, but do unscale the
           NDVI and NDSI values and variances. */
        b1 = input_img->refl_buf[0][samp];
        b2 = input_img->refl_buf[1][samp];
        b3 = input_img->refl_buf[2][samp];
        b4 = input_img->refl_buf[3][samp];
        b5 = input_img->refl_buf[4][samp];
        b7 = input_img->refl_buf[5][samp];
        b1_var = b1_var_arr[samp];
        b2_var = b2_var_arr[samp];
        b4_var = b4_var_arr[samp];
        b5_var = b5_var_arr[samp];
        b7_var = b7_var_arr[samp];
        ndvi = ndvi_arr[samp] * SCALE_FACTOR;
        ndsi = ndsi_arr[samp] * SCALE_FACTOR;
        ndvi_var = ndvi_var_arr[samp] * SCALE_FACTOR;
        ndsi_var = ndsi_var_arr[samp] * SCALE_FACTOR;
        
        /*** Conservative cloud mask model run using the variances ***/
        /* 1st model run */
        /* Rule 0/1: (8646/45, lift 2.1) */
        if (m1_cloud_code == -1 && b7 <= 2077 && b1_var > 308729 && 
            b7_var <= 234934 && ndsi_var > 0.00215485)
            m1_cloud_code = 50;  /* class cloud_free  [0.995] */
        
        /* Rule 0/2: (270/1, lift 2.1) */
        if (m1_cloud_code == -1 && b1 <= 2839 && b3 > 2917 && 
            ndsi_var <= 0.0020372)
            m1_cloud_code = 50;  /* class cloud_free  [0.993] */
        
        /* Rule 0/3: (10298/67, lift 2.1) */
        if (m1_cloud_code == -1 && b7 <= 2077 && ndsi_var > 0.04491)
            m1_cloud_code = 50;  /* class cloud_free  [0.993] */
        
        /* Rule 0/4: (11273/132, lift 2.1) */
        if (m1_cloud_code == -1 && b3 > 1341 && b7 <= 1464 && 
            ndsi_var > 0.00816928)
            m1_cloud_code = 50;  /* class cloud_free  [0.988] */
        
        /* Rule 0/5: (979/11, lift 2.1) */
        if (m1_cloud_code == -1 && b1 <= 2999 && ndsi <= -0.199329 && 
            b1_var <= 86087.4 && ndvi_var > 0.00170841)
            m1_cloud_code = 50;  /* class cloud_free  [0.988] */
        
        /* Rule 0/6: (9044/143, lift 2.1) */
        if (m1_cloud_code == -1 && b5 <= 1005)
            m1_cloud_code = 50;  /* class cloud_free  [0.984] */
        
        /* Rule 0/7: (1837/110, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 2999 && b7 > 1464 && 
            ndvi <= 0.0817003 && ndsi <= -0.0544693 && ndsi_var > 0.000305468)
            m1_cloud_code = 50;  /* class cloud_free  [0.940] */
        
        /* Rule 0/8: (961/69, lift 1.9) */
        if (m1_cloud_code == -1 && b1 <= 2999 && b3 > 2917)
            m1_cloud_code = 50;  /* class cloud_free  [0.927] */
        
        /* Rule 0/9: (342/28, lift 1.9) */
        if (m1_cloud_code == -1 && ndsi <= -0.199329 && b4_var > 5655.39 && 
            b7_var <= 10260.5)
            m1_cloud_code = 50;  /* class cloud_free  [0.916] */
        
        /* Rule 0/10: (2161/200, lift 1.9) */
        if (m1_cloud_code == -1 && b1 <= 2999 && ndvi <= 0.0817003 && 
            ndsi <= -0.0544693)
            m1_cloud_code = 50;  /* class cloud_free  [0.907] */
        
        /* Rule 0/11: (23277/7218, lift 1.4) */
        if (m1_cloud_code == -1 && b1 <= 2999)
            m1_cloud_code = 50;  /* class cloud_free  [0.690] */
        
        /* Rule 0/12: (13760/38, lift 1.9) */
        if (m1_cloud_code == -1 && b1 > 2999 && b7 > 1464 && b1_var <= 308729)
            m1_cloud_code = 100;  /* class cloud  [0.997] */
        
        /* Rule 0/13: (11071/44, lift 1.9) */
        if (m1_cloud_code == -1 && b1 > 2839 && b7 > 1464 && 
            ndsi_var <= 0.0020372)
            m1_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 0/14: (2176/14, lift 1.9) */
        if (m1_cloud_code == -1 && b7 > 1464 && ndvi > 0.0817003 && 
            b1_var > 86087.4 && ndsi_var <= 0.0020372)
            m1_cloud_code = 100;  /* class cloud  [0.993] */
        
        /* Rule 0/15: (1041/9, lift 1.9) */
        if (m1_cloud_code == -1 && b3 <= 2917 && b7 > 1464 && 
            ndvi > 0.0817003 && b4_var <= 5655.39)
            m1_cloud_code = 100;  /* class cloud  [0.990] */
        
        /* Rule 0/16: (10040/128, lift 1.9) */
        if (m1_cloud_code == -1 && b7 > 1464 && ndsi_var <= 0.000305468)
            m1_cloud_code = 100;  /* class cloud  [0.987] */
        
        /* Rule 0/17: (3796/101, lift 1.9) */
        if (m1_cloud_code == -1 && b3 <= 2917 && b7 > 1464 && 
            ndvi > 0.0817003 && ndsi > -0.199329 && ndsi_var <= 0.0020372)
            m1_cloud_code = 100;  /* class cloud  [0.973] */
        
        /* Rule 0/18: (1151/34, lift 1.9) */
        if (m1_cloud_code == -1 && b7 > 1464 && ndvi > 0.169039 && 
            b7_var > 10260.5 && ndvi_var <= 0.00170841 && ndsi_var <= 0.0020372)
            m1_cloud_code = 100;  /* class cloud  [0.970] */
        
        /* Rule 0/19: (18644/561, lift 1.9) */
        if (m1_cloud_code == -1 && b1 > 2999 && b7 > 1464)
            m1_cloud_code = 100;  /* class cloud  [0.970] */
        
        /* Rule 0/20: (13648/467, lift 1.8) */
        if (m1_cloud_code == -1 && b5 > 1005 && ndvi <= 0.410316 && 
            ndsi > -0.117693 && ndsi_var <= 0.0023633)
            m1_cloud_code = 100;  /* class cloud  [0.966] */
        
        /* Rule 0/21: (2291/87, lift 1.8) */
        if (m1_cloud_code == -1 && b7 > 1464 && b7_var > 312503 && 
            ndsi_var <= 0.0109275)
            m1_cloud_code = 100;  /* class cloud  [0.962] */
        
        /* Rule 0/22: (5902/238, lift 1.8) */
        if (m1_cloud_code == -1 && b7 > 1464 && ndvi > 0.0410272 && 
            ndsi > -0.128104 && ndvi_var <= 0.00150992 && ndsi_var <= 0.0109275)
            m1_cloud_code = 100;  /* class cloud  [0.960] */
        
        /* Rule 0/23: (289/15, lift 1.8) */
        if (m1_cloud_code == -1 && b1 > 1639 && b5 > 1005 && b7 <= 1464 && 
            b7_var > 103775 && ndsi_var <= 0.00816928)
            m1_cloud_code = 100;  /* class cloud  [0.945] */
        
        /* Rule 0/24: (10865/622, lift 1.8) */
        if (m1_cloud_code == -1 && b5 > 1005 && ndvi_var <= 0.00017896 && 
            ndsi_var <= 0.0023633)
            m1_cloud_code = 100;  /* class cloud  [0.943] */
        
        /* Rule 0/25: (597/42, lift 1.8) */
        if (m1_cloud_code == -1 && b5 > 1005 && b7 <= 1464 && 
            ndvi <= 0.410316 && b2_var > 14609.2 && ndsi_var <= 0.0023633)
            m1_cloud_code = 100;  /* class cloud  [0.928] */
        
        /* Rule 0/26: (231/21, lift 1.7) */
        if (m1_cloud_code == -1 && b1 > 1536 && b3 <= 1341 && 
            b7_var > 84988.3 && ndvi_var <= 0.00255233 && ndsi_var <= 0.0347677)
            m1_cloud_code = 100;  /* class cloud  [0.906] */
        
        /* Rule 0/27: (942/93, lift 1.7) */
        if (m1_cloud_code == -1 && b5 > 1005 && b5 <= 1538 && 
            ndvi <= 0.410316 && b2_var <= 57818 && ndsi_var <= 0.00416476)
            m1_cloud_code = 100;  /* class cloud  [0.900] */
        
        if (m1_cloud_code == -1)
            m1_cloud_code = 100;


        /* 2nd model run */
        /* Rule 1/1: (2863.5/6.1, lift 2.1) */
        if (m2_cloud_code == -1 && ndsi > 0.863486)
            m2_cloud_code = 50;  /* class cloud_free  [0.998] */
        
        /* Rule 1/2: (524.1/13.6, lift 2.1) */
        if (m2_cloud_code == -1 && b2 <= 3971 && ndsi <= -0.272461 && 
            b7_var > 217782)
            m2_cloud_code = 50;  /* class cloud_free  [0.972] */
        
        /* Rule 1/3: (3592/133.3, lift 2.0) */
        if (m2_cloud_code == -1 && b2 <= 3971 && ndsi > 0.403054)
            m2_cloud_code = 50;  /* class cloud_free  [0.963] */
        
        /* Rule 1/4: (1431.9/108.2, lift 2.0) */
        if (m2_cloud_code == -1 && b1 <= 3503 && b3 > 3076 && 
            b7_var <= 217782 && ndsi_var > 0.000593467)
            m2_cloud_code = 50;  /* class cloud_free  [0.924] */
        
        /* Rule 1/5: (412.1/35.6, lift 1.9) */
        if (m2_cloud_code == -1 && b1 <= 2683 && b3 > 2573 && 
            ndsi_var <= 0.000593467)
            m2_cloud_code = 50;  /* class cloud_free  [0.912] */
        
        /* Rule 1/6: (3358/315.7, lift 1.9) */
        if (m2_cloud_code == -1 && b2 > 3971  && b5 <= 2549 && b4_var > 267311)
            m2_cloud_code = 50;  /* class cloud_free  [0.906] */
        
        /* Rule 1/7: (2239.2/352, lift 1.8) */
        if (m2_cloud_code == -1 && b1 <= 3503  && b7 > 2282 && 
            b7_var <= 217782 && ndsi_var > 0.000593467 && ndsi_var <= 0.0345117)
            m2_cloud_code = 50;  /* class cloud_free  [0.842] */
        
        /* Rule 1/8: (34996.9/15368.2, lift 1.2) */
        if (m2_cloud_code == -1 && b2 <= 3971)
            m2_cloud_code = 50;  /* class cloud_free  [0.561] */
        
        /* Rule 1/9: (9359.6/39.8, lift 1.9) */
        if (m2_cloud_code == -1 && b2 > 3971 && b5 > 2549)
            m2_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 1/10: (192.7, lift 1.9) */
        if (m2_cloud_code == -1 && ndsi <= 0.403054 && b4_var <= 3374.14 && 
            ndsi_var > 0.000593467)
            m2_cloud_code = 100;  /* class cloud  [0.995] */
        
        /* Rule 1/11: (8302/47.1, lift 1.9) */
        if (m2_cloud_code == -1 && b2 > 3971 && ndsi <= 0.863486 && 
            b4_var <= 267311)
            m2_cloud_code = 100;  /* class cloud  [0.994] */
        
        /* Rule 1/12: (4769.8/48.3, lift 1.9) */
        if (m2_cloud_code == -1 && b1 > 3503 && ndsi <= 0.403054 && 
            b4_var > 3374.14 && b7_var <= 217782 && ndsi_var <= 0.0345117)
            m2_cloud_code = 100;  /* class cloud  [0.990] */
        
        /* Rule 1/13: (833.5/44.3, lift 1.8) */
        if (m2_cloud_code == -1 && b1 > 1825  && b2 <= 1540 && 
            ndsi <= 0.403054 && ndsi_var > 0.000593467 && ndsi_var <= 0.0345117)
            m2_cloud_code = 100;  /* class cloud  [0.946] */
        
        /* Rule 1/14: (11430/1036, lift 1.7) */
        if (m2_cloud_code == -1 && b5 > 1090 && ndsi <= 0.863486 && 
            b7_var > 0 && ndsi_var <= 0.000593467)
            m2_cloud_code = 100;  /* class cloud  [0.909] */
        
        /* Rule 1/15: (3499.4/422.9, lift 1.7) */
        if (m2_cloud_code == -1 && b1 > 1825 && b7 > 1167 && b7 <= 2282 && 
            ndsi <= 0.403054 && b1_var > 12517 && ndsi_var <= 0.00544464)
            m2_cloud_code = 100;  /* class cloud  [0.879] */
        
        /* Rule 1/16: (4696.5/581.1, lift 1.7) */
        if (m2_cloud_code == -1 && ndvi <= 0.392095 && ndsi > -0.272461 && 
            ndsi <= 0.863486 && b7_var > 217782 && ndvi_var > 0.000185316 && 
            ndsi_var <= 0.0345117)
            m2_cloud_code = 100;  /* class cloud  [0.876] */
        
        /* Rule 1/17: (5900.9/751, lift 1.6) */
        if (m2_cloud_code == -1 && ndvi <= 0.392095  && ndsi <= 0.863486 && 
            b5_var > 904312)
            m2_cloud_code = 100;  /* class cloud  [0.873] */
        
        /* Rule 1/18: (1387.9/244, lift 1.6) */
        if (m2_cloud_code == -1 && b2 > 1011  && b3 <= 1099 && 
            ndvi <= 0.392095 && ndsi <= 0.403054 && b7_var <= 217782 && 
            ndsi_var <= 0.0345117)
            m2_cloud_code = 100;  /* class cloud  [0.824] */
        
        if (m2_cloud_code == -1)
            m2_cloud_code = 100;


        /* 3rd model run */
        /* Rule 2/1: (408.9, lift 2.1) */
        if (m3_cloud_code == -1 && b5 <= 723 && ndsi_var <= 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.998] */
        
        /* Rule 2/2: (2473.9/26.8, lift 2.1) */
        if (m3_cloud_code == -1 && b3 <= 3829 && ndvi > -0.0104225 && 
            b1_var > 1920000 && b5_var <= 1610000 && b7_var <= 425578 && 
            ndvi_var <= 0.0652907 && ndsi_var > 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.989] */
        
        /* Rule 2/3: (1579.3/20.2, lift 2.1) */
        if (m3_cloud_code == -1 && b1 <= 2587 && ndvi > -0.0104225 && 
            b1_var > 28042.1 && b7_var <= 46730.4 && ndsi_var > 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.987] */
        
        /* Rule 2/4: (1662.5/27.7, lift 2.1) */
        if (m3_cloud_code == -1 && b2 > 1321  && b3 <= 3829 && 
            ndvi > -0.0104225 && ndsi > 0.267579 && ndvi_var <= 0.0652907 && 
            ndsi_var > 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.983] */
        
        /* Rule 2/5: (126.4/2.4, lift 2.0) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b1_var <= 175885 && 
            b2_var > 190105 && ndsi_var <= 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.974] */
        
        /* Rule 2/6: (1403.3/57.5, lift 2.0) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b3 > 2528 && b7_var <= 1860000)
            m3_cloud_code = 50;  /* class cloud_free  [0.958] */
        
        /* Rule 2/7: (1411.1/105.7, lift 1.9) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b7 > 2683 && 
            b1_var <= 175885 && ndsi_var > 0.000125769)
            m3_cloud_code = 50;  /* class cloud_free  [0.925] */
        
        /* Rule 2/8: (1074.5/99, lift 1.9) */
        if (m3_cloud_code == -1 && b1 <= 3218 && b2 > 3032 && b7_var <= 1860000)
            m3_cloud_code = 50;  /* class cloud_free  [0.907] */
        
        /* Rule 2/9: (5184/544, lift 1.9) */
        if (m3_cloud_code == -1 && b1 > 2598 && b7 <= 1208)
            m3_cloud_code = 50;  /* class cloud_free  [0.895] */
        
        /* Rule 2/10: (892.6/97.8, lift 1.9) */
        if (m3_cloud_code == -1 && b2 > 1321 && b3 <= 3829 && 
            b1_var <= 28042.1 && ndvi_var <= 0.0878528 && ndsi_var > 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.890] */
        
        /* Rule 2/11: (6939.7/1213, lift 1.7) */
        if (m3_cloud_code == -1 && b1 <= 2587 && b2 > 1321 && 
            b7_var <= 425578 && ndsi_var > 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.825] */
        
        /* Rule 2/12: (1701.5/298.6, lift 1.7) */
        if (m3_cloud_code == -1 && b3 > 3829 && b1_var > 0 && 
            b7_var <= 162839 && ndsi_var > 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.824] */
        
        /* Rule 2/13: (2671.6/480.1, lift 1.7) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b3 > 1051 && b5 > 1026 && 
            ndvi <= 0.0978049 && b1_var <= 175885 && b5_var > 15651.8 && 
            ndsi_var > 0.000191934 && ndsi_var <= 0.00365532)
            m3_cloud_code = 50;  /* class cloud_free  [0.820] */
        
        /* Rule 2/14: (42439.5/19761.5, lift 1.1) */
        if (m3_cloud_code == -1 && b1_var > 0)
            m3_cloud_code = 50;  /* class cloud_free  [0.534] */
        
        /* Rule 2/15: (1023.8/0.6, lift 1.9) */
        if (m3_cloud_code == -1 && b7_var > 1860000)
            m3_cloud_code = 100;  /* class cloud  [0.998] */
        
        /* Rule 2/16: (267.5, lift 1.9) */
        if (m3_cloud_code == -1 && b3 <= 2528 && b1_var > 0 && 
            ndvi_var <= 3.78e-05)
            m3_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 2/17: (229.7, lift 1.9) */
        if (m3_cloud_code == -1 && b3 <= 2528 && b5 > 723 && b1_var > 0 && 
            b5_var <= 1908.3)
            m3_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 2/18: (126.9, lift 1.9) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b2_var > 2109 && 
            b4_var <= 3585.67 && b5_var > 5639.05 && ndsi_var > 0.000191934 && 
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.992] */
        
        /* Rule 2/19: (2506.1/29.9, lift 1.9) */
        if (m3_cloud_code == -1 && b5 > 723 && ndsi <= 0.881556 && 
            ndvi_var > 3.78e-05 && ndsi_var <= 0.000125769)
            m3_cloud_code = 100;  /* class cloud  [0.988] */
        
        /* Rule 2/20: (1252.6/19.9, lift 1.9) */
        if (m3_cloud_code == -1 && b7 <= 2683 && ndsi <= 0.881556 && 
            b2_var > 2109 && ndsi_var <= 0.000191934)
            m3_cloud_code = 100;  /* class cloud  [0.983] */
        
        /* Rule 2/21: (6791.6/166, lift 1.9) */
        if (m3_cloud_code == -1 && b1 > 3218 && b7 > 1208 && 
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.975] */
        
        /* Rule 2/22: (744.5/20.4, lift 1.9) */
        if (m3_cloud_code == -1 && b5 > 1026 && b7 <= 2683 && 
            ndvi <= 0.0978049 && ndsi <= 0.881556 && b1_var <= 175885 && 
            b5_var > 5639.05 && b5_var <= 15651.8 && ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.971] */
        
        /* Rule 2/23: (5925.2/178.5, lift 1.9) */
        if (m3_cloud_code == -1 && ndsi <= 0.881556  && b1_var <= 0)
            m3_cloud_code = 100;  /* class cloud  [0.970] */
        
        /* Rule 2/24: (471.4/13.7, lift 1.9) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b3 <= 2528 && b5 > 723 && 
            b1_var > 175885 && ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.969] */
        
        /* Rule 2/25: (4937.7/251.3, lift 1.8) */
        if (m3_cloud_code == -1 && b1 > 1604 && b2 > 1321 && b5 > 899 && 
            ndvi <= -0.0104225 && ndsi <= 0.881556 && ndvi_var <= 0.0652907 && 
            ndsi_var <= 0.0378146)
            m3_cloud_code = 100;  /* class cloud  [0.949] */
        
        /* Rule 2/26: (2754.4/162.4, lift 1.8) */
        if (m3_cloud_code == -1 && b1 > 2598 && b2 <= 3032 && b7 > 1208 && 
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.941] */
        
        /* Rule 2/27: (491.2/28.8, lift 1.8) */
        if (m3_cloud_code == -1 && b5 > 899 && ndsi <= 0.881556 && 
            ndvi_var > 0.0652907 && ndsi_var > 0.00365532 && 
            ndsi_var <= 0.0378146)
            m3_cloud_code = 100;  /* class cloud  [0.940] */
        
        /* Rule 2/28: (351.6/21.4, lift 1.8) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b5 > 723 && b5 <= 1026 && 
            b1_var > 0 && b2_var <= 190105 && ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.937] */
        
        /* Rule 2/29: (2898.5/186.1, lift 1.8) */
        if (m3_cloud_code == -1 && ndvi > 0.0978049 && ndsi > -0.228882 && 
            ndsi <= 0.093273 && b2_var > 2109 && b2_var <= 190105 && 
            b5_var > 5639.05 && ndvi_var <= 0.00049353 &&
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.935] */
        
        /* Rule 2/30: (7312.9/623.2, lift 1.8) */
        if (m3_cloud_code == -1 && b1 > 2587 && ndsi <= 0.267579 && 
            b1_var <= 1920000 && b5_var <= 1610000 && ndvi_var <= 0.0652907 && 
            ndsi_var <= 0.0378146)
            m3_cloud_code = 100;  /* class cloud  [0.915] */
        
        /* Rule 2/31: (493.2/45.1, lift 1.7) */
        if (m3_cloud_code == -1 && b1 <= 2598 && b5 > 723 && b7 <= 752 && 
            ndvi <= 0.44504 && b1_var > 0 && b2_var <= 190105 && 
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.907] */
        
        /* Rule 2/32: (3829.7/353.5, lift 1.7) */
        if (m3_cloud_code == -1 && b3 > 3829 && ndsi <= 0.881556 && 
            b7_var > 162839 && ndsi_var <= 0.113723)
            m3_cloud_code = 100;  /* class cloud  [0.907] */
        
        /* Rule 2/33: (10939.1/1102, lift 1.7) */
        if (m3_cloud_code == -1 && b3 > 3829 && ndsi <= 0.881556 && 
            ndsi_var <= 0.113723)
            m3_cloud_code = 100;  /* class cloud  [0.899] */
        
        /* Rule 2/34: (738.8/80.5, lift 1.7) */
        if (m3_cloud_code == -1 && b1 > 1604 && b2 <= 1321 && 
            ndsi_var > 0.00365532 && ndsi_var <= 0.113723)
            m3_cloud_code = 100;  /* class cloud  [0.890] */
        
        /* Rule 2/35: (2149.8/246.3, lift 1.7) */
        if (m3_cloud_code == -1 && b1 > 1604 && b5 > 899 && 
            b5_var <= 1610000 && b7_var > 425578 && ndsi_var <= 0.0378146)
            m3_cloud_code = 100;  /* class cloud  [0.885] */
        
        /* Rule 2/36: (1826.2/221.4, lift 1.7) */
        if (m3_cloud_code == -1 && b5 > 3063 && b7 <= 2683 && 
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.878] */
        
        /* Rule 2/37: (5110.1/634.1, lift 1.7) */
        if (m3_cloud_code == -1 && b3 > 1051 && ndvi > 0.0978049 && 
            ndvi <= 0.44504 && ndsi > -0.228882 && ndsi <= 0.093273 && 
            b2_var > 2109 && b5_var > 5639.05 && b7_var > 9219.45 && 
            ndsi_var <= 0.00365532)
            m3_cloud_code = 100;  /* class cloud  [0.876] */
        
        /* Rule 2/38: (1808.8/399.7, lift 1.5) */
        if (m3_cloud_code == -1 && b1 > 1604 && b5 > 899 && 
            ndsi <= 0.881556 && ndvi_var > 0.0878528 && 
            ndsi_var > 0.00365532 && ndsi_var <= 0.113723)
            m3_cloud_code = 100;  /* class cloud  [0.779] */
        
        /* Rule 2/39: (3338.3/765.1, lift 1.5) */
        if (m3_cloud_code == -1 && b1 > 1885 && b5 > 899 && b5 <= 2807 && 
            b7 <= 2226 && b1_var > 28042.1 && b1_var <= 1920000 && 
            b5_var <= 1610000 && b7_var > 46730.4 && ndsi_var <= 0.0378146)
            m3_cloud_code = 100;  /* class cloud  [0.771] */
        
        /* Rule 2/40: (4731.9/1390.1, lift 1.4) */
        if (m3_cloud_code == -1 && b1 > 1604 && b5 <= 2807 && b7 <= 2226 && 
            b1_var > 28042.1 && b1_var <= 1920000 && b5_var <= 1610000 && 
            b7_var > 46730.4 && ndsi_var <= 0.0378146)
            m3_cloud_code = 100;  /* class cloud  [0.706] */
        
        if (m3_cloud_code == -1)
            m3_cloud_code = 100;


        /* 4th model run */
        /* Rule 3/1: (1616, lift 2.0) */
        if (m4_cloud_code == -1 && ndsi > 0.888811)
            m4_cloud_code = 50;  /* class cloud_free  [0.999] */
        
        /* Rule 3/2: (1937.6/10.3, lift 2.0) */
        if (m4_cloud_code == -1 && b5 <= 4153 && b4_var > 1090000 && 
            b7_var <= 325542)
            m4_cloud_code = 50;  /* class cloud_free  [0.994] */
        
        /* Rule 3/3: (2810.4/27.7, lift 2.0) */
        if (m4_cloud_code == -1 && b5 <= 4153 && ndvi > -0.00044238 && 
            b2_var > 89280.1 && b2_var <= 1.15e+07 && ndsi_var > 0.0355158)
            m4_cloud_code = 50;  /* class cloud_free  [0.990] */
        
        /* Rule 3/4: (509.9/5.1, lift 2.0) */
        if (m4_cloud_code == -1 && b1 <= 2874 && b3 > 2844 && b5 <= 4153 && 
            ndsi_var <= 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.988] */
        
        /* Rule 3/5: (255.6/2.3, lift 2.0) */
        if (m4_cloud_code == -1 && b5 <= 4153 && ndvi > -0.00044238 && 
            b2_var <= 89280.1 && b7_var > 325542)
            m4_cloud_code = 50;  /* class cloud_free  [0.987] */
        
        /* Rule 3/6: (1022.4/19.8, lift 2.0) */
        if (m4_cloud_code == -1 && b1 <= 2499 && b3 > 2484)
            m4_cloud_code = 50;  /* class cloud_free  [0.980] */
        
        /* Rule 3/7: (1331.4/35.4, lift 1.9) */
        if (m4_cloud_code == -1 && b1 <= 2874 && b1_var <= 103107 && 
            ndvi_var > 0.00504752)
            m4_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 3/8: (4721/161.2, lift 1.9) */
        if (m4_cloud_code == -1 && b5 <= 4153 && b5_var <= 1.05e+07 && 
            ndsi_var > 0.0780484)
            m4_cloud_code = 50;  /* class cloud_free  [0.966] */
        
        /* Rule 3/9: (224.5/6.9, lift 1.9) */
        if (m4_cloud_code == -1 && b5 <= 4153 && b1_var > 5.47e+07 && 
            b2_var <= 1.15e+07 && b7_var > 325542 && ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.965] */
        
        /* Rule 3/10: (1647.9/81.9, lift 1.9) */
        if (m4_cloud_code == -1 && b1 <= 1367 && b7_var <= 325542 && 
            ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.950] */
        
        /* Rule 3/11: (529.6/26.4, lift 1.9) */
        if (m4_cloud_code == -1 && b5 <= 1480 && b5_var <= 15695.8 && 
            ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.948] */
        
        /* Rule 3/12: (4600.6/271.2, lift 1.9) */
        if (m4_cloud_code == -1 && b2 > 1673 && b5 <= 1480 && 
            b7_var <= 325542 && ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.941] */
        
        /* Rule 3/13: (747.8/47.3, lift 1.9) */
        if (m4_cloud_code == -1 && b5 > 1480 && b5 <= 4153 && 
            b1_var > 3400000 && b2_var <= 274632 && ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.936] */
        
        /* Rule 3/14: (704.2/49.6, lift 1.9) */
        if (m4_cloud_code == -1 && b1 > 1367 && b2 <= 1673 && 
            ndvi <= 0.00974093 && ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.928] */
        
        /* Rule 3/15: (486.5/35.6, lift 1.9) */
        if (m4_cloud_code == -1 && b1 <= 2499 && ndvi <= -0.0123566 && 
            b1_var <= 103107)
            m4_cloud_code = 50;  /* class cloud_free  [0.925] */
        
        /* Rule 3/16: (3003.1/235.2, lift 1.8) */
        if (m4_cloud_code == -1 && b1 <= 2040 && b5 > 1480 && 
            ndvi <= 0.223642 && b7_var <= 325542 && ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.921] */
        
        /* Rule 3/17: (10526.7/874.7, lift 1.8) */
        if (m4_cloud_code == -1 && b5 <= 4153 && b5_var <= 1.05e+07 && 
            b7_var <= 325542 && ndsi_var > 0.0177706)
            m4_cloud_code = 50;  /* class cloud_free  [0.917] */
        
        /* Rule 3/18: (1419.2/201.5, lift 1.7) */
        if (m4_cloud_code == -1 && b1 <= 2499 && b7 > 2458 && b1_var <= 103107)
            m4_cloud_code = 50;  /* class cloud_free  [0.858] */
        
        /* Rule 3/19: (8128.2/1254.3, lift 1.7) */
        if (m4_cloud_code == -1 && b1 <= 3255  && b2 > 1763 && b5 <= 4153 && 
            ndvi <= 0.223642 && b5_var <= 1.05e+07 && b7_var <= 325542 && 
            ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.846] */
        
        /* Rule 3/20: (4365.8/698.2, lift 1.7) */
        if (m4_cloud_code == -1 && b1 <= 2884 && ndvi > -0.00044238 && 
            ndvi <= 0.0894992 && b2_var > 89280.1 && b5_var <= 1.05e+07 && 
            ndsi_var > 0.00144666)
            m4_cloud_code = 50;  /* class cloud_free  [0.840] */
        
        /* Rule 3/21: (1063.6/187.6, lift 1.6) */
        if (m4_cloud_code == -1 && ndvi > 0.392625)
            m4_cloud_code = 50;  /* class cloud_free  [0.823] */
        
        /* Rule 3/22: (897.1/272.9, lift 1.4) */
        if (m4_cloud_code == -1 && b1 <= 2499 && ndvi > 0.0394758 && 
            ndvi <= 0.107908 && b1_var <= 11846.8 && b4_var > 5354.76)
            m4_cloud_code = 50;  /* class cloud_free  [0.695] */
        
        /* Rule 3/23: (2089/669.5, lift 1.4) */
        if (m4_cloud_code == -1 && ndsi <= -0.201802 && b2_var <= 16676)
            m4_cloud_code = 50;  /* class cloud_free  [0.679] */
        
        /* Rule 3/24: (724.7/301.1, lift 1.2) */
        if (m4_cloud_code == -1 && b1 <= 3446 && b5 > 4153)
            m4_cloud_code = 50;  /* class cloud_free  [0.584] */
        
        /* Rule 3/25: (632.5/0.5, lift 2.0) */
        if (m4_cloud_code == -1 && b1 > 2884  && b5 <= 4153 && 
            b1_var <= 5.47e+07 && b2_var > 89280.1 && b5_var <= 1.05e+07 && 
            b7_var > 325542 && ndsi_var <= 0.0355158)
            m4_cloud_code = 100;  /* class cloud  [0.998] */
        
        /* Rule 3/26: (404.5, lift 2.0) */
        if (m4_cloud_code == -1 && b1 <= 2499 && b3 <= 2484 && 
            b1_var > 103107 && ndsi_var <= 0.00144666)
            m4_cloud_code = 100;  /* class cloud  [0.998] */
        
        /* Rule 3/27: (5612.6/10.7, lift 2.0) */
        if (m4_cloud_code == -1 && b1 > 3446 && b5 > 4153)
            m4_cloud_code = 100;  /* class cloud  [0.998] */
        
        /* Rule 3/28: (746.4/8.1, lift 2.0) */
        if (m4_cloud_code == -1 && b1 > 2040 && b2 <= 1763 && b5 > 1480 && 
            b2_var <= 274632)
            m4_cloud_code = 100;  /* class cloud  [0.988] */
        
        /* Rule 3/29: (282.9/10.7, lift 1.9) */
        if (m4_cloud_code == -1 && b3 <= 2844 && b7 <= 2458 && 
            ndvi <= 0.0394758 && b1_var <= 11846.8 && b4_var > 5354.76 && 
            ndsi_var <= 0.00144666)
            m4_cloud_code = 100;  /* class cloud  [0.959] */
        
        /* Rule 3/30: (1413.5/63.1, lift 1.9) */
        if (m4_cloud_code == -1 && b2_var > 1.15e+07 && b7_var > 325542 && 
            ndsi_var <= 0.0780484)
            m4_cloud_code = 100;  /* class cloud  [0.955] */
        
        /* Rule 3/31: (588.1/30.1, lift 1.9) */
        if (m4_cloud_code == -1 && b7 <= 2458 && ndvi <= 0.392625 && 
            ndsi <= -0.201802 && b2_var > 16676 && ndvi_var <= 0.00504752 && 
            ndsi_var <= 0.00144666)
            m4_cloud_code = 100;  /* class cloud  [0.947] */
        
        /* Rule 3/32: (2585.3/155.5, lift 1.9) */
        if (m4_cloud_code == -1 && b5_var > 1.05e+07)
            m4_cloud_code = 100;  /* class cloud  [0.940] */
        
        /* Rule 3/33: (938.6/67.7, lift 1.9) */
        if (m4_cloud_code == -1 && ndvi <= -0.00044238 && b1_var <= 5.47e+07 && 
            b2_var <= 1.15e+07 && b7_var > 325542 && ndsi_var <= 0.0780484)
            m4_cloud_code = 100;  /* class cloud  [0.927] */
        
        /* Rule 3/34: (1250.1/124, lift 1.8) */
        if (m4_cloud_code == -1 && ndvi > 0.0894992 && b2_var > 89280.1 && 
            b7_var > 325542 && ndsi_var <= 0.0355158)
            m4_cloud_code = 100;  /* class cloud  [0.900] */
        
        /* Rule 3/35: (5470.8/578.9, lift 1.8) */
        if (m4_cloud_code == -1 && ndvi > 0.107908 && ndvi <= 0.392625 && 
            ndsi > -0.201802 && ndvi_var <= 0.00504752 && 
            ndsi_var <= 0.00144666)
            m4_cloud_code = 100;  /* class cloud  [0.894] */
        
        /* Rule 3/36: (2035.7/366.4, lift 1.6) */
        if (m4_cloud_code == -1 && b1 > 1367 && b2 <= 1673 && b5 <= 1480 && 
            ndvi > 0.00974093 && b4_var <= 1090000 && b5_var > 15695.8 && 
            ndsi_var <= 0.0177706)
            m4_cloud_code = 100;  /* class cloud  [0.820] */
        
        /* Rule 3/37: (2207.5/402.6, lift 1.6) */
        if (m4_cloud_code == -1 && b1 > 1367 && b5 > 1480 && b2_var > 274632 && 
            ndsi_var <= 0.0177706)
            m4_cloud_code = 100;  /* class cloud  [0.817] */
        
        /* Rule 3/38: (28846/10368.5, lift 1.3) */
        if (m4_cloud_code == -1 && b1 > 1367 && b4_var <= 1090000 && 
            ndsi_var <= 0.0177706)
            m4_cloud_code = 100;  /* class cloud  [0.641] */
        
        if (m4_cloud_code == -1)
            m4_cloud_code = 50;

        
        /* 5th model run */
        /* Rule 4/1: (713.8, lift 2.1) */
        if (m5_cloud_code == -1 && b2 <= 4222 && b7 <= 1298 && 
            ndsi > 0.438529 && ndsi_var <= 0.0284024)
            m5_cloud_code = 50;  /* class cloud_free  [0.999] */
        
        /* Rule 4/2: (1328.1, lift 2.1) */
        if (m5_cloud_code == -1 && ndsi > 0.888811)
            m5_cloud_code = 50;  /* class cloud_free  [0.999] */
        
        /* Rule 4/3: (1293.4, lift 2.1) */
        if (m5_cloud_code == -1 && ndsi_var > 0.136021)
            m5_cloud_code = 50;  /* class cloud_free  [0.999] */
        
        /* Rule 4/4: (225.7, lift 2.1) */
        if (m5_cloud_code == -1 && b1 <= 2756 && b3 > 2842 && 
            ndsi_var <= 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.996] */
        
        /* Rule 4/5: (864.9/7.2, lift 2.1) */
        if (m5_cloud_code == -1 && b1 <= 2162  && b3 > 2095 && 
            ndsi_var > 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.991] */
        
        /* Rule 4/6: (219.2/5.1, lift 2.0) */
        if (m5_cloud_code == -1 && b2 <= 4222  && b7 <= 1298 && 
            ndvi <= 0.110778 && b2_var <= 3129.56)
            m5_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 4/7: (798.1/21.3, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 1187)
            m5_cloud_code = 50;  /* class cloud_free  [0.972] */
        
        /* Rule 4/8: (672.1/27.5, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 2162 && b5 > 1649 && 
            ndvi <= 0.0932514 && b1_var <= 170233 && b4_var > 5354.76 && 
            ndsi_var > 0.000488092 && ndsi_var <= 0.00313462)
            m5_cloud_code = 50;  /* class cloud_free  [0.958] */
        
        /* Rule 4/9: (830.5/39.5, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 3362 && b2 > 3250)
            m5_cloud_code = 50;  /* class cloud_free  [0.951] */
        
        /* Rule 4/10: (1202.4/69.1, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 2162 && b5 > 1819 && 
            b1_var <= 170233 && ndvi_var > 0.00355646 && 
            ndsi_var > 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.942] */
        
        /* Rule 4/11: (4352.5/423, lift 1.9) */
        if (m5_cloud_code == -1 && b1 <= 2162 && b5 > 1649 && 
            ndsi_var > 0.00313462)
            m5_cloud_code = 50;  /* class cloud_free  [0.903] */
        
        /* Rule 4/12: (1488/161, lift 1.9) */
        if (m5_cloud_code == -1 && ndvi <= 0.0564516 && ndsi <= -0.115053)
            m5_cloud_code = 50;  /* class cloud_free  [0.891] */
        
        /* Rule 4/13: (2569.8/299.6, lift 1.9) */
        if (m5_cloud_code == -1 && b1 <= 2162 && b5 > 1819 && 
            ndvi <= 0.236486 && b4_var > 5354.76 && b7_var <= 121237 && 
            ndsi_var > 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.883] */
        
        /* Rule 4/14: (1296.5/153.8, lift 1.8) */
        if (m5_cloud_code == -1 && b5 > 1770 && b7 <= 1298)
            m5_cloud_code = 50;  /* class cloud_free  [0.881] */
        
        /* Rule 4/15: (1612.4/206.9, lift 1.8) */
        if (m5_cloud_code == -1 && b1 > 2162 && b7 <= 1505 && 
            ndvi > 0.0749564 && b4_var > 5354.76 && ndsi_var > 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.871] */
        
        /* Rule 4/16: (4377.1/800.4, lift 1.7) */
        if (m5_cloud_code == -1 && b1 > 3362 && b5 <= 2173 && 
            ndsi_var > 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.817] */
        
        /* Rule 4/17: (302.5/77.1, lift 1.6) */
        if (m5_cloud_code == -1 && ndvi > 0.307271 && ndsi <= -0.115053 && 
            ndsi_var <= 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.744] */
        
        /* Rule 4/18: (28646.7/10725.7, lift 1.3) */
        if (m5_cloud_code == -1 && b1 <= 3362 && ndsi_var > 0.000488092)
            m5_cloud_code = 50;  /* class cloud_free  [0.626] */
        
        /* Rule 4/19: (1608.5, lift 1.9) */
        if (m5_cloud_code == -1 && b1 > 2756 && ndvi > 0.0564516 && 
            ndsi_var <= 0.000488092)
            m5_cloud_code = 100;  /* class cloud  [0.999] */
        
        /* Rule 4/20: (227.1, lift 1.9) */
        if (m5_cloud_code == -1 && b5 <= 1770 && ndvi > 0.110778 && 
            ndsi_var <= 0.000350481)
            m5_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 4/21: (193.5, lift 1.9) */
        if (m5_cloud_code == -1 && b3 <= 2095 && b7 > 1298 && 
            b1_var > 170233 && ndsi_var <= 0.00313462)
            m5_cloud_code = 100;  /* class cloud  [0.995] */
        
        /* Rule 4/22: (877.9/10, lift 1.9) */
        if (m5_cloud_code == -1 && b1 > 2162 && b3 <= 2003 && b7 > 1505 && 
            ndvi > 0.0749564 && ndsi_var <= 0.136021)
            m5_cloud_code = 100;  /* class cloud  [0.988] */
        
        /* Rule 4/23: (8938.3/275.6, lift 1.9) */
        if (m5_cloud_code == -1 && b1 > 3362 && b5 > 2173 && 
            ndsi_var <= 0.136021)
            m5_cloud_code = 100;  /* class cloud  [0.969] */
        
        /* Rule 4/24: (323.2/9.4, lift 1.9) */
        if (m5_cloud_code == -1 && b1 > 2162 && b2 <= 1898 && b7 > 1298 && 
            ndvi <= 0.0749564)
            m5_cloud_code = 100;  /* class cloud  [0.968] */
        
        /* Rule 4/25: (3254.9/166.5, lift 1.8) */
        if (m5_cloud_code == -1 && b2 > 4222 && ndsi <= 0.888811 && 
            ndvi_var <= 0.00194084 && ndsi_var <= 0.0284024)
            m5_cloud_code = 100;  /* class cloud  [0.949] */
        
        /* Rule 4/26: (1165.7/97.5, lift 1.8) */
        if (m5_cloud_code == -1 && b1 <= 3362 && b5 <= 1649 && b7 > 1298 && 
            ndsi_var <= 0.136021)
            m5_cloud_code = 100;  /* class cloud  [0.916] */
        
        /* Rule 4/27: (5129.8/488.6, lift 1.8) */
        if (m5_cloud_code == -1 && b1 > 2376 && b7 > 1505 && 
            ndvi > 0.0749564 && b4_var > 5354.76 && ndsi_var <= 0.136021)
            m5_cloud_code = 100;  /* class cloud  [0.905] */
        
        /* Rule 4/28: (126.5/12.1, lift 1.7) */
        if (m5_cloud_code == -1 && b1 > 1187 && b7 <= 1298 && 
            ndsi <= 0.438529 && b2_var > 356345 && ndvi_var <= 0.00194084 && 
            ndsi_var <= 0.0284024)
            m5_cloud_code = 100;  /* class cloud  [0.898] */
        
        /* Rule 4/29: (7965.1/822.8, lift 1.7) */
        if (m5_cloud_code == -1 && b7 > 1298 && ndvi <= 0.307271 && 
            ndsi_var <= 0.000488092)
            m5_cloud_code = 100;  /* class cloud  [0.897] */
        
        /* Rule 4/30: (2268.3/362.9, lift 1.6) */
        if (m5_cloud_code == -1 && b1 > 1187 && b2 <= 4222 && b5 <= 1770 && 
            b7 > 1192 && ndsi <= 0.438529 && b2_var <= 356345 && 
            ndvi_var <= 0.00194084 && ndsi_var <= 0.0284024)
            m5_cloud_code = 100;  /* class cloud  [0.840] */
        
        /* Rule 4/31: (3018/556.2, lift 1.6) */
        if (m5_cloud_code == -1 && b1 > 2629 && b2 <= 3250 && b5 > 1649 && 
            b7 > 1298 && b7 <= 2681 && ndsi_var <= 0.136021)
            m5_cloud_code = 100;  /* class cloud  [0.815] */
        
        /* Rule 4/32: (2135.4/411.8, lift 1.6) */
        if (m5_cloud_code == -1 && b1 > 1187 && b3 <= 1873 && b5 <= 1770 && 
            ndvi > 0.110778 && ndvi_var <= 0.00194084 && ndsi_var <= 0.0284024)
            m5_cloud_code = 100;  /* class cloud  [0.807] */
        
        /* Rule 4/33: (10111.7/2158.9, lift 1.5) */
        if (m5_cloud_code == -1 && b7 > 1298 && ndvi > 0.0932514 && 
            ndsi_var <= 0.00313462)
            m5_cloud_code = 100;  /* class cloud  [0.786] */
        
        /* Rule 4/34: (1258.5/278.7, lift 1.5) */
        if (m5_cloud_code == -1 && b1 > 1187 && b2 <= 1276 && b5 <= 1770 && 
            b7 <= 1192 && ndsi <= 0.438529 && b2_var <= 356345 && 
            ndvi_var <= 0.00194084 && ndsi_var <= 0.0284024)
            m5_cloud_code = 100;  /* class cloud  [0.778] */
        
        if (m5_cloud_code == -1)
            m5_cloud_code = 50;
 
        /* Use the maximum cloud code as the cloud value from this set of
           model runs */
        conserv_cloud_code = m1_cloud_code;
        if (m2_cloud_code > conserv_cloud_code)
            conserv_cloud_code = m2_cloud_code;
        if (m3_cloud_code > conserv_cloud_code)
            conserv_cloud_code = m3_cloud_code;
        if (m4_cloud_code > conserv_cloud_code)
            conserv_cloud_code = m4_cloud_code;
        if (m5_cloud_code > conserv_cloud_code)
            conserv_cloud_code = m5_cloud_code;


        /*** Conservative cloud mask model run without the variances ***/
        /* Initialize the cloud mask */
        m1_lim_cloud_code = -1;
        m2_lim_cloud_code = -1;
        m3_lim_cloud_code = -1;
        m4_lim_cloud_code = -1;
        m5_lim_cloud_code = -1;

        /* Implement rule-based model for cases where variance calculation is
           not possible, such as SLC-off areas and areas near the edge of the
           image */
        /* 1st model run */
        /* Rule 0/1: (10764/235, lift 2.1) */
        if (m1_lim_cloud_code == -1 && b3 > 1424 && b7 <= 1187)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.978] */
        
        /* Rule 0/2: (753/16, lift 2.1) */
        if (m1_lim_cloud_code == -1 && b1 <= 2999 && ndvi <= 0.158153 &&
            ndsi <= -0.212327)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.977] */
        
        /* Rule 0/3: (447/10, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 2079 && b3 > 2080)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.976] */
        
        /* Rule 0/4: (1763/48, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 > 6251 && b4 <= 5012 && b7 <= 2077)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.972] */
        
        /* Rule 0/5: (960/27, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 2483 && b3 > 2408 &&
            ndsi > -0.212327)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.971] */
        
        /* Rule 0/6: (728/29, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 2999 && b2 > 2771 &&
            ndvi <= 0.0912175)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.959] */
        
        /* Rule 0/7: (4653/195, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b4 > 4172 && b7 <= 1464)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.958] */
        
        /* Rule 0/8: (5826/279, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 3532 && b3 > 1727 && b4 > 1928 &&
            b7 <= 1464)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.952] */
        
        /* Rule 0/9: (1917/93, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 2999 && b2 > 1970 && b5 <= 1937 &&
            b7 <= 1619 && ndvi <= 0.0912175)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.951] */
        
        /* Rule 0/10: (1198/63, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 > 2999 && b1 <= 3589 && b4 > 3588 &&
            b7 <= 2077)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.947] */
        
        /* Rule 0/11: (1293/71, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 2999 && b4 > 3173 &&
            ndsi > -0.0444104)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.944] */
        
        /* Rule 0/12: (2832/157, lift 2.0) */
        if (m1_lim_cloud_code == -1 && b1 <= 1600 && b4 <= 1928)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.944] */
        
        /* Rule 0/13: (2380/185, lift 1.9) */
        if (m1_lim_cloud_code == -1 && b4 <= 2729 && ndsi <= -0.212327)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.922] */
        
        /* Rule 0/14: (2208/176, lift 1.9) */
        if (m1_lim_cloud_code == -1 && b1 <= 2638 && b5 > 1937 &&
            ndvi <= 0.0912175)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.920] */
        
        /* Rule 0/15: (1742/156, lift 1.9) */
        if (m1_lim_cloud_code == -1 && b7 <= 1799 && ndsi <= -0.212327)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.910] */
        
        /* Rule 0/16: (1724/180, lift 1.9) */
        if (m1_lim_cloud_code == -1 && b1 <= 1993 && b5 > 1937 && b7 > 1464 &&
            ndvi <= 0.210062)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.895] */
        
        /* Rule 0/17: (23277/7218, lift 1.4) */
        if (m1_lim_cloud_code == -1 && b1 <= 2999)
            m1_lim_cloud_code = 50;  /* class cloud_free  [0.690] */
        
        /* Rule 0/18: (3593/38, lift 1.9) */
        if (m1_lim_cloud_code == -1 && b1 > 2483 && ndvi > 0.0912175 &&
            ndsi <= -0.0444104)
            m1_lim_cloud_code = 100;  /* class cloud  [0.989] */
        
        /* Rule 0/19: (5613/71, lift 1.9) */
        if (m1_lim_cloud_code == -1 && b7 > 1266 && ndsi > 0.491986)
            m1_lim_cloud_code = 100;  /* class cloud  [0.987] */
        
        /* Rule 0/20: (1100/125, lift 1.7) */
        if (m1_lim_cloud_code == -1 && b7 > 1464 && ndvi > 0.210062 &&
            ndsi > -0.212327)
            m1_lim_cloud_code = 100;  /* class cloud  [0.886] */
        
        /* Rule 0/21: (267/39, lift 1.6) */
        if (m1_lim_cloud_code == -1 && b1 > 1600 && b3 <= 1299 && b4 <= 1928 &&
            b7 > 919 && b7 <= 1464)
            m1_lim_cloud_code = 100;  /* class cloud  [0.851] */
        
        /* Rule 0/22: (193/35, lift 1.6) */
        if (m1_lim_cloud_code == -1 && b1 > 1451 && b3 <= 1139 && b7 <= 919 &&
            ndvi > 0.0602801 && ndvi <= 0.421493 && ndsi <= 0.169811)
            m1_lim_cloud_code = 100;  /* class cloud  [0.815] */
        
        /* Rule 0/23: (32298/7826, lift 1.4) */
        if (m1_lim_cloud_code == -1 && b1 > 1802 && b7 > 919)
            m1_lim_cloud_code = 100;  /* class cloud  [0.758] */
        
        if (m1_lim_cloud_code == -1)
            m1_lim_cloud_code = 100;


        /* 2nd model run */
        /* Rule 1/1: (801.4/3.8, lift 2.1) */
        if (m2_lim_cloud_code == -1 && b4 > 6829 && b5 <= 3315)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.994] */
        
        /* Rule 1/2: (657/4.5, lift 2.1) */
        if (m2_lim_cloud_code == -1 && ndvi > 0.44504)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.992] */
        
        /* Rule 1/3: (220.2/1.5, lift 2.1) */
        if (m2_lim_cloud_code == -1 && b1 <= 3544 && b4 > 4212 && b5 > 1136 &&
            b5 <= 3315)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.989] */
        
        /* Rule 1/4: (6122.8/130.4, lift 2.1) */
        if (m2_lim_cloud_code == -1 && b3 > 1258 && b5 <= 1136)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.979] */
        
        /* Rule 1/5: (988.6/45.5, lift 2.1) */
        if (m2_lim_cloud_code == -1 && b1 <= 3544 && b5 > 1136 && b7 <= 1067 &&
            ndsi > 0.143936)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.953] */
        
        /* Rule 1/6: (3721.5/193.8, lift 2.0) */
        if (m2_lim_cloud_code == -1 && b1 > 6251 && b7 <= 1683)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.948] */
        
        /* Rule 1/7: (504.2/47.8, lift 1.9) */
        if (m2_lim_cloud_code == -1 && b1 > 1555 && b1 <= 2321 && b4 <= 3520 &&
            b5 > 2467 && b7 <= 1726)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.904] */
        
        /* Rule 1/8: (1057/103.5, lift 1.9) */
        if (m2_lim_cloud_code == -1 && b1 <= 3544 && b2 > 3236 && b5 <= 3315)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.901] */
        
        /* Rule 1/9: (671.1/76.4, lift 1.9) */
        if (m2_lim_cloud_code == -1 && b1 <= 2835 && b3 > 2538 &&
            ndsi > -0.106196 && ndsi <= 0.143936)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.885] */
        
        /* Rule 1/10: (650.1/75.1, lift 1.9) */
        if (m2_lim_cloud_code == -1 && b4 <= 4001 && b5 > 3315 &&
            ndvi <= 0.0617448)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.883] */
        
        /* Rule 1/11: (1812.6/334.2, lift 1.8) */
        if (m2_lim_cloud_code == -1 && b1 <= 2082 && b3 > 1839 && b4 <= 3520)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.815] */
        
        /* Rule 1/12: (31120.8/12936.8, lift 1.3) */
        if (m2_lim_cloud_code == -1 && b4 <= 4001)
            m2_lim_cloud_code = 50;  /* class cloud_free  [0.584] */
        
        /* Rule 1/13: (11094.4/118.4, lift 1.8) */
        if (m2_lim_cloud_code == -1 && b4 > 4001 && b5 > 3315)
            m2_lim_cloud_code = 100;  /* class cloud  [0.989] */
        
        /* Rule 1/14: (519/12.9, lift 1.8) */
        if (m2_lim_cloud_code == -1 && b1 > 2043 && b2 <= 1763 &&
            ndsi <= 0.143936)
            m2_lim_cloud_code = 100;  /* class cloud  [0.973] */
        
        /* Rule 1/15: (12178.2/334.3, lift 1.8) */
        if (m2_lim_cloud_code == -1 && b1 > 3544 && b7 > 1683)
            m2_lim_cloud_code = 100;  /* class cloud  [0.972] */
        
        /* Rule 1/16: (3402.2/106.4, lift 1.8) */
        if (m2_lim_cloud_code == -1 && b1 > 2599 && b5 > 3315 &&
            ndvi > 0.0617448)
            m2_lim_cloud_code = 100;  /* class cloud  [0.968] */
        
        /* Rule 1/17: (9675.1/785.7, lift 1.7) */
        if (m2_lim_cloud_code == -1 && b1 > 2835 && ndsi <= 0.143936)
            m2_lim_cloud_code = 100;  /* class cloud  [0.919] */
        
        /* Rule 1/18: (1278.7/104.9, lift 1.7) */
        if (m2_lim_cloud_code == -1 && b1 > 3264 && b2 <= 3334 && b7 > 1067)
            m2_lim_cloud_code = 100;  /* class cloud  [0.917] */
        
        /* Rule 1/19: (564.8/53.3, lift 1.7) */
        if (m2_lim_cloud_code == -1 && b1 > 2151 && b2 <= 1898 && b5 > 1537)
            m2_lim_cloud_code = 100;  /* class cloud  [0.904] */
        
        /* Rule 1/20: (228.7/24.9, lift 1.7) */
        if (m2_lim_cloud_code == -1 && b1 > 2082 && b3 <= 1994 &&
            ndvi > 0.179141 && ndsi <= 0.143936)
            m2_lim_cloud_code = 100;  /* class cloud  [0.888] */
        
        /* Rule 1/21: (4430.9/528.7, lift 1.6) */
        if (m2_lim_cloud_code == -1 && b1 > 3544 && b1 <= 6251 && b4 <= 6829 &&
            b5 > 1136)
            m2_lim_cloud_code = 100;  /* class cloud  [0.881] */
        
        /* Rule 1/22: (179.6/32.6, lift 1.5) */
        if (m2_lim_cloud_code == -1 && b2 <= 2204 && b4 > 3520 && b5 <= 3315 &&
            ndvi <= 0.44504 && ndsi <= 0.143936)
            m2_lim_cloud_code = 100;  /* class cloud  [0.815] */
        
        /* Rule 1/23: (481.3/131.1, lift 1.4) */
        if (m2_lim_cloud_code == -1 && b1 > 1800 && b1 <= 2082 && b3 <= 1839 &&
            b5 <= 3315 && b7 > 1342 && ndvi > 0.179141 && ndvi <= 0.207702)
            m2_lim_cloud_code = 100;  /* class cloud  [0.727] */
        
        /* Rule 1/24: (4787.7/1392.3, lift 1.3) */
        if (m2_lim_cloud_code == -1 && b1 > 2321 && b2 <= 3236 && b4 <= 4212 &&
            b5 <= 3315 && ndsi <= 0.143936)
            m2_lim_cloud_code = 100;  /* class cloud  [0.709] */
        
        /* Rule 1/25: (2379.8/758.4, lift 1.3) */
        if (m2_lim_cloud_code == -1 && b1 > 1555 && b1 <= 1800 && b3 <= 1839 &&
            b5 > 1537 && ndvi > 0.179141 && ndvi <= 0.44504 &&
            ndsi <= -0.0472716)
            m2_lim_cloud_code = 100;  /* class cloud  [0.681] */
        
        /* Rule 1/26: (2563.3/868.5, lift 1.2) */
        if (m2_lim_cloud_code == -1 && b1 > 1182 && b5 > 1136 && b5 <= 1537 &&
            b7 > 542 && ndvi <= 0.44504 && ndsi <= 0.143936)
            m2_lim_cloud_code = 100;  /* class cloud  [0.661] */
        
        /* Rule 1/27: (897.2/366.5, lift 1.1) */
        if (m2_lim_cloud_code == -1 && b3 <= 1258 && b5 <= 1136 && b7 > 542 &&
            ndvi <= 0.44504)
            m2_lim_cloud_code = 100;  /* class cloud  [0.591] */
        
        if (m2_lim_cloud_code == -1)
            m2_lim_cloud_code = 50;  /* Default class: cloud_free */


        /* 3rd model run */
        /* Rule 2/1: (286.3, lift 2.0) */
        if (m3_lim_cloud_code == -1 && b1 <= 2497 && b3 > 2156 && b3 <= 2276 &&
            b7 <= 1626 && ndvi > 0.0650248)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.997] */
        
        /* Rule 2/2: (2084.7/6.1, lift 2.0) */
        if (m3_lim_cloud_code == -1 && ndsi > 0.881556)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.997] */
        
        /* Rule 2/3: (347/7.9, lift 2.0) */
        if (m3_lim_cloud_code == -1 && b2 > 3653 && b7 <= 1626 &&
            ndvi > 0.0650248)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.974] */
        
        /* Rule 2/4: (642.5/16.7, lift 2.0) */
        if (m3_lim_cloud_code == -1 && b1 <= 2497 && b3 > 1979 && b7 <= 1626 &&
            ndsi <= 0.0546875)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 2/5: (3129.8/84.8, lift 2.0) */
        if (m3_lim_cloud_code == -1 && b2 > 1011 && b7 <= 486)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 2/6: (264.7/8.4, lift 2.0) */
        if (m3_lim_cloud_code == -1 && b1 <= 1966 && b5 > 1707 && b7 <= 1626 &&
            ndvi > 0.0650248 && ndvi <= 0.12763)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.965] */
        
        /* Rule 2/7: (152/6.1, lift 1.9) */
        if (m3_lim_cloud_code == -1 && b5 > 2507 && b7 <= 1626 &&
            ndvi > 0.0650248 && ndsi > -0.254112)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.954] */
        
        /* Rule 2/8: (1358.7/66.7, lift 1.9) */
        if (m3_lim_cloud_code == -1 && b2 <= 1011)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.950] */
        
        /* Rule 2/9: (7074/402.6, lift 1.9) */
        if (m3_lim_cloud_code == -1 && b3 > 1265 && b7 <= 1008)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.943] */
        
        /* Rule 2/10: (918.4/59.4, lift 1.9) */
        if (m3_lim_cloud_code == -1 && b3 <= 1265 && ndvi > 0.403095)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.934] */
        
        /* Rule 2/11: (2111.8/153.8, lift 1.9) */
        if (m3_lim_cloud_code == -1 && b1 <= 1725 && b7 <= 1286 &&
            ndvi <= 0.139127)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.927] */
        
        /* Rule 2/12: (1275.7/107.6, lift 1.9) */
        if (m3_lim_cloud_code == -1 && b1 <= 1725 && b3 > 1205 && b4 <= 2458 &&
            b7 <= 1286)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.915] */
        
        /* Rule 2/13: (3639.5/324.7, lift 1.8) */
        if (m3_lim_cloud_code == -1 && b1 > 1842 && b2 <= 4235 && b4 > 1742 &&
            b7 <= 1286 && ndvi <= 0.0979757)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.911] */
        
        /* Rule 2/14: (1170.6/111.3, lift 1.8) */
        if (m3_lim_cloud_code == -1 && b2 <= 4235 && b3 > 1265 && b5 > 1747 &&
            b7 <= 1286)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.904] */
        
        /* Rule 2/15: (126.8/12.2, lift 1.8) */
        if (m3_lim_cloud_code == -1 && b4 > 6150 && b5 <= 3315 &&
            ndsi <= 0.49737)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.898] */
        
        /* Rule 2/16: (1897.2/223.5, lift 1.8) */
        if (m3_lim_cloud_code == -1 && b5 <= 3315 && ndsi > 0.368125 &&
            ndsi <= 0.49737)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.882] */
        
        /* Rule 2/17: (775.9/102.5, lift 1.8) */
        if (m3_lim_cloud_code == -1 && b1 <= 3284 && b3 > 2276 && b5 > 1707 &&
            b7 <= 1626 && ndvi > 0.0650248)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.867] */
        
        /* Rule 2/18: (441.9/66.5, lift 1.7) */
        if (m3_lim_cloud_code == -1 && b4 <= 3074 && b5 > 3315)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.848] */
        
        /* Rule 2/19: (972/154.3, lift 1.7) */
        if (m3_lim_cloud_code == -1 && b1 <= 3418 && b5 > 3315 &&
            ndvi <= 0.079329)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.841] */
        
        /* Rule 2/20: (1000.4/161.7, lift 1.7) */
        if (m3_lim_cloud_code == -1 && b3 <= 1205 && b7 <= 1286 &&
            ndsi <= -0.0934991)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.838] */
        
        /* Rule 2/21: (168.1/31.5, lift 1.6) */
        if (m3_lim_cloud_code == -1 && b5 > 2885 && b7 <= 2112 &&
            ndvi > 0.0650248 && ndsi > -0.254112)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.809] */
        
        /* Rule 2/22: (3523.5/696, lift 1.6) */
        if (m3_lim_cloud_code == -1 && b1 <= 1966 && b3 > 1220 && b4 <= 2470 &&
            b7 <= 1458)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.802] */
        
        /* Rule 2/23: (795.3/160.6, lift 1.6) */
        if (m3_lim_cloud_code == -1 && b1 <= 1966 && b4 > 2470 && b4 <= 3217 &&
            b7 > 1286 && b7 <= 1626)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.797] */
        
        /* Rule 2/24: (557.4/120.5, lift 1.6) */
        if (m3_lim_cloud_code == -1 && b5 > 2528 && b7 <= 1805 &&
            ndvi > 0.0650248 && ndsi > -0.254112)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.783] */
        
        /* Rule 2/25: (8172.8/1851.2, lift 1.6) */
        if (m3_lim_cloud_code == -1 && b1 <= 3376 && ndvi <= 0.0650248)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.773] */
        
        /* Rule 2/26: (1800.2/430.1, lift 1.5) */
        if (m3_lim_cloud_code == -1 && b4 <= 2856 && ndsi <= -0.254112)
            m3_lim_cloud_code = 50;  /* class cloud_free  [0.761] */
        
        /* Rule 2/27: (308.2/41.3, lift 1.7) */
        if (m3_lim_cloud_code == -1 && b1 > 1725 && b3 <= 1265 && b7 > 486)
            m3_lim_cloud_code = 100;  /* class cloud  [0.864] */
        
        /* Rule 2/28: (278.4/47.5, lift 1.6) */
        if (m3_lim_cloud_code == -1 && b1 <= 1725 && b2 > 1011 && b3 <= 1265 &&
            b4 > 2458 && b7 > 486 && b7 <= 1286 && ndvi <= 0.403095)
            m3_lim_cloud_code = 100;  /* class cloud  [0.827] */
        
        /* Rule 2/29: (1011.9/222.6, lift 1.5) */
        if (m3_lim_cloud_code == -1 && b2 > 1011 && b5 <= 1707 && b7 > 1286 &&
            ndsi <= 0.368125)
            m3_lim_cloud_code = 100;  /* class cloud  [0.779] */
        
        /* Rule 2/30: (515.1/136.5, lift 1.5) */
        if (m3_lim_cloud_code == -1 && b1 > 1842 && b4 <= 1742 && b5 <= 1747 &&
            b7 > 1008)
            m3_lim_cloud_code = 100;  /* class cloud  [0.734] */
        
        /* Rule 2/31: (46915.3/22120.7, lift 1.0) */
        if (m3_lim_cloud_code == -1 && ndsi <= 0.881556)
            m3_lim_cloud_code = 100;  /* class cloud  [0.528] */
        
        if (m3_lim_cloud_code == -1)
            m3_lim_cloud_code = 100;  /* Default class: cloud */


        /* 4th model run */
        /* Rule 3/1: (709.4/15.7, lift 1.9) */
        if (m4_lim_cloud_code == -1 && b1 <= 2118 && b2 > 2012)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.976] */
        
        /* Rule 3/2: (196.5/4.4, lift 1.8) */
        if (m4_lim_cloud_code == -1 && b1 <= 3234 && b2 > 3108 && b5 > 3315)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 3/3: (1753.6/96.9, lift 1.8) */
        if (m4_lim_cloud_code == -1 && b1 <= 2467 && b3 > 2364)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.944] */
        
        /* Rule 3/4: (1047.4/69.7, lift 1.8) */
        if (m4_lim_cloud_code == -1 && b1 <= 2654 && b5 > 3315 &&
            ndvi <= 0.1719)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.933] */
        
        /* Rule 3/5: (2283.3/179.7, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b1 <= 2467 && b2 > 1983 && b5 > 1706 &&
            ndvi <= 0.114312)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.921] */
        
        /* Rule 3/6: (1552.5/128.6, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b1 <= 2762 && b3 > 2437 &&
            ndsi > -0.115228)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.917] */
        
        /* Rule 3/7: (3784.7/387.6, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b3 > 1665 && b7 <= 1347 &&
            ndvi > -0.0240739 && ndvi <= 0.131339)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.897] */
        
        /* Rule 3/8: (3288.9/361.2, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b1 > 6251 && b5 <= 3112 && b7 <= 2077)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.890] */
        
        /* Rule 3/9: (1266.7/159.1, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b1 <= 3160 && b3 > 2834 && b5 <= 3112)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.874] */
        
        /* Rule 3/10: (1701.1/238.9, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b4 > 5404 && b5 <= 3112 && b7 <= 2077)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.859] */
        
        /* Rule 3/11: (388.7/55.3, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 <= 2467 && b5 > 1706 && b7 > 1347 &&
            ndvi > 0.114312 && ndsi > -0.0444104)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.856] */
        
        /* Rule 3/12: (902.2/145.9, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 <= 2118 && b4 <= 3249 && b5 > 2406 &&
            b7 <= 1726)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.838] */
        
        /* Rule 3/13: (2501.9/446.9, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 <= 2118 && b4 <= 2223 && b5 > 1706)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.821] */
        
        /* Rule 3/14: (2864.7/535.5, lift 1.5) */
        if (m4_lim_cloud_code == -1 && b1 <= 1873 && b5 > 1138 && b7 <= 1347 &&
            ndvi <= 0.317016)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.813] */
        
        /* Rule 3/15: (1899.1/420.7, lift 1.5) */
        if (m4_lim_cloud_code == -1 && b1 <= 3544 && b2 > 3043 && b5 <= 3112)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.778] */
        
        /* Rule 3/16: (1598/394.3, lift 1.4) */
        if (m4_lim_cloud_code == -1 && b1 <= 2654 && b5 > 3315)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.753] */
        
        /* Rule 3/17: (737.3/186, lift 1.4) */
        if (m4_lim_cloud_code == -1 && b1 <= 1684 && b2 > 1342 && b5 <= 2406 &&
            b7 > 1347)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.747] */
        
        /* Rule 3/18: (38683.4/14494.6, lift 1.2) */
        if (m4_lim_cloud_code == -1 && b5 <= 3315)
            m4_lim_cloud_code = 50;  /* class cloud_free  [0.625] */
        
        /* Rule 3/19: (7194.6, lift 2.1) */
        if (m4_lim_cloud_code == -1 && b1 > 3544 && b7 > 2077)
            m4_lim_cloud_code = 100;   /* class cloud  [1.000] */
        
        /* Rule 3/20: (323/3.8, lift 2.1) */
        if (m4_lim_cloud_code == -1 && b1 > 2467 && b5 <= 3315 &&
            ndsi <= -0.115228)
            m4_lim_cloud_code = 100;  /* class cloud  [0.985] */
        
        /* Rule 3/21: (563.4/23.1, lift 2.0) */
        if (m4_lim_cloud_code == -1 && b1 > 3160 && b2 <= 3043 && b7 > 1347)
            m4_lim_cloud_code = 100;  /* class cloud  [0.957] */
        
        /* Rule 3/22: (8425.6/406.2, lift 2.0) */
        if (m4_lim_cloud_code == -1 && b1 > 2654 && b5 > 3315)
            m4_lim_cloud_code = 100;  /* class cloud  [0.952] */
        
        /* Rule 3/23: (8942.7/654, lift 2.0) */
        if (m4_lim_cloud_code == -1 && b1 > 3544 && b7 > 1347 &&
            ndvi > -0.626148)
            m4_lim_cloud_code = 100;  /* class cloud  [0.927] */
        
        /* Rule 3/24: (5373.8/411.2, lift 2.0) */
        if (m4_lim_cloud_code == -1 && b5 > 1157 && ndvi > -0.564211 &&
            ndvi <= -0.0240739 && ndsi <= 0.881468)
            m4_lim_cloud_code = 100;  /* class cloud  [0.923] */
        
        /* Rule 3/25: (1531.5/144.8, lift 1.9) */
        if (m4_lim_cloud_code == -1 && b1 > 2762 && b3 <= 2834 && b7 > 1347)
            m4_lim_cloud_code = 100;  /* class cloud  [0.905] */
        
        /* Rule 3/26: (115.6/13, lift 1.9) */
        if (m4_lim_cloud_code == -1 && b1 > 1526 && b2 <= 1342 && b5 <= 3315 &&
            b7 > 1347 && ndvi > 0.114312)
            m4_lim_cloud_code = 100;  /* class cloud  [0.881] */
        
        /* Rule 3/27: (1386.2/198.3, lift 1.8) */
        if (m4_lim_cloud_code == -1 && b1 > 2467 && b3 <= 2437 && b7 > 1347)
            m4_lim_cloud_code = 100;  /* class cloud  [0.856] */
        
        /* Rule 3/28: (107.2/23.2, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b5 > 1136 && b5 <= 1138)
            m4_lim_cloud_code = 100;  /* class cloud  [0.778] */
        
        /* Rule 3/29: (1097.6/247.2, lift 1.7) */
        if (m4_lim_cloud_code == -1 && b1 > 2123 && b2 <= 1983 && b7 > 1347)
            m4_lim_cloud_code = 100;  /* class cloud  [0.774] */
        
        /* Rule 3/30: (4230.3/989.1, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 > 2118 && b7 > 1347 &&
            ndvi > 0.114312 && ndsi <= -0.0444104)
            m4_lim_cloud_code = 100;  /* class cloud  [0.766] */
        
        /* Rule 3/31: (483.4/114.7, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 > 1526 && b2 <= 2012 && b4 > 3249 &&
            b5 <= 3315 && b7 > 1347 && ndvi <= 0.44504)
            m4_lim_cloud_code = 100;  /* class cloud  [0.762] */
        
        /* Rule 3/32: (516.7/123.6, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 > 1526 && b1 <= 2467 && b3 <= 2364 &&
            b5 <= 1706 && b7 > 1347)
            m4_lim_cloud_code = 100;  /* class cloud  [0.760] */
        
        /* Rule 3/33: (1601.7/412.5, lift 1.6) */
        if (m4_lim_cloud_code == -1 && b1 > 1873 && b3 <= 1665 && b5 > 1157 &&
            ndvi > -0.0240739 && ndvi <= 0.317016 && ndsi > -0.184974)
            m4_lim_cloud_code = 100;  /* class cloud  [0.742] */
        
        /* Rule 3/34: (2091.4/633.1, lift 1.5) */
        if (m4_lim_cloud_code == -1 && b1 > 1684 && b4 > 2223 && b5 <= 2406 &&
            b7 > 1347 && ndvi > 0.114312 && ndvi <= 0.44504 &&
            ndsi <= -0.0800831)
            m4_lim_cloud_code = 100;  /* class cloud  [0.697] */
        
        /* Rule 3/35: (2531.2/772, lift 1.5) */
        if (m4_lim_cloud_code == -1 && b1 > 1526 && b3 <= 2364 && b5 <= 2626 &&
            b7 > 1726 && ndvi > 0.114312)
            m4_lim_cloud_code = 100;  /* class cloud  [0.695] */
        
        /* Rule 3/36: (781.6/292.1, lift 1.3) */
        if (m4_lim_cloud_code == -1 && b7 > 1086 && ndvi > 0.131339 &&
            ndsi > 0.0375777 && ndsi <= 0.416984)
            m4_lim_cloud_code = 100;  /* class cloud  [0.626] */
        
        /* Rule 3/37: (4079.4/1565.3, lift 1.3) */
        if (m4_lim_cloud_code == -1 && b1 > 1526 && b3 <= 2364 && b5 <= 2406 &&
            b7 > 1347 && ndvi > 0.114312 && ndvi <= 0.44504)
            m4_lim_cloud_code = 100;  /* class cloud  [0.616] */

        if (m4_lim_cloud_code == -1)
            m4_lim_cloud_code = 50;  /* Default class: cloud_free */

        
        /* 5th model run */
        /* Rule 4/1: (400.1, lift 2.1) */
        if (m5_lim_cloud_code == -1 && b7 <= 1461 && ndsi <= -0.224928)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.998] */
        
        /* Rule 4/2: (815.6/0.4, lift 2.1) */
        if (m5_lim_cloud_code == -1 && b1 <= 2601 && b3 > 2260 && b7 <= 1623)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.998] */
        
        /* Rule 4/3: (273, lift 2.1) */
        if (m5_lim_cloud_code == -1 && b1 <= 1511 && b2 > 1317 && b5 > 1516)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.996] */
        
        /* Rule 4/4: (1665.9/5.6, lift 2.1) */
        if (m5_lim_cloud_code == -1 && ndsi > 0.85271)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.996] */
        
        /* Rule 4/5: (677.9/4.5, lift 2.1) */
        if (m5_lim_cloud_code == -1 && b1 <= 2601 && ndvi <= 0.107894 &&
            ndsi <= -0.177533)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.992] */
        
        /* Rule 4/6: (1275.5/13, lift 2.1) */
        if (m5_lim_cloud_code == -1 && b1 <= 1645 && ndvi <= 0.107894)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.989] */
        
        /* Rule 4/7: (586.3/6.4, lift 2.1) */
        if (m5_lim_cloud_code == -1 && b2 <= 1027 && b5 > 796)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.987] */
        
        /* Rule 4/8: (2953.5/46.3, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b5 <= 796)
            m5_lim_cloud_code = 50;   /* class cloud_free  [0.984] */
        
        /* Rule 4/9: (1042.7/19.5, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 <= 2601 && b3 > 2569)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.980] */
        
        /* Rule 4/10: (988.3/19.6, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 <= 1802 && b2 > 1292 &&
            ndvi <= 0.107894)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.979] */
        
        /* Rule 4/11: (268/4.7, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b2 > 1387 && b3 <= 1461 &&
            ndvi <= 0.412224 && ndsi <= -0.22841)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.979] */
        
        /* Rule 4/12: (719.9/16, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 <= 2005 && b3 > 1949)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.976] */
        
        /* Rule 4/13: (565.2/15.4, lift 2.0) */
        if (m5_lim_cloud_code == -1 && ndvi > 0.412224)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.971] */
        
        /* Rule 4/14: (3293.8/99.1, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b7 <= 579)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.970] */
        
        /* Rule 4/15: (966.4/33.4, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 <= 3037 && b2 > 2885 && b5 <= 4374)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.964] */
        
        /* Rule 4/16: (744.4/28.7, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b2 <= 1317 && b5 > 1516 &&
            ndvi > 0.107894)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.960] */
        
        /* Rule 4/17: (864.7/38.2, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 > 2601 && b1 <= 3037 &&
            ndsi > 0.143982)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.955] */
        
        /* Rule 4/18: (1667.8/114.9, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b1 <= 3037 && b3 > 2909)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.931] */
        
        /* Rule 4/19: (441.2/30.5, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b4 > 7015 && b5 <= 4374)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.929] */
        
        /* Rule 4/20: (129.7/8.6, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b4 <= 4251 && b5 > 4374)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.927] */
        
        /* Rule 4/21: (1596/121.5, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b1 <= 2247 && b2 > 1763 && b5 > 1493 &&
            ndvi <= 0.107894)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.923] */
        
        /* Rule 4/22: (6175.3/491.4, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b3 > 1223 && b7 <= 1084 &&
            ndsi > -0.0853392)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.920] */
        
        /* Rule 4/23: (3766.2/366.4, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b2 <= 3959 && ndsi > 0.267478)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.902] */
        
        /* Rule 4/24: (414.9/46.3, lift 1.8) */
        if (m5_lim_cloud_code == -1 && b5 <= 1516 && ndsi <= -0.0853392)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.887] */
        
        /* Rule 4/25: (3150.8/363.8, lift 1.8) */
        if (m5_lim_cloud_code == -1 && b1 <= 2005 && b5 > 1516 &&
            ndvi <= 0.175138)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.884] */
        
        /* Rule 4/26: (1868.5/321.3, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b1 <= 2005 && b4 <= 2856 &&
            ndsi <= -0.22841)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.828] */
        
        /* Rule 4/27: (126.6/21.5, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b2 > 3959 && b5 <= 4374 &&
            ndvi > 0.0799109)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.825] */
        
        /* Rule 4/28: (2695.3/480.2, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b1 <= 2601 && b2 > 1906 &&
            ndvi <= 0.107894 && ndsi > -0.115822)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.822] */
        
        /* Rule 4/29: (3038/650.2, lift 1.6) */
        if (m5_lim_cloud_code == -1 && b1 <= 2394 && b5 > 2533 &&
            ndvi <= 0.201976)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.786] */
        
        /* Rule 4/30: (13065.9/2877.5, lift 1.6) */
        if (m5_lim_cloud_code == -1 && b3 > 1267 && b7 <= 1461 &&
            ndvi <= 0.412224)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.780] */
        
        /* Rule 4/31: (3319.6/846.9, lift 1.5) */
        if (m5_lim_cloud_code == -1 && b1 <= 2392 && b2 > 1698 && b5 > 1516 &&
            b7 <= 1623 && ndvi <= 0.412224)
            m5_lim_cloud_code = 50;  /* class cloud_free  [0.745] */
        
        /* Rule 4/32: (477.2, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 > 2005 && b2 <= 1698 && b5 > 1516)
            m5_lim_cloud_code = 100;  /* class cloud  [0.998] */
        
        /* Rule 4/33: (261.5, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 > 2247 && b2 <= 1906 && b7 > 1082)
            m5_lim_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 4/34: (546.4/1.6, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 > 1645 && b2 <= 1292 && b5 > 796)
            m5_lim_cloud_code = 100;  /* class cloud  [0.995] */
        
        /* Rule 4/35: (286.3/3.2, lift 2.0) */
        if (m5_lim_cloud_code == -1 && b1 > 2040 && b2 <= 1763 && b7 > 1082 &&
            ndvi <= 0.107894)
            m5_lim_cloud_code = 100;  /* class cloud  [0.986] */
        
        /* Rule 4/36: (262.7/10.2, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b5 <= 1516 && b7 > 1084 &&
            ndvi > 0.107894 && ndsi > -0.0853392 && ndsi <= 0.267478)
            m5_lim_cloud_code = 100;  /* class cloud  [0.958] */
        
        /* Rule 4/37: (375.3/16.8, lift 1.9) */
        if (m5_lim_cloud_code == -1 && b1 > 2392 && b3 <= 2260 && b5 > 1516 &&
            ndvi > 0.107894)
            m5_lim_cloud_code = 100;  /* class cloud  [0.953] */
        
        /* Rule 4/38: (7438.5/586.8, lift 1.8) */
        if (m5_lim_cloud_code == -1 && b2 > 3959 && b5 > 796 &&
            ndvi <= 0.0799109 && ndsi <= 0.85271)
            m5_lim_cloud_code = 100;  /* class cloud  [0.921] */
        
        /* Rule 4/39: (863.1/107.7, lift 1.8) */
        if (m5_lim_cloud_code == -1 && b1 > 2247 && b3 <= 2569 &&
            ndsi > -0.177533 && ndsi <= -0.115822)
            m5_lim_cloud_code = 100;  /* class cloud  [0.874] */
        
        /* Rule 4/40: (198.3/30.9, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b1 > 1511 && b2 > 1317 && b2 <= 1387 &&
            b7 > 1461 && ndvi > 0.175138)
            m5_lim_cloud_code = 100;  /* class cloud  [0.841] */
        
        /* Rule 4/41: (2902.8/463.9, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b1 > 2601 && b2 <= 2885 && b3 <= 2909 &&
            ndsi <= 0.143982)
            m5_lim_cloud_code = 100;  /* class cloud  [0.840] */
        
        /* Rule 4/42: (1471.7/245.2, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b1 > 1802 && b2 <= 1584 && b5 > 796 &&
            ndsi > -0.177533 && ndsi <= 0.267478)
            m5_lim_cloud_code = 100;  /* class cloud  [0.833] */
        
        /* Rule 4/43: (369.8/61.3, lift 1.7) */
        if (m5_lim_cloud_code == -1 && b1 > 2040 && b1 <= 2601 && b3 <= 2569 &&
            b5 <= 1493 && b7 > 1082 && ndsi <= 0.267478)
            m5_lim_cloud_code = 100;  /* class cloud  [0.832] */
        
        /* Rule 4/44: (3386.7/901.7, lift 1.5) */
        if (m5_lim_cloud_code == -1 && b1 > 1511 && b7 > 1461 &&
            ndvi > 0.175138 && ndsi > -0.22841)
            m5_lim_cloud_code = 100;  /* class cloud  [0.734] */
        
        /* Rule 4/45: (44465.5/20067, lift 1.1) */
        if (m5_lim_cloud_code == -1 && b5 > 796 && ndsi <= 0.85271)
            m5_lim_cloud_code = 100;  /* class cloud  [0.549] */ 
        
        if (m5_lim_cloud_code == -1)
            m5_lim_cloud_code = 100;  /* Default class: cloud_free */

        /* Use the maximum cloud code as the cloud value from this set of
           model runs */
        limited_cloud_code = m1_lim_cloud_code;
        if (m2_lim_cloud_code > limited_cloud_code)
            limited_cloud_code = m2_lim_cloud_code;
        if (m3_lim_cloud_code > limited_cloud_code)
            limited_cloud_code = m3_lim_cloud_code;
        if (m4_lim_cloud_code > limited_cloud_code)
            limited_cloud_code = m4_lim_cloud_code;
        if (m5_lim_cloud_code > limited_cloud_code)
            limited_cloud_code = m5_lim_cloud_code;

        /* Use the conservative and limited cloud codes to determine the
           revised cloud codes */
        if (conserv_cloud_code == 50)
            /* Original fmask identified cloud, but upon further inspection
               it looks like it's not cloudy */
            rev_cloud_mask[samp] = 0;
        else
            /* Original fmask identified cloud, and it appears to have been
               correct */
            rev_cloud_mask[samp] = 4;

        if (limited_cloud_code == 50)
            /* Original fmask identified cloud, but upon further inspection it
               looks like it's not cloudy */
            rev_lim_cloud_mask[samp] = 0;
        else
            /* Original fmask identified cloud, and it appears to have been
               correct */
            rev_lim_cloud_mask[samp] = 4;
    }  /* for samp */
}

