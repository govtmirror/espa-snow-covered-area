#include <math.h>
#include <string.h>
#include "revised_cloud_mask.h"

/******************************************************************************
MODULE:  rule_based_model

PURPOSE:  Runs the rule-based model on the input line of data, including the
reflectance bands, NDVI, and NDSI.  Five model runs are made per model.

RETURN VALUE:
Type = none

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
---------     ---------------  -------------------------------------
5/21/2014     Gail Schmidt     Original Development
12/2/2015     Gail Schmidt     Modified the output value for the cloud mask
12/2/2015     Gail Schmidt     Removed the variance-based model run

NOTES:
  1. Input and output arrays are 1D arrays of size nsamps.
  2. This algorithm was provided by David Selkowitz, USGS Alaska Science Center.
  3. The algorithm uses the scaled reflectance values as-is.
  4. The algorithm expects the NDVI and NDSI to be unscaled.
******************************************************************************/
void rule_based_model
(
    Input_t *input_img,     /* I: pointer to input data structure containing
                                  the scaled reflectance and cloud mask
                                  buffers (reflectance bands are scaled) */
    int nsamps,             /* I: number of samples in the input arrays */
    uint8 *rev_cloud_mask   /* O: revised cloud mask */
)
{
    int samp;             /* current sample being processed */
    int16 m1_cloud_code;  /* cloud mask for the 1st Rule-based model */
    int16 m2_cloud_code;  /* cloud mask for the 2nd Rule-based model */
    int16 m3_cloud_code;  /* cloud mask for the 3rd Rule-based model */
    int16 m4_cloud_code;  /* cloud mask for the 4th Rule-based model */
    int16 m5_cloud_code;  /* cloud mask for the 5th Rule-based model */
    int16 limited_cloud_code; /* maximum of the model runs for the
                                 limited (non-variance) cloud code */
    int16 b1, b2, b3, b4, b5, b7;  /* pixel values for the reflectance bands */
    double ndvi, ndsi;    /* pixel values for the NDVI and NDSI bands */

    /* Initialize the revised cloud mask to all zeros (no cloud) */
    memset (rev_cloud_mask, 0, nsamps * sizeof (uint8));

    /* Loop through the pixels in the array and run the Rule-based models.
       Any pixel which is not flagged as cloudy in the input cfmask will be
       skipped. */
    for (samp = 0; samp < nsamps; samp++)
    {
        /* If this isn't a cloudy pixel in the cfmask then skip to the next
           pixel */
        if (input_img->cfmask_buf[samp] != CFMASK_CLOUD)
            continue;

        /* Initialize the cloud mask */
        m1_cloud_code = -1;
        m2_cloud_code = -1;
        m3_cloud_code = -1;
        m4_cloud_code = -1;
        m5_cloud_code = -1;

        /* Set up the individual pixel values. Don't unscale the reflectance
           values or the reflectance variance values. */
        b1 = input_img->refl_buf[0][samp];
        b2 = input_img->refl_buf[1][samp];
        b3 = input_img->refl_buf[2][samp];
        b4 = input_img->refl_buf[3][samp];
        b5 = input_img->refl_buf[4][samp];
        b7 = input_img->refl_buf[5][samp];
        ndvi = make_index (b4, b3, input_img->refl_fill,
            input_img->refl_saturate_val);
        ndsi = make_index (b2, b5, input_img->refl_fill,
            input_img->refl_saturate_val);

        /* Initialize the cloud mask */
        m1_cloud_code = -1;
        m2_cloud_code = -1;
        m3_cloud_code = -1;
        m4_cloud_code = -1;
        m5_cloud_code = -1;

        /* Implement rule-based model for cases where variance calculation is
           not possible, such as SLC-off areas and areas near the edge of the
           image */
        /* 1st model run */
        /* Rule 0/1: (10764/235, lift 2.1) */
        if (m1_cloud_code == -1 && b3 > 1424 && b7 <= 1187)
            m1_cloud_code = 50;  /* class cloud_free  [0.978] */
        
        /* Rule 0/2: (753/16, lift 2.1) */
        if (m1_cloud_code == -1 && b1 <= 2999 && ndvi <= 0.158153 &&
            ndsi <= -0.212327)
            m1_cloud_code = 50;  /* class cloud_free  [0.977] */
        
        /* Rule 0/3: (447/10, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 2079 && b3 > 2080)
            m1_cloud_code = 50;  /* class cloud_free  [0.976] */
        
        /* Rule 0/4: (1763/48, lift 2.0) */
        if (m1_cloud_code == -1 && b1 > 6251 && b4 <= 5012 && b7 <= 2077)
            m1_cloud_code = 50;  /* class cloud_free  [0.972] */
        
        /* Rule 0/5: (960/27, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 2483 && b3 > 2408 &&
            ndsi > -0.212327)
            m1_cloud_code = 50;  /* class cloud_free  [0.971] */
        
        /* Rule 0/6: (728/29, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 2999 && b2 > 2771 &&
            ndvi <= 0.0912175)
            m1_cloud_code = 50;  /* class cloud_free  [0.959] */
        
        /* Rule 0/7: (4653/195, lift 2.0) */
        if (m1_cloud_code == -1 && b4 > 4172 && b7 <= 1464)
            m1_cloud_code = 50;  /* class cloud_free  [0.958] */
        
        /* Rule 0/8: (5826/279, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 3532 && b3 > 1727 && b4 > 1928 &&
            b7 <= 1464)
            m1_cloud_code = 50;  /* class cloud_free  [0.952] */
        
        /* Rule 0/9: (1917/93, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 2999 && b2 > 1970 && b5 <= 1937 &&
            b7 <= 1619 && ndvi <= 0.0912175)
            m1_cloud_code = 50;  /* class cloud_free  [0.951] */
        
        /* Rule 0/10: (1198/63, lift 2.0) */
        if (m1_cloud_code == -1 && b1 > 2999 && b1 <= 3589 && b4 > 3588 &&
            b7 <= 2077)
            m1_cloud_code = 50;  /* class cloud_free  [0.947] */
        
        /* Rule 0/11: (1293/71, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 2999 && b4 > 3173 &&
            ndsi > -0.0444104)
            m1_cloud_code = 50;  /* class cloud_free  [0.944] */
        
        /* Rule 0/12: (2832/157, lift 2.0) */
        if (m1_cloud_code == -1 && b1 <= 1600 && b4 <= 1928)
            m1_cloud_code = 50;  /* class cloud_free  [0.944] */
        
        /* Rule 0/13: (2380/185, lift 1.9) */
        if (m1_cloud_code == -1 && b4 <= 2729 && ndsi <= -0.212327)
            m1_cloud_code = 50;  /* class cloud_free  [0.922] */
        
        /* Rule 0/14: (2208/176, lift 1.9) */
        if (m1_cloud_code == -1 && b1 <= 2638 && b5 > 1937 &&
            ndvi <= 0.0912175)
            m1_cloud_code = 50;  /* class cloud_free  [0.920] */
        
        /* Rule 0/15: (1742/156, lift 1.9) */
        if (m1_cloud_code == -1 && b7 <= 1799 && ndsi <= -0.212327)
            m1_cloud_code = 50;  /* class cloud_free  [0.910] */
        
        /* Rule 0/16: (1724/180, lift 1.9) */
        if (m1_cloud_code == -1 && b1 <= 1993 && b5 > 1937 && b7 > 1464 &&
            ndvi <= 0.210062)
            m1_cloud_code = 50;  /* class cloud_free  [0.895] */
        
        /* Rule 0/17: (23277/7218, lift 1.4) */
        if (m1_cloud_code == -1 && b1 <= 2999)
            m1_cloud_code = 50;  /* class cloud_free  [0.690] */
        
        /* Rule 0/18: (3593/38, lift 1.9) */
        if (m1_cloud_code == -1 && b1 > 2483 && ndvi > 0.0912175 &&
            ndsi <= -0.0444104)
            m1_cloud_code = 100;  /* class cloud  [0.989] */
        
        /* Rule 0/19: (5613/71, lift 1.9) */
        if (m1_cloud_code == -1 && b7 > 1266 && ndsi > 0.491986)
            m1_cloud_code = 100;  /* class cloud  [0.987] */
        
        /* Rule 0/20: (1100/125, lift 1.7) */
        if (m1_cloud_code == -1 && b7 > 1464 && ndvi > 0.210062 &&
            ndsi > -0.212327)
            m1_cloud_code = 100;  /* class cloud  [0.886] */
        
        /* Rule 0/21: (267/39, lift 1.6) */
        if (m1_cloud_code == -1 && b1 > 1600 && b3 <= 1299 && b4 <= 1928 &&
            b7 > 919 && b7 <= 1464)
            m1_cloud_code = 100;  /* class cloud  [0.851] */
        
        /* Rule 0/22: (193/35, lift 1.6) */
        if (m1_cloud_code == -1 && b1 > 1451 && b3 <= 1139 && b7 <= 919 &&
            ndvi > 0.0602801 && ndvi <= 0.421493 && ndsi <= 0.169811)
            m1_cloud_code = 100;  /* class cloud  [0.815] */
        
        /* Rule 0/23: (32298/7826, lift 1.4) */
        if (m1_cloud_code == -1 && b1 > 1802 && b7 > 919)
            m1_cloud_code = 100;  /* class cloud  [0.758] */
        
        if (m1_cloud_code == -1)
            m1_cloud_code = 100;


        /* 2nd model run */
        /* Rule 1/1: (801.4/3.8, lift 2.1) */
        if (m2_cloud_code == -1 && b4 > 6829 && b5 <= 3315)
            m2_cloud_code = 50;  /* class cloud_free  [0.994] */
        
        /* Rule 1/2: (657/4.5, lift 2.1) */
        if (m2_cloud_code == -1 && ndvi > 0.44504)
            m2_cloud_code = 50;  /* class cloud_free  [0.992] */
        
        /* Rule 1/3: (220.2/1.5, lift 2.1) */
        if (m2_cloud_code == -1 && b1 <= 3544 && b4 > 4212 && b5 > 1136 &&
            b5 <= 3315)
            m2_cloud_code = 50;  /* class cloud_free  [0.989] */
        
        /* Rule 1/4: (6122.8/130.4, lift 2.1) */
        if (m2_cloud_code == -1 && b3 > 1258 && b5 <= 1136)
            m2_cloud_code = 50;  /* class cloud_free  [0.979] */
        
        /* Rule 1/5: (988.6/45.5, lift 2.1) */
        if (m2_cloud_code == -1 && b1 <= 3544 && b5 > 1136 && b7 <= 1067 &&
            ndsi > 0.143936)
            m2_cloud_code = 50;  /* class cloud_free  [0.953] */
        
        /* Rule 1/6: (3721.5/193.8, lift 2.0) */
        if (m2_cloud_code == -1 && b1 > 6251 && b7 <= 1683)
            m2_cloud_code = 50;  /* class cloud_free  [0.948] */
        
        /* Rule 1/7: (504.2/47.8, lift 1.9) */
        if (m2_cloud_code == -1 && b1 > 1555 && b1 <= 2321 && b4 <= 3520 &&
            b5 > 2467 && b7 <= 1726)
            m2_cloud_code = 50;  /* class cloud_free  [0.904] */
        
        /* Rule 1/8: (1057/103.5, lift 1.9) */
        if (m2_cloud_code == -1 && b1 <= 3544 && b2 > 3236 && b5 <= 3315)
            m2_cloud_code = 50;  /* class cloud_free  [0.901] */
        
        /* Rule 1/9: (671.1/76.4, lift 1.9) */
        if (m2_cloud_code == -1 && b1 <= 2835 && b3 > 2538 &&
            ndsi > -0.106196 && ndsi <= 0.143936)
            m2_cloud_code = 50;  /* class cloud_free  [0.885] */
        
        /* Rule 1/10: (650.1/75.1, lift 1.9) */
        if (m2_cloud_code == -1 && b4 <= 4001 && b5 > 3315 &&
            ndvi <= 0.0617448)
            m2_cloud_code = 50;  /* class cloud_free  [0.883] */
        
        /* Rule 1/11: (1812.6/334.2, lift 1.8) */
        if (m2_cloud_code == -1 && b1 <= 2082 && b3 > 1839 && b4 <= 3520)
            m2_cloud_code = 50;  /* class cloud_free  [0.815] */
        
        /* Rule 1/12: (31120.8/12936.8, lift 1.3) */
        if (m2_cloud_code == -1 && b4 <= 4001)
            m2_cloud_code = 50;  /* class cloud_free  [0.584] */
        
        /* Rule 1/13: (11094.4/118.4, lift 1.8) */
        if (m2_cloud_code == -1 && b4 > 4001 && b5 > 3315)
            m2_cloud_code = 100;  /* class cloud  [0.989] */
        
        /* Rule 1/14: (519/12.9, lift 1.8) */
        if (m2_cloud_code == -1 && b1 > 2043 && b2 <= 1763 &&
            ndsi <= 0.143936)
            m2_cloud_code = 100;  /* class cloud  [0.973] */
        
        /* Rule 1/15: (12178.2/334.3, lift 1.8) */
        if (m2_cloud_code == -1 && b1 > 3544 && b7 > 1683)
            m2_cloud_code = 100;  /* class cloud  [0.972] */
        
        /* Rule 1/16: (3402.2/106.4, lift 1.8) */
        if (m2_cloud_code == -1 && b1 > 2599 && b5 > 3315 &&
            ndvi > 0.0617448)
            m2_cloud_code = 100;  /* class cloud  [0.968] */
        
        /* Rule 1/17: (9675.1/785.7, lift 1.7) */
        if (m2_cloud_code == -1 && b1 > 2835 && ndsi <= 0.143936)
            m2_cloud_code = 100;  /* class cloud  [0.919] */
        
        /* Rule 1/18: (1278.7/104.9, lift 1.7) */
        if (m2_cloud_code == -1 && b1 > 3264 && b2 <= 3334 && b7 > 1067)
            m2_cloud_code = 100;  /* class cloud  [0.917] */
        
        /* Rule 1/19: (564.8/53.3, lift 1.7) */
        if (m2_cloud_code == -1 && b1 > 2151 && b2 <= 1898 && b5 > 1537)
            m2_cloud_code = 100;  /* class cloud  [0.904] */
        
        /* Rule 1/20: (228.7/24.9, lift 1.7) */
        if (m2_cloud_code == -1 && b1 > 2082 && b3 <= 1994 &&
            ndvi > 0.179141 && ndsi <= 0.143936)
            m2_cloud_code = 100;  /* class cloud  [0.888] */
        
        /* Rule 1/21: (4430.9/528.7, lift 1.6) */
        if (m2_cloud_code == -1 && b1 > 3544 && b1 <= 6251 && b4 <= 6829 &&
            b5 > 1136)
            m2_cloud_code = 100;  /* class cloud  [0.881] */
        
        /* Rule 1/22: (179.6/32.6, lift 1.5) */
        if (m2_cloud_code == -1 && b2 <= 2204 && b4 > 3520 && b5 <= 3315 &&
            ndvi <= 0.44504 && ndsi <= 0.143936)
            m2_cloud_code = 100;  /* class cloud  [0.815] */
        
        /* Rule 1/23: (481.3/131.1, lift 1.4) */
        if (m2_cloud_code == -1 && b1 > 1800 && b1 <= 2082 && b3 <= 1839 &&
            b5 <= 3315 && b7 > 1342 && ndvi > 0.179141 && ndvi <= 0.207702)
            m2_cloud_code = 100;  /* class cloud  [0.727] */
        
        /* Rule 1/24: (4787.7/1392.3, lift 1.3) */
        if (m2_cloud_code == -1 && b1 > 2321 && b2 <= 3236 && b4 <= 4212 &&
            b5 <= 3315 && ndsi <= 0.143936)
            m2_cloud_code = 100;  /* class cloud  [0.709] */
        
        /* Rule 1/25: (2379.8/758.4, lift 1.3) */
        if (m2_cloud_code == -1 && b1 > 1555 && b1 <= 1800 && b3 <= 1839 &&
            b5 > 1537 && ndvi > 0.179141 && ndvi <= 0.44504 &&
            ndsi <= -0.0472716)
            m2_cloud_code = 100;  /* class cloud  [0.681] */
        
        /* Rule 1/26: (2563.3/868.5, lift 1.2) */
        if (m2_cloud_code == -1 && b1 > 1182 && b5 > 1136 && b5 <= 1537 &&
            b7 > 542 && ndvi <= 0.44504 && ndsi <= 0.143936)
            m2_cloud_code = 100;  /* class cloud  [0.661] */
        
        /* Rule 1/27: (897.2/366.5, lift 1.1) */
        if (m2_cloud_code == -1 && b3 <= 1258 && b5 <= 1136 && b7 > 542 &&
            ndvi <= 0.44504)
            m2_cloud_code = 100;  /* class cloud  [0.591] */
        
        if (m2_cloud_code == -1)
            m2_cloud_code = 50;  /* Default class: cloud_free */


        /* 3rd model run */
        /* Rule 2/1: (286.3, lift 2.0) */
        if (m3_cloud_code == -1 && b1 <= 2497 && b3 > 2156 && b3 <= 2276 &&
            b7 <= 1626 && ndvi > 0.0650248)
            m3_cloud_code = 50;  /* class cloud_free  [0.997] */
        
        /* Rule 2/2: (2084.7/6.1, lift 2.0) */
        if (m3_cloud_code == -1 && ndsi > 0.881556)
            m3_cloud_code = 50;  /* class cloud_free  [0.997] */
        
        /* Rule 2/3: (347/7.9, lift 2.0) */
        if (m3_cloud_code == -1 && b2 > 3653 && b7 <= 1626 &&
            ndvi > 0.0650248)
            m3_cloud_code = 50;  /* class cloud_free  [0.974] */
        
        /* Rule 2/4: (642.5/16.7, lift 2.0) */
        if (m3_cloud_code == -1 && b1 <= 2497 && b3 > 1979 && b7 <= 1626 &&
            ndsi <= 0.0546875)
            m3_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 2/5: (3129.8/84.8, lift 2.0) */
        if (m3_cloud_code == -1 && b2 > 1011 && b7 <= 486)
            m3_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 2/6: (264.7/8.4, lift 2.0) */
        if (m3_cloud_code == -1 && b1 <= 1966 && b5 > 1707 && b7 <= 1626 &&
            ndvi > 0.0650248 && ndvi <= 0.12763)
            m3_cloud_code = 50;  /* class cloud_free  [0.965] */
        
        /* Rule 2/7: (152/6.1, lift 1.9) */
        if (m3_cloud_code == -1 && b5 > 2507 && b7 <= 1626 &&
            ndvi > 0.0650248 && ndsi > -0.254112)
            m3_cloud_code = 50;  /* class cloud_free  [0.954] */
        
        /* Rule 2/8: (1358.7/66.7, lift 1.9) */
        if (m3_cloud_code == -1 && b2 <= 1011)
            m3_cloud_code = 50;  /* class cloud_free  [0.950] */
        
        /* Rule 2/9: (7074/402.6, lift 1.9) */
        if (m3_cloud_code == -1 && b3 > 1265 && b7 <= 1008)
            m3_cloud_code = 50;  /* class cloud_free  [0.943] */
        
        /* Rule 2/10: (918.4/59.4, lift 1.9) */
        if (m3_cloud_code == -1 && b3 <= 1265 && ndvi > 0.403095)
            m3_cloud_code = 50;  /* class cloud_free  [0.934] */
        
        /* Rule 2/11: (2111.8/153.8, lift 1.9) */
        if (m3_cloud_code == -1 && b1 <= 1725 && b7 <= 1286 &&
            ndvi <= 0.139127)
            m3_cloud_code = 50;  /* class cloud_free  [0.927] */
        
        /* Rule 2/12: (1275.7/107.6, lift 1.9) */
        if (m3_cloud_code == -1 && b1 <= 1725 && b3 > 1205 && b4 <= 2458 &&
            b7 <= 1286)
            m3_cloud_code = 50;  /* class cloud_free  [0.915] */
        
        /* Rule 2/13: (3639.5/324.7, lift 1.8) */
        if (m3_cloud_code == -1 && b1 > 1842 && b2 <= 4235 && b4 > 1742 &&
            b7 <= 1286 && ndvi <= 0.0979757)
            m3_cloud_code = 50;  /* class cloud_free  [0.911] */
        
        /* Rule 2/14: (1170.6/111.3, lift 1.8) */
        if (m3_cloud_code == -1 && b2 <= 4235 && b3 > 1265 && b5 > 1747 &&
            b7 <= 1286)
            m3_cloud_code = 50;  /* class cloud_free  [0.904] */
        
        /* Rule 2/15: (126.8/12.2, lift 1.8) */
        if (m3_cloud_code == -1 && b4 > 6150 && b5 <= 3315 &&
            ndsi <= 0.49737)
            m3_cloud_code = 50;  /* class cloud_free  [0.898] */
        
        /* Rule 2/16: (1897.2/223.5, lift 1.8) */
        if (m3_cloud_code == -1 && b5 <= 3315 && ndsi > 0.368125 &&
            ndsi <= 0.49737)
            m3_cloud_code = 50;  /* class cloud_free  [0.882] */
        
        /* Rule 2/17: (775.9/102.5, lift 1.8) */
        if (m3_cloud_code == -1 && b1 <= 3284 && b3 > 2276 && b5 > 1707 &&
            b7 <= 1626 && ndvi > 0.0650248)
            m3_cloud_code = 50;  /* class cloud_free  [0.867] */
        
        /* Rule 2/18: (441.9/66.5, lift 1.7) */
        if (m3_cloud_code == -1 && b4 <= 3074 && b5 > 3315)
            m3_cloud_code = 50;  /* class cloud_free  [0.848] */
        
        /* Rule 2/19: (972/154.3, lift 1.7) */
        if (m3_cloud_code == -1 && b1 <= 3418 && b5 > 3315 &&
            ndvi <= 0.079329)
            m3_cloud_code = 50;  /* class cloud_free  [0.841] */
        
        /* Rule 2/20: (1000.4/161.7, lift 1.7) */
        if (m3_cloud_code == -1 && b3 <= 1205 && b7 <= 1286 &&
            ndsi <= -0.0934991)
            m3_cloud_code = 50;  /* class cloud_free  [0.838] */
        
        /* Rule 2/21: (168.1/31.5, lift 1.6) */
        if (m3_cloud_code == -1 && b5 > 2885 && b7 <= 2112 &&
            ndvi > 0.0650248 && ndsi > -0.254112)
            m3_cloud_code = 50;  /* class cloud_free  [0.809] */
        
        /* Rule 2/22: (3523.5/696, lift 1.6) */
        if (m3_cloud_code == -1 && b1 <= 1966 && b3 > 1220 && b4 <= 2470 &&
            b7 <= 1458)
            m3_cloud_code = 50;  /* class cloud_free  [0.802] */
        
        /* Rule 2/23: (795.3/160.6, lift 1.6) */
        if (m3_cloud_code == -1 && b1 <= 1966 && b4 > 2470 && b4 <= 3217 &&
            b7 > 1286 && b7 <= 1626)
            m3_cloud_code = 50;  /* class cloud_free  [0.797] */
        
        /* Rule 2/24: (557.4/120.5, lift 1.6) */
        if (m3_cloud_code == -1 && b5 > 2528 && b7 <= 1805 &&
            ndvi > 0.0650248 && ndsi > -0.254112)
            m3_cloud_code = 50;  /* class cloud_free  [0.783] */
        
        /* Rule 2/25: (8172.8/1851.2, lift 1.6) */
        if (m3_cloud_code == -1 && b1 <= 3376 && ndvi <= 0.0650248)
            m3_cloud_code = 50;  /* class cloud_free  [0.773] */
        
        /* Rule 2/26: (1800.2/430.1, lift 1.5) */
        if (m3_cloud_code == -1 && b4 <= 2856 && ndsi <= -0.254112)
            m3_cloud_code = 50;  /* class cloud_free  [0.761] */
        
        /* Rule 2/27: (308.2/41.3, lift 1.7) */
        if (m3_cloud_code == -1 && b1 > 1725 && b3 <= 1265 && b7 > 486)
            m3_cloud_code = 100;  /* class cloud  [0.864] */
        
        /* Rule 2/28: (278.4/47.5, lift 1.6) */
        if (m3_cloud_code == -1 && b1 <= 1725 && b2 > 1011 && b3 <= 1265 &&
            b4 > 2458 && b7 > 486 && b7 <= 1286 && ndvi <= 0.403095)
            m3_cloud_code = 100;  /* class cloud  [0.827] */
        
        /* Rule 2/29: (1011.9/222.6, lift 1.5) */
        if (m3_cloud_code == -1 && b2 > 1011 && b5 <= 1707 && b7 > 1286 &&
            ndsi <= 0.368125)
            m3_cloud_code = 100;  /* class cloud  [0.779] */
        
        /* Rule 2/30: (515.1/136.5, lift 1.5) */
        if (m3_cloud_code == -1 && b1 > 1842 && b4 <= 1742 && b5 <= 1747 &&
            b7 > 1008)
            m3_cloud_code = 100;  /* class cloud  [0.734] */
        
        /* Rule 2/31: (46915.3/22120.7, lift 1.0) */
        if (m3_cloud_code == -1 && ndsi <= 0.881556)
            m3_cloud_code = 100;  /* class cloud  [0.528] */
        
        if (m3_cloud_code == -1)
            m3_cloud_code = 100;  /* Default class: cloud */


        /* 4th model run */
        /* Rule 3/1: (709.4/15.7, lift 1.9) */
        if (m4_cloud_code == -1 && b1 <= 2118 && b2 > 2012)
            m4_cloud_code = 50;  /* class cloud_free  [0.976] */
        
        /* Rule 3/2: (196.5/4.4, lift 1.8) */
        if (m4_cloud_code == -1 && b1 <= 3234 && b2 > 3108 && b5 > 3315)
            m4_cloud_code = 50;  /* class cloud_free  [0.973] */
        
        /* Rule 3/3: (1753.6/96.9, lift 1.8) */
        if (m4_cloud_code == -1 && b1 <= 2467 && b3 > 2364)
            m4_cloud_code = 50;  /* class cloud_free  [0.944] */
        
        /* Rule 3/4: (1047.4/69.7, lift 1.8) */
        if (m4_cloud_code == -1 && b1 <= 2654 && b5 > 3315 &&
            ndvi <= 0.1719)
            m4_cloud_code = 50;  /* class cloud_free  [0.933] */
        
        /* Rule 3/5: (2283.3/179.7, lift 1.7) */
        if (m4_cloud_code == -1 && b1 <= 2467 && b2 > 1983 && b5 > 1706 &&
            ndvi <= 0.114312)
            m4_cloud_code = 50;  /* class cloud_free  [0.921] */
        
        /* Rule 3/6: (1552.5/128.6, lift 1.7) */
        if (m4_cloud_code == -1 && b1 <= 2762 && b3 > 2437 &&
            ndsi > -0.115228)
            m4_cloud_code = 50;  /* class cloud_free  [0.917] */
        
        /* Rule 3/7: (3784.7/387.6, lift 1.7) */
        if (m4_cloud_code == -1 && b3 > 1665 && b7 <= 1347 &&
            ndvi > -0.0240739 && ndvi <= 0.131339)
            m4_cloud_code = 50;  /* class cloud_free  [0.897] */
        
        /* Rule 3/8: (3288.9/361.2, lift 1.7) */
        if (m4_cloud_code == -1 && b1 > 6251 && b5 <= 3112 && b7 <= 2077)
            m4_cloud_code = 50;  /* class cloud_free  [0.890] */
        
        /* Rule 3/9: (1266.7/159.1, lift 1.7) */
        if (m4_cloud_code == -1 && b1 <= 3160 && b3 > 2834 && b5 <= 3112)
            m4_cloud_code = 50;  /* class cloud_free  [0.874] */
        
        /* Rule 3/10: (1701.1/238.9, lift 1.6) */
        if (m4_cloud_code == -1 && b4 > 5404 && b5 <= 3112 && b7 <= 2077)
            m4_cloud_code = 50;  /* class cloud_free  [0.859] */
        
        /* Rule 3/11: (388.7/55.3, lift 1.6) */
        if (m4_cloud_code == -1 && b1 <= 2467 && b5 > 1706 && b7 > 1347 &&
            ndvi > 0.114312 && ndsi > -0.0444104)
            m4_cloud_code = 50;  /* class cloud_free  [0.856] */
        
        /* Rule 3/12: (902.2/145.9, lift 1.6) */
        if (m4_cloud_code == -1 && b1 <= 2118 && b4 <= 3249 && b5 > 2406 &&
            b7 <= 1726)
            m4_cloud_code = 50;  /* class cloud_free  [0.838] */
        
        /* Rule 3/13: (2501.9/446.9, lift 1.6) */
        if (m4_cloud_code == -1 && b1 <= 2118 && b4 <= 2223 && b5 > 1706)
            m4_cloud_code = 50;  /* class cloud_free  [0.821] */
        
        /* Rule 3/14: (2864.7/535.5, lift 1.5) */
        if (m4_cloud_code == -1 && b1 <= 1873 && b5 > 1138 && b7 <= 1347 &&
            ndvi <= 0.317016)
            m4_cloud_code = 50;  /* class cloud_free  [0.813] */
        
        /* Rule 3/15: (1899.1/420.7, lift 1.5) */
        if (m4_cloud_code == -1 && b1 <= 3544 && b2 > 3043 && b5 <= 3112)
            m4_cloud_code = 50;  /* class cloud_free  [0.778] */
        
        /* Rule 3/16: (1598/394.3, lift 1.4) */
        if (m4_cloud_code == -1 && b1 <= 2654 && b5 > 3315)
            m4_cloud_code = 50;  /* class cloud_free  [0.753] */
        
        /* Rule 3/17: (737.3/186, lift 1.4) */
        if (m4_cloud_code == -1 && b1 <= 1684 && b2 > 1342 && b5 <= 2406 &&
            b7 > 1347)
            m4_cloud_code = 50;  /* class cloud_free  [0.747] */
        
        /* Rule 3/18: (38683.4/14494.6, lift 1.2) */
        if (m4_cloud_code == -1 && b5 <= 3315)
            m4_cloud_code = 50;  /* class cloud_free  [0.625] */
        
        /* Rule 3/19: (7194.6, lift 2.1) */
        if (m4_cloud_code == -1 && b1 > 3544 && b7 > 2077)
            m4_cloud_code = 100;   /* class cloud  [1.000] */
        
        /* Rule 3/20: (323/3.8, lift 2.1) */
        if (m4_cloud_code == -1 && b1 > 2467 && b5 <= 3315 &&
            ndsi <= -0.115228)
            m4_cloud_code = 100;  /* class cloud  [0.985] */
        
        /* Rule 3/21: (563.4/23.1, lift 2.0) */
        if (m4_cloud_code == -1 && b1 > 3160 && b2 <= 3043 && b7 > 1347)
            m4_cloud_code = 100;  /* class cloud  [0.957] */
        
        /* Rule 3/22: (8425.6/406.2, lift 2.0) */
        if (m4_cloud_code == -1 && b1 > 2654 && b5 > 3315)
            m4_cloud_code = 100;  /* class cloud  [0.952] */
        
        /* Rule 3/23: (8942.7/654, lift 2.0) */
        if (m4_cloud_code == -1 && b1 > 3544 && b7 > 1347 &&
            ndvi > -0.626148)
            m4_cloud_code = 100;  /* class cloud  [0.927] */
        
        /* Rule 3/24: (5373.8/411.2, lift 2.0) */
        if (m4_cloud_code == -1 && b5 > 1157 && ndvi > -0.564211 &&
            ndvi <= -0.0240739 && ndsi <= 0.881468)
            m4_cloud_code = 100;  /* class cloud  [0.923] */
        
        /* Rule 3/25: (1531.5/144.8, lift 1.9) */
        if (m4_cloud_code == -1 && b1 > 2762 && b3 <= 2834 && b7 > 1347)
            m4_cloud_code = 100;  /* class cloud  [0.905] */
        
        /* Rule 3/26: (115.6/13, lift 1.9) */
        if (m4_cloud_code == -1 && b1 > 1526 && b2 <= 1342 && b5 <= 3315 &&
            b7 > 1347 && ndvi > 0.114312)
            m4_cloud_code = 100;  /* class cloud  [0.881] */
        
        /* Rule 3/27: (1386.2/198.3, lift 1.8) */
        if (m4_cloud_code == -1 && b1 > 2467 && b3 <= 2437 && b7 > 1347)
            m4_cloud_code = 100;  /* class cloud  [0.856] */
        
        /* Rule 3/28: (107.2/23.2, lift 1.7) */
        if (m4_cloud_code == -1 && b5 > 1136 && b5 <= 1138)
            m4_cloud_code = 100;  /* class cloud  [0.778] */
        
        /* Rule 3/29: (1097.6/247.2, lift 1.7) */
        if (m4_cloud_code == -1 && b1 > 2123 && b2 <= 1983 && b7 > 1347)
            m4_cloud_code = 100;  /* class cloud  [0.774] */
        
        /* Rule 3/30: (4230.3/989.1, lift 1.6) */
        if (m4_cloud_code == -1 && b1 > 2118 && b7 > 1347 &&
            ndvi > 0.114312 && ndsi <= -0.0444104)
            m4_cloud_code = 100;  /* class cloud  [0.766] */
        
        /* Rule 3/31: (483.4/114.7, lift 1.6) */
        if (m4_cloud_code == -1 && b1 > 1526 && b2 <= 2012 && b4 > 3249 &&
            b5 <= 3315 && b7 > 1347 && ndvi <= 0.44504)
            m4_cloud_code = 100;  /* class cloud  [0.762] */
        
        /* Rule 3/32: (516.7/123.6, lift 1.6) */
        if (m4_cloud_code == -1 && b1 > 1526 && b1 <= 2467 && b3 <= 2364 &&
            b5 <= 1706 && b7 > 1347)
            m4_cloud_code = 100;  /* class cloud  [0.760] */
        
        /* Rule 3/33: (1601.7/412.5, lift 1.6) */
        if (m4_cloud_code == -1 && b1 > 1873 && b3 <= 1665 && b5 > 1157 &&
            ndvi > -0.0240739 && ndvi <= 0.317016 && ndsi > -0.184974)
            m4_cloud_code = 100;  /* class cloud  [0.742] */
        
        /* Rule 3/34: (2091.4/633.1, lift 1.5) */
        if (m4_cloud_code == -1 && b1 > 1684 && b4 > 2223 && b5 <= 2406 &&
            b7 > 1347 && ndvi > 0.114312 && ndvi <= 0.44504 &&
            ndsi <= -0.0800831)
            m4_cloud_code = 100;  /* class cloud  [0.697] */
        
        /* Rule 3/35: (2531.2/772, lift 1.5) */
        if (m4_cloud_code == -1 && b1 > 1526 && b3 <= 2364 && b5 <= 2626 &&
            b7 > 1726 && ndvi > 0.114312)
            m4_cloud_code = 100;  /* class cloud  [0.695] */
        
        /* Rule 3/36: (781.6/292.1, lift 1.3) */
        if (m4_cloud_code == -1 && b7 > 1086 && ndvi > 0.131339 &&
            ndsi > 0.0375777 && ndsi <= 0.416984)
            m4_cloud_code = 100;  /* class cloud  [0.626] */
        
        /* Rule 3/37: (4079.4/1565.3, lift 1.3) */
        if (m4_cloud_code == -1 && b1 > 1526 && b3 <= 2364 && b5 <= 2406 &&
            b7 > 1347 && ndvi > 0.114312 && ndvi <= 0.44504)
            m4_cloud_code = 100;  /* class cloud  [0.616] */

        if (m4_cloud_code == -1)
            m4_cloud_code = 50;  /* Default class: cloud_free */

        
        /* 5th model run */
        /* Rule 4/1: (400.1, lift 2.1) */
        if (m5_cloud_code == -1 && b7 <= 1461 && ndsi <= -0.224928)
            m5_cloud_code = 50;  /* class cloud_free  [0.998] */
        
        /* Rule 4/2: (815.6/0.4, lift 2.1) */
        if (m5_cloud_code == -1 && b1 <= 2601 && b3 > 2260 && b7 <= 1623)
            m5_cloud_code = 50;  /* class cloud_free  [0.998] */
        
        /* Rule 4/3: (273, lift 2.1) */
        if (m5_cloud_code == -1 && b1 <= 1511 && b2 > 1317 && b5 > 1516)
            m5_cloud_code = 50;  /* class cloud_free  [0.996] */
        
        /* Rule 4/4: (1665.9/5.6, lift 2.1) */
        if (m5_cloud_code == -1 && ndsi > 0.85271)
            m5_cloud_code = 50;  /* class cloud_free  [0.996] */
        
        /* Rule 4/5: (677.9/4.5, lift 2.1) */
        if (m5_cloud_code == -1 && b1 <= 2601 && ndvi <= 0.107894 &&
            ndsi <= -0.177533)
            m5_cloud_code = 50;  /* class cloud_free  [0.992] */
        
        /* Rule 4/6: (1275.5/13, lift 2.1) */
        if (m5_cloud_code == -1 && b1 <= 1645 && ndvi <= 0.107894)
            m5_cloud_code = 50;  /* class cloud_free  [0.989] */
        
        /* Rule 4/7: (586.3/6.4, lift 2.1) */
        if (m5_cloud_code == -1 && b2 <= 1027 && b5 > 796)
            m5_cloud_code = 50;  /* class cloud_free  [0.987] */
        
        /* Rule 4/8: (2953.5/46.3, lift 2.0) */
        if (m5_cloud_code == -1 && b5 <= 796)
            m5_cloud_code = 50;   /* class cloud_free  [0.984] */
        
        /* Rule 4/9: (1042.7/19.5, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 2601 && b3 > 2569)
            m5_cloud_code = 50;  /* class cloud_free  [0.980] */
        
        /* Rule 4/10: (988.3/19.6, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 1802 && b2 > 1292 &&
            ndvi <= 0.107894)
            m5_cloud_code = 50;  /* class cloud_free  [0.979] */
        
        /* Rule 4/11: (268/4.7, lift 2.0) */
        if (m5_cloud_code == -1 && b2 > 1387 && b3 <= 1461 &&
            ndvi <= 0.412224 && ndsi <= -0.22841)
            m5_cloud_code = 50;  /* class cloud_free  [0.979] */
        
        /* Rule 4/12: (719.9/16, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 2005 && b3 > 1949)
            m5_cloud_code = 50;  /* class cloud_free  [0.976] */
        
        /* Rule 4/13: (565.2/15.4, lift 2.0) */
        if (m5_cloud_code == -1 && ndvi > 0.412224)
            m5_cloud_code = 50;  /* class cloud_free  [0.971] */
        
        /* Rule 4/14: (3293.8/99.1, lift 2.0) */
        if (m5_cloud_code == -1 && b7 <= 579)
            m5_cloud_code = 50;  /* class cloud_free  [0.970] */
        
        /* Rule 4/15: (966.4/33.4, lift 2.0) */
        if (m5_cloud_code == -1 && b1 <= 3037 && b2 > 2885 && b5 <= 4374)
            m5_cloud_code = 50;  /* class cloud_free  [0.964] */
        
        /* Rule 4/16: (744.4/28.7, lift 2.0) */
        if (m5_cloud_code == -1 && b2 <= 1317 && b5 > 1516 &&
            ndvi > 0.107894)
            m5_cloud_code = 50;  /* class cloud_free  [0.960] */
        
        /* Rule 4/17: (864.7/38.2, lift 2.0) */
        if (m5_cloud_code == -1 && b1 > 2601 && b1 <= 3037 &&
            ndsi > 0.143982)
            m5_cloud_code = 50;  /* class cloud_free  [0.955] */
        
        /* Rule 4/18: (1667.8/114.9, lift 1.9) */
        if (m5_cloud_code == -1 && b1 <= 3037 && b3 > 2909)
            m5_cloud_code = 50;  /* class cloud_free  [0.931] */
        
        /* Rule 4/19: (441.2/30.5, lift 1.9) */
        if (m5_cloud_code == -1 && b4 > 7015 && b5 <= 4374)
            m5_cloud_code = 50;  /* class cloud_free  [0.929] */
        
        /* Rule 4/20: (129.7/8.6, lift 1.9) */
        if (m5_cloud_code == -1 && b4 <= 4251 && b5 > 4374)
            m5_cloud_code = 50;  /* class cloud_free  [0.927] */
        
        /* Rule 4/21: (1596/121.5, lift 1.9) */
        if (m5_cloud_code == -1 && b1 <= 2247 && b2 > 1763 && b5 > 1493 &&
            ndvi <= 0.107894)
            m5_cloud_code = 50;  /* class cloud_free  [0.923] */
        
        /* Rule 4/22: (6175.3/491.4, lift 1.9) */
        if (m5_cloud_code == -1 && b3 > 1223 && b7 <= 1084 &&
            ndsi > -0.0853392)
            m5_cloud_code = 50;  /* class cloud_free  [0.920] */
        
        /* Rule 4/23: (3766.2/366.4, lift 1.9) */
        if (m5_cloud_code == -1 && b2 <= 3959 && ndsi > 0.267478)
            m5_cloud_code = 50;  /* class cloud_free  [0.902] */
        
        /* Rule 4/24: (414.9/46.3, lift 1.8) */
        if (m5_cloud_code == -1 && b5 <= 1516 && ndsi <= -0.0853392)
            m5_cloud_code = 50;  /* class cloud_free  [0.887] */
        
        /* Rule 4/25: (3150.8/363.8, lift 1.8) */
        if (m5_cloud_code == -1 && b1 <= 2005 && b5 > 1516 &&
            ndvi <= 0.175138)
            m5_cloud_code = 50;  /* class cloud_free  [0.884] */
        
        /* Rule 4/26: (1868.5/321.3, lift 1.7) */
        if (m5_cloud_code == -1 && b1 <= 2005 && b4 <= 2856 &&
            ndsi <= -0.22841)
            m5_cloud_code = 50;  /* class cloud_free  [0.828] */
        
        /* Rule 4/27: (126.6/21.5, lift 1.7) */
        if (m5_cloud_code == -1 && b2 > 3959 && b5 <= 4374 &&
            ndvi > 0.0799109)
            m5_cloud_code = 50;  /* class cloud_free  [0.825] */
        
        /* Rule 4/28: (2695.3/480.2, lift 1.7) */
        if (m5_cloud_code == -1 && b1 <= 2601 && b2 > 1906 &&
            ndvi <= 0.107894 && ndsi > -0.115822)
            m5_cloud_code = 50;  /* class cloud_free  [0.822] */
        
        /* Rule 4/29: (3038/650.2, lift 1.6) */
        if (m5_cloud_code == -1 && b1 <= 2394 && b5 > 2533 &&
            ndvi <= 0.201976)
            m5_cloud_code = 50;  /* class cloud_free  [0.786] */
        
        /* Rule 4/30: (13065.9/2877.5, lift 1.6) */
        if (m5_cloud_code == -1 && b3 > 1267 && b7 <= 1461 &&
            ndvi <= 0.412224)
            m5_cloud_code = 50;  /* class cloud_free  [0.780] */
        
        /* Rule 4/31: (3319.6/846.9, lift 1.5) */
        if (m5_cloud_code == -1 && b1 <= 2392 && b2 > 1698 && b5 > 1516 &&
            b7 <= 1623 && ndvi <= 0.412224)
            m5_cloud_code = 50;  /* class cloud_free  [0.745] */
        
        /* Rule 4/32: (477.2, lift 2.0) */
        if (m5_cloud_code == -1 && b1 > 2005 && b2 <= 1698 && b5 > 1516)
            m5_cloud_code = 100;  /* class cloud  [0.998] */
        
        /* Rule 4/33: (261.5, lift 2.0) */
        if (m5_cloud_code == -1 && b1 > 2247 && b2 <= 1906 && b7 > 1082)
            m5_cloud_code = 100;  /* class cloud  [0.996] */
        
        /* Rule 4/34: (546.4/1.6, lift 2.0) */
        if (m5_cloud_code == -1 && b1 > 1645 && b2 <= 1292 && b5 > 796)
            m5_cloud_code = 100;  /* class cloud  [0.995] */
        
        /* Rule 4/35: (286.3/3.2, lift 2.0) */
        if (m5_cloud_code == -1 && b1 > 2040 && b2 <= 1763 && b7 > 1082 &&
            ndvi <= 0.107894)
            m5_cloud_code = 100;  /* class cloud  [0.986] */
        
        /* Rule 4/36: (262.7/10.2, lift 1.9) */
        if (m5_cloud_code == -1 && b5 <= 1516 && b7 > 1084 &&
            ndvi > 0.107894 && ndsi > -0.0853392 && ndsi <= 0.267478)
            m5_cloud_code = 100;  /* class cloud  [0.958] */
        
        /* Rule 4/37: (375.3/16.8, lift 1.9) */
        if (m5_cloud_code == -1 && b1 > 2392 && b3 <= 2260 && b5 > 1516 &&
            ndvi > 0.107894)
            m5_cloud_code = 100;  /* class cloud  [0.953] */
        
        /* Rule 4/38: (7438.5/586.8, lift 1.8) */
        if (m5_cloud_code == -1 && b2 > 3959 && b5 > 796 &&
            ndvi <= 0.0799109 && ndsi <= 0.85271)
            m5_cloud_code = 100;  /* class cloud  [0.921] */
        
        /* Rule 4/39: (863.1/107.7, lift 1.8) */
        if (m5_cloud_code == -1 && b1 > 2247 && b3 <= 2569 &&
            ndsi > -0.177533 && ndsi <= -0.115822)
            m5_cloud_code = 100;  /* class cloud  [0.874] */
        
        /* Rule 4/40: (198.3/30.9, lift 1.7) */
        if (m5_cloud_code == -1 && b1 > 1511 && b2 > 1317 && b2 <= 1387 &&
            b7 > 1461 && ndvi > 0.175138)
            m5_cloud_code = 100;  /* class cloud  [0.841] */
        
        /* Rule 4/41: (2902.8/463.9, lift 1.7) */
        if (m5_cloud_code == -1 && b1 > 2601 && b2 <= 2885 && b3 <= 2909 &&
            ndsi <= 0.143982)
            m5_cloud_code = 100;  /* class cloud  [0.840] */
        
        /* Rule 4/42: (1471.7/245.2, lift 1.7) */
        if (m5_cloud_code == -1 && b1 > 1802 && b2 <= 1584 && b5 > 796 &&
            ndsi > -0.177533 && ndsi <= 0.267478)
            m5_cloud_code = 100;  /* class cloud  [0.833] */
        
        /* Rule 4/43: (369.8/61.3, lift 1.7) */
        if (m5_cloud_code == -1 && b1 > 2040 && b1 <= 2601 && b3 <= 2569 &&
            b5 <= 1493 && b7 > 1082 && ndsi <= 0.267478)
            m5_cloud_code = 100;  /* class cloud  [0.832] */
        
        /* Rule 4/44: (3386.7/901.7, lift 1.5) */
        if (m5_cloud_code == -1 && b1 > 1511 && b7 > 1461 &&
            ndvi > 0.175138 && ndsi > -0.22841)
            m5_cloud_code = 100;  /* class cloud  [0.734] */
        
        /* Rule 4/45: (44465.5/20067, lift 1.1) */
        if (m5_cloud_code == -1 && b5 > 796 && ndsi <= 0.85271)
            m5_cloud_code = 100;  /* class cloud  [0.549] */ 
        
        if (m5_cloud_code == -1)
            m5_cloud_code = 100;  /* Default class: cloud_free */

        /* Use the maximum cloud code as the cloud value from this set of
           model runs */
        limited_cloud_code = m1_cloud_code;
        if (m2_cloud_code > limited_cloud_code)
            limited_cloud_code = m2_cloud_code;
        if (m3_cloud_code > limited_cloud_code)
            limited_cloud_code = m3_cloud_code;
        if (m4_cloud_code > limited_cloud_code)
            limited_cloud_code = m4_cloud_code;
        if (m5_cloud_code > limited_cloud_code)
            limited_cloud_code = m5_cloud_code;

        /* Use the cloud codes to determine the revised cloud codes */
        if (limited_cloud_code == 50)
            /* Original fmask identified cloud, but upon further inspection it
               looks like it's not cloudy */
            rev_cloud_mask[samp] = 0;
        else
            /* Original fmask identified cloud, and it appears to have been
               correct */
            rev_cloud_mask[samp] = OUT_CLOUD;
    }  /* for samp */
}

