#include "sca.h"

/******************************************************************************
MODULE:  write_envi_hdr

PURPOSE:  Writes the ENVI header to the specified file using the input info
provided.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred generating the header file
SUCCESS         Header file was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/8/2013    Gail Schmidt     Original Development

NOTES:
  1. It's assumed the header file will be for unsigned byte products and
     therefore an ENVI data type of 1.
******************************************************************************/
int write_envi_hdr
(
    char *hdr_file,     /* I: name of header file to be generated */
    Input_t *toa_input, /* I: input structure for both the TOA reflectance
                              and brightness temperature products */
    Space_def_t *space_def /* I: spatial definition information */
)
{
    char FUNC_NAME[] = "write_envi_hdr";   /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char *bin_file=NULL;     /* name of raw binary file */
    FILE *hdr_fptr=NULL;     /* file pointer to the ENVI header file */

    /* Open the header file */
    hdr_fptr = fopen (hdr_file, "w");
    if (hdr_fptr == NULL)
    {
        sprintf (errmsg, "Error opening %s for write access.", hdr_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Verify the projection is UTM and sphere is WGS-84 */
    if (space_def->proj_num != 1)
    {
        sprintf (errmsg, "Error UTM projection code (1) expected.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (space_def->sphere != 12)
    {
        sprintf (errmsg, "Error WGS-84 sphere code (12) expected.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }


    /* Create the name of the associated raw binary file by replacing the
       .hdr of the header filename with .bin for the output binary file */
    bin_file = strdup (hdr_file);
    strncpy (&bin_file[strlen(hdr_file)-3], "bin", 3);

    /* Write the header to the file */
    fprintf (hdr_fptr,
        "ENVI\n"
        "description = {%s}\n"
        "samples = %d\n"
        "lines   = %d\n"
        "bands   = 1\n"
        "header offset = 0\n"
        "file type = ENVI Standard\n"
        "data type = 1\n"
        "interleave = bsq\n"
        "byte order = 0\n"
        "map info = {UTM, 1, 1, %f, %f, %f, %f, %d, North, WGS-84}\n"
        "band names = {Band 1}\n", bin_file, toa_input->nsamps,
            toa_input->nlines, space_def->ul_corner.x, space_def->ul_corner.y,
            space_def->pixel_size, space_def->pixel_size, space_def->zone);

    /* Close the header file and free pointers */
    fclose (hdr_fptr);
    if (bin_file)
        free (bin_file);

    /* Successful completion */
    return (SUCCESS);
}
