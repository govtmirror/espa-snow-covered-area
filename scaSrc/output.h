#ifndef OUTPUT_H
#define OUTPUT_H

#include "mystring.h"
#include "input.h"
#include "space.h"

/* Define the number of SDS that will be output to the HDF-EOS file.  The
   actual SDS names are defined at the top of scene_based_sca.c. */
#define NUM_OUT_SDS 5

/* Structure for the 'output' data type */
typedef struct {
  char *file_name;      /* Output file name */
  bool open;            /* Flag to indicate whether output file is open 
                           for access; 'true' = open, 'false' = not open */
  int nband;            /* Number of output image bands */
  Img_coord_int_t size; /* Output image size */
  int32 sds_file_id;    /* SDS file id */
  Myhdf_sds_t sds[NUM_OUT_SDS]; /* SDS data structures for image data */
  uint8 *buf[NUM_OUT_SDS]; /* Output data buffer */
} Output_t;

/* Prototypes */
int create_output
(
    char *file_name    /* I: name of HDF file to be created */
);

Output_t *open_output
(
    char *file_name,                /* I: name of output HDF file */
    int nband,                      /* I: number of image bands (SDSs) to be
                                          created */
    char *sds_names[NUM_OUT_SDS],   /* I: array of SDS names for each band */
    int nlines,                     /* I: number of lines in image */
    int nsamps                      /* I: number of samples in image */
);

int close_output
(
    Output_t *this    /* I/O: Output data structure to close */
);

int free_output
(
    Output_t *this    /* I/O: Output data structure to free */
);

int put_output_line
(
    Output_t *this,    /* I: Output data structure; buf contains the line to
                             be written */
    int iband,         /* I: current band to be written (0-based) */
    int iline,         /* I: current line to be written (0-based) */
    int nlines         /* I: number of lines to be written */
);

int put_metadata
(
    Output_t *this,      /* I: Output data structure */
    int nband,           /* I: number of bands to write */
    char *band_names[NUM_OUT_SDS],  /* I: band names to write */
    char *QA_on[NUM_OUT_SDS],  /* I: metadata info for the current band "on" */
    char *QA_off[NUM_OUT_SDS], /* I: metadata info for the current band "off" */
    Input_meta_t *meta         /* I: metadata to be written */
);

#endif
