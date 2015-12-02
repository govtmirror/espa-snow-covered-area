#ifndef OUTPUT_H
#define OUTPUT_H

#include "common.h"
#include "input.h"

#define MAX_DATE_LEN (28)

/* Define the number of bands that might be output to the file */
#define MAX_OUT_BANDS NUM_CM

/* Define some of the constants to use in the output data products */
#define FILL_VALUE -9999
#define CFMASK_FILL_VALUE 255
#define SATURATE_VALUE 20000
#define FLOAT_TO_INT 10000.0
#define SCALE_FACTOR 0.0001

/* Structure for the 'output' data type */
typedef struct {
  bool open;            /* Flag to indicate whether output file is open;
                           'true' = open, 'false' = not open */
  int nband;            /* Number of output image bands */
  int nlines;           /* Number of output lines */
  int nsamps;           /* Number of output samples */
  Espa_internal_meta_t metadata;  /* Metadata container to hold the band
                           metadata for the output bands; global metadata
                           won't be valid */
  FILE *fp_bin[MAX_OUT_BANDS];  /* File pointer for binary files */
} Output_t;

/* Prototypes */
Output_t *open_output
(
    Espa_internal_meta_t *in_meta,  /* I: input metadata structure */
    Input_t *input,                 /* I: input reflectance band data */
    int nband,                      /* I: number of bands to be created */
    char short_names[][STR_SIZE],   /* I: array of short names for new bands */
    char long_names[][STR_SIZE],    /* I: array of long names for new bands */
    char data_units[][STR_SIZE],    /* I: array of data units for new bands */
    bool toa         /* I: are we processing TOA reflectance data, otherwise
                           process surface reflectance data */
);

int close_output
(
    Output_t *this    /* I/O: Output data structure to close */
);

int free_output
(
    Output_t *this    /* I/O: Output data structure to free */
);

int put_output_lines
(
    Output_t *this,    /* I: Output data structure; buf contains the line to
                             be written */
    void *buf,         /* I: buffer to be written */
    int iband,         /* I: current band to be written (0-based) */
    int iline,         /* I: current line to be written (0-based) */
    int nlines,        /* I: number of lines to be written */
    int nbytes         /* I: number of bytes per pixel in this band */
);

int get_output_lines
(
    Output_t *this,  /* I: pointer to output data structure */
    int iband,       /* I: current band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int nbytes,      /* I: number of bytes per pixel in this band */
    void *buf        /* I: pointer to the buffer to be returned */
);

char *upper_case_str
(
    char *str    /* I: string to be converted to upper case */
);

#endif
