#ifndef _DATE_H
#define _DATE_H

#include "bool.h"
#include "mystring.h"
#include "error_handler.h"

/* Date/time type definition */
#define MAX_DATE_LEN (28)

#define DATE_FORMAT_DATEA_TIME_STR  "yyyy-mm-ddThh:mm:ss.ssssssZ"
#define DATE_FORMAT_DATEB_TIME_STR  "yyyy-dddThh:mm:ss.ssssssZ"
#define DATE_FORMAT_DATEA_STR       "yyyy-mm-dd"
#define DATE_FORMAT_DATEB_STR       "yyyy-ddd"
#define DATE_FORMAT_TIME_STR        "hh:mm:ss.ssssss"

typedef enum {
    DATE_FORMAT_DATEA_TIME,  /* yyyy-mm-ddThh:mm:ss.ssssssZ */
    DATE_FORMAT_DATEB_TIME,  /* yyyy-dddThh:mm:ss.ssssssZ */
    DATE_FORMAT_DATEA,       /* yyyy-mm-dd */
    DATE_FORMAT_DATEB,       /* yyyy-ddd */
    DATE_FORMAT_TIME         /* hh:mm:ss.ssssss */
} Date_format_t;

typedef struct {
    bool valid;              /* is this structure valid / populated? */
    int year;                /* year */
    int doy;                 /* day of year */
    int month;               /* month */
    int day;                 /* day of month */
    int hour;                /* hour */
    int minute;              /* minute */
    double second;           /* second */
    long jday2000;           /* julian day circa Jan. 1, 2000 */
    double sod;              /* seconds of day */
} Date_t;

/* Prototypes */
int date_init
(
    char *s,               /* I: date string to use for initialization */
    Date_format_t iformat, /* I: which date format is input? */
    Date_t *this           /* I/O: pointer to date structure to be populated */
);

int format_date
(
    Date_t *this,           /* I: pointer to the date structure */
    Date_format_t iformat,  /* I: which date format should be output? */
    char *s                 /* O: correctly-formatted output date string */
);

#endif
