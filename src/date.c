#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "date.h"

/******************************************************************************
MODULE:  date_init

PURPOSE:  Initializes the date structure based on the input string and
specified format

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred processing the date
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int date_init
(
    char *s,               /* I: date string to use for initialization */
    Date_format_t iformat, /* I: which date format is input? */
    Date_t *this           /* I/O: pointer to date structure to be populated */
)
{
    char FUNC_NAME[] = "date_init";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    char *date = NULL;      /* pointer to the date section of the string */
    char *time = NULL;      /* pointer to the time section of the string */
    bool leap;              /* is this a leap year */
    int year1;              /* number of years since 1900 */
    int len;                /* length of input date string */
    int nday[12] = {31, 29, 31, 30,  31,  30,  31,  31,  30,  31,  30,  31};
                   /* number of days in each month */
    int idoy[12] = { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336};
                   /* starting DOY for each month */
    int jleap, idoy_nonleap;
  
    /* Validate the format parameter */
    this->valid = false;
    if (iformat != DATE_FORMAT_DATEA_TIME && 
        iformat != DATE_FORMAT_DATEB_TIME &&
        iformat != DATE_FORMAT_DATEA &&
        iformat != DATE_FORMAT_DATEB)
    {
        strcpy (errmsg, "Invalid date format parameter");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Check for validity in the string length and the format type then
       start parsing the string based on the format type */
    len = strlen (s);
    if (iformat == DATE_FORMAT_DATEA_TIME)
    {
        if (len < 20 || len > 27) 
        {
            sprintf (errmsg, "Invalid date/time string length. Expected %s.",
                DATE_FORMAT_DATEA_TIME_STR);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (s[10] != 'T' || s[len-1] != 'Z')
        {
            sprintf (errmsg, "Invalid date/time format. Expected %s.",
                DATE_FORMAT_DATEA_TIME_STR);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        date = &s[0];
        time = &s[11];
    }
    else if (iformat == DATE_FORMAT_DATEB_TIME)
    {
        if (len < 18 || len > 25) 
        {
            sprintf (errmsg, "Invalid date/time string length. Expected %s.",
                DATE_FORMAT_DATEB_TIME_STR);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (s[8] != 'T' || s[len-1] != 'Z')
        {
            sprintf (errmsg, "Invalid date/time format. Expected %s.",
                DATE_FORMAT_DATEB_TIME_STR);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        date = &s[0];
        time = &s[9];
    }
    else if (iformat == DATE_FORMAT_DATEA)
    {
        if (len != 10) 
        {
            sprintf (errmsg, "Invalid date/time string length. Expected %s.",
                DATE_FORMAT_DATEA_STR);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        date = s;
    }
    else if (iformat == DATE_FORMAT_DATEB)
    {
        if (len != 8) 
        {
            sprintf (errmsg, "Invalid date/time string length. Expected %s.",
                DATE_FORMAT_DATEB_STR);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        date = s;
    }
  
    if (iformat == DATE_FORMAT_DATEA_TIME  ||
        iformat == DATE_FORMAT_DATEA)
    {
        if (sscanf (date, "%4d-%2d-%2d", &this->year, &this->month, &this->day)
            != 3) 
        {
            sprintf (errmsg, "Invalid date format.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (this->year < 1900 || this->year > 2400) 
        {
            sprintf (errmsg, "Invalid year: %d.", this->year);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (this->month < 1 || this->month > 12) 
        {
            sprintf (errmsg, "Invalid month: %d.", this->month);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (this->day < 1 || this->day > nday[this->month-1])
        {
            sprintf (errmsg, "Invalid day of month: %d.", this->day);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        this->doy = this->day + idoy[this->month - 1] - 1;
    }
    else
    {
        if (sscanf(date, "%4d-%3d", &this->year, &this->doy) != 2) 
        {
            sprintf (errmsg, "Invalid date format.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (this->year < 1900 || this->year > 2400) 
        {
            sprintf (errmsg, "Invalid year: %d.", this->year);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (this->doy < 1 || this->doy > 366) 
        {
            sprintf (errmsg, "Invalid day of year: %d.", this->doy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }
  
    /* Check if the year is a leap year */
    if (this->year % 4 == 0 && (this->year % 100 != 0 || this->year % 400 == 0))
        leap = true;
    else
        leap = false;

    if (iformat == DATE_FORMAT_DATEA_TIME  ||
        iformat == DATE_FORMAT_DATEA)
    {
        if ((this->month == 2) && !leap && (this->day > 28))
        {
            sprintf (errmsg, "Invalid day of month in February for a non-leap "
                "year: %d.", this->month);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (!leap && (this->month > 2))
            this->doy--;
    }
    else
    {
        if (leap)
        {
            for (this->month = 0; this->month < 12; this->month++)
                if (this->doy < idoy[this->month])
                    break;
        }
        else
        {
            if (this->doy > 365) 
            {
                sprintf (errmsg, "Invalid day of year: %d.", this->doy);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            for (this->month = 0; this->month < 12; this->month++)
            {
                idoy_nonleap = idoy[this->month];
                if (this->month > 1)
                    idoy_nonleap--;
                if (this->doy < idoy_nonleap)
                    break;
            }
        }
    }
  
    /* Convert to Julian days ca. 2000 (1 = Jan. 1, 2000) */
    year1 = this->year - 1900;
    if (year1 > 0)
    {
        jleap = (year1 - 1) / 4;
        if (this->year > 2100)
            jleap -= (this->year - 2001) / 100;
    }
    else
    {
        jleap = 0;
    }
    this->jday2000 = (year1 * 365) + jleap + this->doy;
    this->jday2000 -= 36524;
  
    /* Parse and check time */
    if (time != (char *)NULL)
    {
        if (sscanf (time, "%2d:%2d:%lf", &this->hour, &this->minute,
            &this->second) != 3)
        {
            sprintf (errmsg, "Invalid time format.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }
    else
    {
        this->hour = this->minute = 0;
        this->second = 0.0;
    }
    if (this->hour < 0 || this->hour > 23) 
    {
        sprintf (errmsg, "Invalid hour: %d", this->hour);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (this->minute < 0 || this->minute > 59)
    {
        sprintf (errmsg, "Invalid minute: %d", this->minute);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (this->second < 0.0 || this->second > 59.999999)
    {
        sprintf (errmsg, "Invalid second: %f", this->second);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Convert to seconds of day */
    this->sod = (((this->hour * 60.0) + this->minute) * 60.0) + this->second;
  
    /* Date structure is now valid */
    this->valid = true;
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  format_date

PURPOSE:  Formats an output string for the date in the date structure based
on the desired date format

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred processing the date
SUCCESS    Successful completion

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
******************************************************************************/
int format_date
(
    Date_t *this,           /* I: pointer to the date structure */
    Date_format_t iformat,  /* I: which date format should be output? */
    char *s                 /* O: correctly-formatted output date string */
)
{
    char FUNC_NAME[] = "format_date";   /* function name */
    char errmsg[STR_SIZE];        /* error message */

    if (this == (Date_t *) NULL || !this->valid)
    {
        sprintf (errmsg, "Invalid date structure or structure not initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    
    if (iformat == DATE_FORMAT_DATEA_TIME)
        sprintf (s, "%4d-%02d-%02dT%02d:%02d:%09.6fZ", this->year,
            this->month, this->day, this->hour, this->minute, this->second);
    else if (iformat == DATE_FORMAT_DATEB_TIME)
        sprintf (s, "%4d-%03dT%02d:%02d:%09.6fZ", this->year, this->doy, 
      	    this->hour, this->minute, this->second);
    else if (iformat == DATE_FORMAT_DATEA)
        sprintf (s, "%4d-%02d-%02d", this->year, this->month, this->day);
    else if (iformat == DATE_FORMAT_DATEB)
        sprintf (s, "%4d-%03d", this->year, this->doy);
    else if (iformat == DATE_FORMAT_TIME)
        sprintf (s, "%02d:%02d:%09.6f", this->hour, this->minute, this->second);
    else 
    {
        sprintf (errmsg, "Invalid format parameter.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    
    return (SUCCESS);
}
