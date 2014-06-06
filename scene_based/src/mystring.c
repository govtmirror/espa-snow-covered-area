#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "mystring.h"

/******************************************************************************
MODULE:  dup_string

PURPOSE:  Duplicates a string, including allocating memory for the returned
pointer.

RETURN VALUE:
Type = char *
Value      Description
-----      -----------
NULL       Error occurred duplicating the string or the input string was NULL
non-NULL   Successful completion
 

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/3/2012    Gail Schmidt     Original Development (based on input routines
                             from the LEDAPS lndsr application)

NOTES:
  1. This routine allocates memory for the returned pointer.  It is up to
     the caller to free that memory when done using the pointer.
******************************************************************************/
char *dup_string
(
    char *string    /* I: string to be duplicated */
)
{
    char FUNC_NAME[] = "dup_string";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int len;          /* length of the string */
    char *s = NULL;   /* string pointer to use and then return to caller */
  
    /* Make sure the input string isn't NULL */
    if (string == (char *) NULL)
        return ((char *) NULL);

    /* Get the string length */
    len = strlen (string);
    if (len < 0) 
    {
        strcpy (errmsg, "Invalid string length");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    /* Allocate memory for the length plus one for the string terminator */
    s = (char *) calloc ((len + 1), sizeof(char));
    if (s == (char *)NULL) 
    {
        strcpy (errmsg, "Unable to allocate memory for the string.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Copy the string */
    if (strcpy (s, string) != s)
    {
        free (s);
        strcpy (errmsg, "Unable to copy the string.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    return (s);
}
