#ifndef _ERROR_HANDLER_H_
#define _ERROR_HANDLER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bool.h"

/* Set up default global defines */
#define SUCCESS 0
#define ERROR 1

/* Prototypes */
void error_handler
(
    bool error_flag,  /* I: true for errors, false for warnings */
    char *module,     /* I: calling module name */
    char *errmsg      /* I: error message to be printed, without ending EOL */
);

#endif
