#ifndef _MYSTRING_H
#define _MYSTRING_H

#include <stdio.h>
#include "error_handler.h"

#define SCA_VERSION "1.1.2"

#define STR_SIZE 1024

/* Key string type definition */
typedef struct {
  int key;
  char *string;
} Key_string_t;

/* Prototypes */
char *dup_string
(
    char *string    /* I: string to be duplicated */
);

#endif

