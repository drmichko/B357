#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "basistools.h"

typedef unsigned char uchar;
int getinfoparms( int num, size_t *orbsize, size_t *orb, size_t *last );
uchar *loadZiptarget(char *dirname, int target, basis_t *b );
#endif
