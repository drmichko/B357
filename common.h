#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "basistools.h"

extern int dimen ;
extern int target ;
extern int goal ;
extern int nbproc;
extern int optw;

typedef unsigned char uchar;
int getinfoparms( int num, size_t *orbsize, size_t *orb, size_t *last );
boole getfullinfos( int num, size_t *stabSize, size_t *orb, size_t *last, aglGroup *g );
uchar *loadZiptarget(char *dirname, int target, basis_t *b );
void option( int argc, char** argv );
#endif
