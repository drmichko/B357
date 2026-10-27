#define _POSIX_C_SOURCE 200809L     /* avant tout #include pour avoir effet (getline) */

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>

#include "orbitData.h"
#include "boolean.h"
#include "orbitools.h"

#define M      7
#define N      (1 << M)

    int dimen = 7;
    int target = 0;
    int goal = 38;
    int nbproc = 0;
    int optwrite = 0;
    char *dirname = "./data";
    basis_t base;
    uint64_t total = 0, card;

    int verb = 0;

uchar *loadZipname(char *dirname, int target, basis_t *b)
{
    FILE *srcFile;
    size_t lus;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;
    size = size / 8;
    uchar *res;
    sprintf(fn, "%s/%d-%d-%d-%d.zip", dirname, target, b->s, b->t, b->m);
    srcFile = fopen(fn, "r");
    if ( ! srcFile ) {
	perror( fn );
	exit( 1 );
    } 
    res = calloc(size, sizeof(uchar));
    if (!res) {
        perror("allocation table");
        exit(1);
    }

    if (!srcFile) {
        perror(fn);
        exit(1);
    }
    lus = fread(res, 1, size, srcFile);
    if (lus != size) {
        printf("load problem:-(\n");
        printf("lus=%ld spaceSize=%ld\n", (long) lus, (long) size);
        exit(1);
    }
    fclose(srcFile);
    return res;
}



static size_t doit( int num, uint64_t order, boole h, aglGroup g)
{
    uchar *zip = loadZipname(dirname, num , &base);
    assert( zip != NULL );
    vector last  = 1ULL << 32;
    while ( ! zip[last] ) last--;
const uint64_t *w = (const uint64_t *)zip;
uint64_t count = 0;
#pragma omp parallel for reduction(+:count) schedule(static)
for (size_t i = 0;  i < (1ULL << 29); i++)
    count += __builtin_popcountll(w[i]);
free( zip );       
printf("#num=%3d stabSize=%ld orbits=%ld last=%ld", num, order, count, last );
panf( stdout, h );
paglGroup( stdout, g );
printf("\nstabSize=%ld\n\n", order);
return count;
}

int main(int argc, char *argv[])
{
    int opt;
    nbproc =  omp_get_max_threads();
    card   =  aglCardinality( 7 );
    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:g:ip:")) != -1) {
        switch (opt) {
        case 'w': optwrite++;break;
        case 'd': dirname = strdup( optarg ); break;
        case 'c': break;
        case 'p': nbproc = atoi( optarg );break;
        case 'm': dimen = atoi(optarg); break;
        case 'f': break;
        case 'g': goal = atoi(optarg); break;
        case 't': target = atoi(optarg); break;
        default:  exit(0);
        }
    }

    printf("\nUsing %d threads...\n", nbproc );

    char *fn = "../boole/data/B-4-5-7.dat";
    FILE *src = fopen( fn ,"r");
    if (!src) { perror(  fn ); exit(1); }

    initboole(dimen);
    initagldim(dimen);

    base = monomialBasis(3, 3, 7);
    int item = 0;
    boole h;
    size_t grpSize;
    aglGroup grp = NULL;
    while ( ( h = loadaglboolesize(src, &grp, &grpSize) )){
		item++;
		if (  item == 179  ) 
			total += doit( item, grpSize, h , grp );
		aglfreeGroup( grp );
		free(h);
    }
    fclose(src);
    printf("\n\nclass number=%ld\n", total );
    fflush(stdout);
    return 0;
}
