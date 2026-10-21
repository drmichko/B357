
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "boolean.h"
#include "orbitools.h"
#include "vectorial.h"
#include "common.h"


void checklistname(int num, char* dirname, char *name, basis_t * b  )
{
    FILE *dstFile;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;

    size_t stab, orb, last;
    if ( ! getinfoparms( num, &stab, &orb, &last ) ) exit(1);
    sprintf(fn, "%s/%s-%d-%d-%d.lst", dirname, name, b->s, b->t, b->m);
    
    dstFile = fopen(fn, "wb");
    if (!dstFile) {
        perror(fn);
        exit(1);
    }
    uint64_t count = 0;
    uint64_t old = 0; 
    uint64_t delta = 0; 

    for (vector v = 0; v < size; v ++) {
            if  ( b->table[v]  ) {
	    	if  ( v - old > delta ) delta = v - old;
		count++;
		old = v;
	    }
	    if ( count == orb ) break;
    }
   
    printf("num=%d delta= %ld\n", num, delta );

}




int main(int argc, char *argv[])
{

    initboole( 7 );
    char cname[64];
    int item = atoi( argv[1] );
    sprintf(cname, "%d",  item  );

    size_t limite = 1ULL << 32, orb, stab, last;
    if ( ! getinfoparms( item , &stab, &orb, &last ) ) 
	    exit(1);
    basis_t base = monomialBasis( 3, 3, ffdimen);
    char *dirname="data";
    loadBasename(cname, &base);
    checklistname( item, "data", cname, &base );

    return 0;
}
