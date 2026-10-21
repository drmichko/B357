
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


void savelistname(int num, char* dirname, char *name, basis_t * b  )
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
    uint64_t *table = calloc( orb, sizeof(vector) );
    for (vector v = 0; v < size; v ++) {
            if ( b->table[v]  ) {
                table[ count++] = v;
	    }
	    if ( count == orb ) break;
    }
   
size_t w = fwrite(table, sizeof table[0], count , dstFile);
if (w != count ) { perror("fwrite"); exit(1); }
    fclose(dstFile);
    printf("\n#count:%ld\n", count );
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
    if ( orb * 8 >=  limite ) {
	    printf("large : num=%d orb=%ld size=%ld\n", item, orb, orb*8  );
	    exit(0);
    }
    char fn[ 128 ];
    basis_t base = monomialBasis( 3, 3, ffdimen);
    char *dirname="data";
    sprintf(fn, "%s/%s-%d-%d-%d.lst", dirname, cname, base.s, base.t, base.m);
    if ( access( fn , F_OK) == 0) {
            puts("already done...");
            exit(0);
    }

    loadBasename(cname, &base);
    savelistname( item, "data", cname, &base );

    return 0;
}
