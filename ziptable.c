
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "boolean.h"
#include "orbitools.h"
#include "vectorial.h"


void saveZipname(char* dirname, char *name, basis_t * b  )
{
    FILE *dstFile;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;
    sprintf(fn, "%s/%s-%d-%d-%d.zip", dirname, name, b->s, b->t, b->m);
    dstFile = fopen(fn, "w");
    if (!dstFile) {
        perror(fn);
        exit(1);
    }
    uint64_t count = 0;
    uchar *z = calloc(size / 8, 1);
    for (vector v = 0; v < size; v += 8) {
        uchar val = 0;
        for (unsigned int i = 0; i < 8; i++)
            if ( b->table[v + i]  ) {
                val |= ((uchar) 1) << i;
		count++;
	    }
        z[v / 8] = val;
    }
    size = size / 8;
    if (fwrite(z, 1, size, dstFile) != size)
        printf("problem\n");
    fclose(dstFile);
    free(z);
    printf("\n#count:%ld", count );
}




int main(int argc, char *argv[])
{

    initboole( 7 );

    basis_t base = monomialBasis( 3, 3, ffdimen);
    char cname[64];
    sprintf(cname, "%d",  atoi( argv[1] ) );
    loadBasename(cname, &base);
    saveZipname( "data", cname, &base );

    return 0;
}
