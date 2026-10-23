
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
#define BUF 65536

void write5( FILE *f, uint64_t *tab, size_t n )
{
    unsigned char *buf = malloc( 5 * BUF );
    size_t i, k, m, done = 0;

    while ( done < n ) {
        m = n - done;
        if ( m > BUF ) m = BUF;
        for ( i = 0; i < m; i++ )
            for ( k = 0; k < 5; k++ )
                buf[5*i + k] = ( tab[done + i] >> ( 8*k ) ) & 0xFF;
        if ( fwrite( buf, 5, m, f ) != m ) { perror( "fwrite" ); exit( 1 ); }
        done += m;
    }
    free( buf );
}



void savelistname(int num, char* dirname, char *name, basis_t * b  )
{
    FILE *dstFile;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;

    size_t stab, orb, last;
    if ( ! getinfoparms( num, &stab, &orb, &last ) ) exit(1);
    sprintf(fn, "%s/%s-%d-%d-%d.lst5", dirname, name, b->s, b->t, b->m);
    
    dstFile = fopen(fn, "wb");
    if (!dstFile) {
        perror(fn);
        exit(1);
    }
    uint64_t count = 0;
    uchar     *table = calloc( orb,  5  );
    for (vector v = 0; v < size  && count < orb ; v ++) {
            if ( b->table[v]  ) {
		vector tmp    = v;
		int pos=0;
		table[ pos++] =  tmp &  255; tmp >>= 8;
		table[ pos++] =  tmp &  255; tmp >>= 8;
		table[ pos++] =  tmp &  255; tmp >>= 8;
		table[ pos++] =  tmp &  255; tmp >>= 8;
		table[ pos++] =  tmp &  255; tmp >>= 8;
                count++;
	    }
    }
   
size_t w = fwrite(table, 5 , count , dstFile);
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
    if ( orb * 5 >=  limite ) {
	    printf("large : num=%d orb=%ld size=%ld\n", item, orb, orb*8  );
	    exit(0);
    }
    char fn[ 128 ];
    basis_t base = monomialBasis( 3, 3, ffdimen);
    char *dirname="list";
    sprintf(fn, "%s/%s-%d-%d-%d.lst5", dirname, cname, base.s, base.t, base.m);
    if ( access( fn , F_OK) == 0) {
            puts("already done...");
            exit(0);
    }

    loadBasename(cname, &base);
    savelistname( item, "list", cname, &base );
    puts( fn );
    return 0;
}
