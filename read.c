




#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <inttypes.h>
#include "boolean.h"
#include "orbitools.h"
#include "vectorial.h"
#include "common.h"
#define BUF 65536

#define BLK 65536

// lit tout le fichier .lst ; appelle traite( v ) pour chaque mot
void read5( char *fn )
{
    FILE *src;
    unsigned char *buf;
    size_t i, n, k;
    uint64_t v, count = 0;

    src = fopen( fn, "rb" );
    if ( ! src ) { perror( fn ); exit( 1 ); }
    buf = malloc( 5 * BLK );
    if ( ! buf ) { perror( "malloc" ); exit( 1 ); }

    while ( ( n = fread( buf, 5, BLK, src ) ) > 0 ) {
        for ( i = 0; i < n; i++ ) {
            v = 0;
            for ( k = 0; k < 5; k++ )
                v |= (uint64_t) buf[5*i + k] << ( 8*k );
            //traite( v );
            count++;
        }
    }
    if ( ferror( src ) ) { perror( "fread" ); exit( 1 ); }
    free( buf );
    fclose( src );
    printf( "\n#lus:%" PRIu64 "\n", count );
}



int main(int argc, char *argv[])
{

    initboole( 7 );
    int item = atoi( argv[1] );

    size_t  orb, stab, last;
    if ( ! getinfoparms( item , &stab, &orb, &last ) ) 
	    exit(1);
    char fn[ 128 ];
    sprintf(fn, "list/%d-3-3-7.lst5",  item);
    read5( fn );
    return 0;
}
