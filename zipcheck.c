
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


void checkzip(int num, basis_t * b  )
{
    size_t stab, orb, last;
    if ( ! getinfoparms( num, &stab, &orb, &last ) ) exit(1);
    
    uint64_t count = 0;
    uint64_t *W = ( uint64_t *)  b->table;
    for (size_t i = 0;  i < (1ULL << 29); i++)
            count += __builtin_popcountll( W[i] );
    
   
    if ( count != orb )
	    printf("num=%d %ld %ld\n", num, count , orb );
    else
	    printf("num=%d is ok \n", num );

}




int main(int argc, char *argv[])
{

    initboole( 7 );
    int item = atoi( argv[1] );

    size_t  orb, stab, last;
    if ( ! getinfoparms( item , &stab, &orb, &last ) ) 
	    exit(1);
    basis_t base = monomialBasis( 3, 3, ffdimen);
    uchar * zip = loadZiptarget( "../tools/data", item , &base);
    base.table = zip;
    checkzip( item, &base );

    return 0;
}
