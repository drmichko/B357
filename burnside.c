
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"
#include <omp.h>

#include "boolean.h"
#include "orbitools.h"
#include <omp.h>

#include "boolean.h"
#include "orbitools.h"
#include "common.h"


#include <stdatomic.h>   


int checks( boole f, agl L, int r )
{
    boole h = getboole( );
    for ( shortvec x = 0; x < ffsize; x++)
        h[x] = f[  aglImage( x,  L  ) ] ^ f[x];
    if( degree(h) >  r ) {
                printf("\ndegrees : f=%d  b=%d\n", degree(f),  degree(h) );
		panf( stdout, f );
		panf( stdout, h );
		puts("\n");
            return 0;
    }
    free( h );
return 1;
}
int checkstab( boole f, aglGroup g, int r )
{
while ( g ) {
    boole h = getboole( );
    for ( shortvec x = 0; x < ffsize; x++)
        h[x] = f[  aglImage(x, g->per ) ] ^ f[x];
    if( degree(h) >  r ) {
                printf("\ndegrees : f=%d  b=%d\n", degree(f),  degree(h) );
		panf( stdout, f );
		panf( stdout, h );
		puts("\n");
            return 0;
    }
    free( h );
    g = g -> next;
}

return 1;
}

int aglrankvector (vector *s, int dim )
{   
int i, l;
vector  bit, x[ dim ]; 
for( i = 0;  i < dim; i++ )
	x[i] = s[i];
i=0;
while ( i < dim ) {
  if ( x[i] ) {
        bit = x[i];
        bit = bit ^ ( bit & ( bit - 1) );
        for( l = 0; l < dim; l++)
          if ( l != i )
                if ( bit & x[l] ) x[l] ^= x[i];
        i++;
   } else {
        dim = dim - 1;
        x[i] = x[ dim ];
   }
 }
return dim;
}


uint64_t nbfixed(const uint64_t *s, int dim)   // dim = 35, s a 36 éléments
{
    const uint64_t mask = (1ULL << dim) - 1;
    uint64_t basis[64] = {0};                  // basis[h] : vecteur de bit de poids fort h
    int rank = 0;

    for (int i = 0; i < dim; i++) {
        uint64_t v = (s[i] & mask) ^ (1ULL << i);   // colonne i de L+I
        while (v) {
            int h = 63 - __builtin_clzll(v);
            if (!basis[h]) { basis[h] = v; rank++; break; }
            v ^= basis[h];
        }
    }

    uint64_t b = s[dim] & mask;                // f(0)
    while (b) {
        int h = 63 - __builtin_clzll(b);
        if (!basis[h]) return 0;               // b hors de l'image : aucun point fixe
        b ^= basis[h];
    }
    return 1ULL << (dim - rank);
}

aglVector myboundaryAction(boole f, agl s, basis_t * b)
{
    // assuming  f o s = f + lower degree terms...
    int j;
    shortvec x;
    aglVector res = getaglVector( b->mat.nbl  + 1 );
    boole h = calloc( b->q, sizeof(uchar) );

    for (x = 0; x < b->q; x++)
        h[x] = f[  aglImage(x, s) ] ^ f[x];


    assert( degree(h) <=   b->t );
	projboole(b->s, b->t, h );
    res[ b->mat.nbl ] = booleVector( h , b );

    for (j = 0; j < b->mat.nbl; j++) {
        for (x = 0; x < b->q ; x++)
            h[x] = b->mat.fct[j][aglImage(x, s)];
	projboole( b->s, b->t, h );
        res[j] = booleVector(h, b );
        }

    free(h);
    return res;
}



    printf("\nthread=%d\n", nbthread );

    initboole(  dimen );
    initagldim( dimen );

    uint64_t stabSize, orb, last;
    aglGroup grp = NULL;

    boole f = getfullinfos( target, &stabSize, &orb, &last, &grp );
    if ( ! f ) {
	    puts("infos?");
	    return 1;
    }
   
    basis_t base = monomialBasis( 3, 3, 7 );
    initBrowse( &base );


    if ( stabSize > 1024 ) {
	    puts("too large");
	    exit ( 0 );
    }

    exit(0);
    FILE* dst = stdout;

    double debut = omp_get_wtime();
 
    panf( dst , f );
    fprintf( dst, "\nitem=%d degree=%d order=%ld\n", target , degree(f), stabSize );
    //paglGroup( stderr , g );
    //assert( degree( f ) > 4  );
    //assert(  checkstab( f, g,  3) == 1 );
    agl*  table = devGroup( grp, &stabSize );
    aglVectorGroup G = aglBoundaryGroupAction( f, grp, &base );
    size_t  count = 0;

    aglVector *L = calloc( stabSize, sizeof( aglVector ) );
    for( size_t i = 0;  i < stabSize; i++ ) {
        assert(  checks( f, table[i],  3) == 1 );
    	L[i]   = myboundaryAction(  f, table[ i ], &base );
	vector tmp = nbfixed( L[i] , 35  );
	count+= tmp;
    }
    assert( count % stabSize  == 0 );  
    size_t burnside = count / stabSize;
    fprintf( stderr, "\n#orbits=%ld\n", burnside );

omp_set_num_threads(nbthread);

}
