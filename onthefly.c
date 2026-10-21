
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


int main(int argc, char *argv[])
{  

    char *fn = "../boole/data/B-4-5-7.dat";


    aglGroup g = NULL;
    int opt, optw = 0; 
    int dimen = 7;
    int target = 0;
    int nbthread = omp_get_max_threads();;
    while ((opt = getopt(argc, argv, "i:f:p:wb:d:m:r:t:")) != -1) {
	switch (opt) {
	case 'd':
	   break;
	case 'p':
	    nbthread = atoi( optarg );
	    break;
	case 'w':
	   	optw++; 
	break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    break;
	case 't':
	    target=atoi( optarg );
	    break;
	default:
	    exit(0);
	}

    }


    printf("\nthread=%d\n", nbthread );
    FILE *src = fopen( fn, "r" );
    initboole(  dimen );
    initagldim( dimen );

    boole f;
    int item = 0;
    uint64_t size = 0;
    basis_t base = monomialBasis( 3, 3, 7 );
    initBrowse( &base );
    while ((f = loadaglboolesize( src , &g, &size ))) {
		item++;
	        //panf( stderr, f );	
		//fprintf(stderr,  "\nitem=%d degree=%d order=%ld\n", item, degree(f), size );
		//paglGroup( stderr, g );
                if ( item == target ) 
			break;
		aglfreeGroup( g );
		free(f);
		g=NULL;
    }
    if ( item != target ){
	    exit(0);
    }
    if ( size > 1024 ) {
	    puts("too large");
	    exit ( 0 );
    }
    fclose(src);

    FILE* dst = stdout;

    double debut = omp_get_wtime();
 
    panf( dst , f );
    fprintf( dst, "\nitem=%d degree=%d order=%ld\n", item, degree(f), size );
    paglGroup( stderr , g );
    assert( degree( f ) > 4  );
    assert(  checkstab( f, g,  3) == 1 );
    agl*  table = devGroup( g, &size );
    aglVectorGroup G = aglBoundaryGroupAction( f, g, &base );
    size_t  count = 0;

    aglVector *L = calloc( size, sizeof( aglVector ) );
    for( size_t i = 0;  i < size; i++ ) {
        assert(  checks( f, table[i],  3) == 1 );
    	L[i]   = myboundaryAction(  f, table[ i ], &base );
	vector tmp = nbfixed( L[i] , 35  );
	count+= tmp;
    }
    assert( count % size  == 0 );  
    size_t burnside = count / size;
    fprintf( stderr, "\n#orbits=%ld\n", burnside );

omp_set_num_threads(nbthread);


#pragma omp for schedule(guided) 
    for (vector v = 0; v < base.size; v++) {
        if (base.table[v]) continue;
        int is_rep = 1;
        for (size_t i = 0; i < size; i++) {
            vector y = aglVectorImages(v, L[i], 35);
            if (y < v) { is_rep = 0; break; }
            if (y > v) base.table[y] = 1;
        }
        if (is_rep) {
            base.table[v] = 2;
        }
    }

size_t reps = 0;

    for( vector v = 0; v < base.size; v++ ){
            base.table[v] =  ( base.table[v] == 2 );
            reps += base.table[v] ;
    }

    assert( reps == burnside );
    double fin = omp_get_wtime();
    double temps = (double)(fin - debut);
    fprintf(dst, "\nruntime : %.1f secondes\n", temps);


    aglfreeGroup( g );
    aglVectorGroupFree(  G );
    free( base.table );
    
    free( f );
    fclose(dst);
    return 0;
}
