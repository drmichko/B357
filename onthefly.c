
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

    FILE* dst = stdout;

    double debut = omp_get_wtime();
 
    panf( dst , f );
    fprintf( dst, "\nitem=%d degree=%d order=%ld\n", target , degree(f), stabSize );
    
    size_t tmp;
    agl* table = devGroup(  grp, &tmp );
    assert ( tmp == stabSize ); 
    aglVector *L = calloc( stabSize, sizeof( aglVector ) );
    for (size_t i = 0; i < stabSize; i++) 
	    L[i] = boundaryAction( f , table[i], &base );

    // here we checked using Burnside formula...
omp_set_num_threads(nbthread);

#pragma omp parallel for schedule(guided)
    for (vector v = 0; v < base.size; v++) {
        if (base.table[v]) continue;
        int is_rep = 1;
        for (size_t i = 0; i < stabSize; i++) {
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

    assert( reps == orb);
    double fin = omp_get_wtime();
    double temps = (double)(fin - debut);
    fprintf(dst, "\nruntime : %.1f secondes\n", temps);


    aglfreeGroup( grp );
    free( base.table );
    
    free( f );
    fclose(dst);
    return 0;
}
