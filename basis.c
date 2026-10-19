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

int main( void )
{
    int dimen = 7;
    initboole(dimen);
    initagldim(dimen);

    puts("standard basis of B(3,3,7):" );
    basis_t base = monomialBasis(3, 3, 7);
    for( int i = 0; i < base.mat.nbl; i++ ){
	panf( stdout, base.mat.fct[i] );
    }
    return 0;
}
