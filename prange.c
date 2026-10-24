/*
 * nlq7_omp.c : test NQ(f) >= T pour des fonctions booléennes à 7 variables,
 * version parallèle OpenMP.
 * Table de vérité : f(x) = (tt[x >> 6] >> (x & 63)) & 1.
 *
 * Compilation : gcc -O3 -march=native -fopenmp -o nlq7_omp nlq7_omp.c [+ vos autres .c]
 * Exécution   : OMP_NUM_THREADS=16 ./nlq7_omp -t <target> -g <goal>
 */
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
#include "common.h"

#define M      7
#define N      (1 << M)
#define NMONO  (M + M * (M - 1) / 2)   /* 28 monômes de degré 1 ou 2 */

static uint64_t mono[NMONO][2];

/* À appeler une fois, AVANT la région parallèle. */
static void build_monomials(void)
{
    int k = 0;
    for (int i = 0; i < M; i++)
        for (int j = i; j < M; j++, k++) {      /* j == i : monôme x_i */
            mono[k][0] = mono[k][1] = 0;
            for (int x = 0; x < N; x++)
                if (((x >> i) & (x >> j)) & 1)
                    mono[k][x >> 6] |= 1ULL << (x & 63);
        }
}

/*
 * Retourne 1 si NQ(f) >= T, 0 sinon.
 * Échec dès qu'une g de degré <= 2 vérifie d(f,g) < T ou d(f,g) > N - T.
 * Parcours de Gray des 2^28 fonctions de degré <= 2 ; pure (thread-safe).
 */
static int nq_at_least_m7(const uint64_t tt[2], int T)
{
    const int dmin = T;
    const int dmax = N - T;

    uint64_t c0 = tt[0], c1 = tt[1];           /* f + g courant */

    for (uint32_t i = 1; ; i++) {
        int d = __builtin_popcountll(c0) + __builtin_popcountll(c1);
        if (d < dmin || d > dmax) return 0;    /* sortie anticipée */

        if (i >> NMONO) return 1;              /* 2^28 termes parcourus */
        const uint64_t *m = mono[__builtin_ctz(i)];
        c0 ^= m[0];
        c1 ^= m[1];
    }
}


    #include <string.h>   /* memcpy */


static inline uint64_t zip_word(const uint8_t *zip, uint64_t b, uint64_t nbytes)
{
    uint64_t w = 0;                               /* octets absents => 0 */
    uint64_t off = 8 * b;
    size_t n = (nbytes - off < 8) ? (size_t)(nbytes - off) : 8;
    memcpy(&w, zip + off, n);                     /* little-endian */
    return w;
}


int main(int argc, char *argv[])
{
    nbproc =  omp_get_max_threads();
	option( argc, argv );

    printf("\nUsing no %d  (threads: %d)\n", target, nbproc );
    initboole(dimen);
    initagldim(dimen);

    if (!target) return 0;

    uint64_t stabSize, orb, last;
    aglGroup grp = NULL;

    boole h = getfullinfos( target, &stabSize, &orb, &last, &grp );
    if ( ! h ) {
            puts("infos?");
            return 1;
    }

    basis_t base = monomialBasis(3, 3, 7);


    uchar *zip = loadZiptarget("data", target, &base);
    uint64_t orbits = 0;
    uint64_t * W = ( uint64_t *) zip;
    assert( zip != NULL );
    size_t limite = (1ULL << 32)  - 1;
    while ( ! zip[ limite ]  ) limite--;
    for( size_t i = 0; i <  (1ULL << 29 ) ; i++ )
            orbits += __builtin_popcountll( W[i]);
    assert( orbits == orb );
    printf("\nlimites = %ld %ld\n", limite, last );

    fflush(stdout);


    build_monomials();

    const int nbl = (int) base.mat.nbl;
    uint64_t (*fct64)[2] = calloc((size_t) nbl, sizeof *fct64);
    if (!fct64) { perror("fct64"); exit(1); }
    for (int j = 0; j < nbl; j++)
        for (int x = 0; x < ffsize; x++)
            if (base.mat.fct[j][x])
                fct64[j][x >> 6] |= 1ULL << (x & 63);

    uint64_t h64[2] = {0, 0};
    for (int x = 0; x < ffsize; x++)
        if (h[x]) h64[x >> 6] |= 1ULL << (x & 63);

    /* --- Solutions : on stocke les v, affichage séquentiel et trié après --- */
    size_t cap = 1024, nsol = 0;
    uint64_t *sol = malloc(cap * sizeof *sol);
    if (!sol) { perror("sol"); exit(1); }

    size_t count = 0;
    double t0 = omp_get_wtime();

   

const uint64_t nbytes =  last + 1;     
const uint64_t nwords = (nbytes + 7) / 8;

#pragma omp parallel for schedule(dynamic, 64) reduction(+:count)
for (uint64_t b = 0; b < nwords; b++) {
    uint64_t w = zip_word((const uint8_t *)zip, b, nbytes);

    while (  w   ) {
        int t = __builtin_ctzll(w);
        w &= w - 1;
        vector v = (vector)(64 * b + t);

        count++;

        uint64_t tt[2] = { h64[0], h64[1] };
        for (int j = 0; j < nbl; j++)
            if (BITVECTOR(j, v)) {
                tt[0] ^= fct64[j][0];
                tt[1] ^= fct64[j][1];
            }

        if (nq_at_least_m7(tt, goal)) {
            #pragma omp critical(sol)
            {
                if (nsol == cap) {
                    cap *= 2;
                    sol = realloc(sol, cap * sizeof *sol);
                    if (!sol) { perror("realloc"); exit(1); }
                }
                sol[nsol++] = (uint64_t) v;
            }
        }
    }
}




    printf("candidat:%ld\n", nsol );
    double temps = omp_get_wtime() - t0;
    char name[128];
    sprintf(name, "results/res-%d.txt", target);
    FILE *dst = stdout;
    if ( optw ) dst= fopen( name, "a" );;
    /* --- Affichage séquentiel, dans l'ordre des v (comme la version série) --- */
    fprintf(dst, "\nUsing no %d  (threads: %d)\n", target, nbproc );
    fprintf(dst, "\n#%zu members have NL_2=38 :\n", nsol);
    for (size_t i = 0; i < nsol; i++) {
        boole f = vectortoboole((vector) sol[i], &base);
        for (int x = 0; x < ffsize; x++) f[x] ^= h[x];
        panf( dst, f);
        free(f);
    }
    fflush(dst);

    fprintf(dst, "\n#orbit=%" PRIu64 " count=%zu  %zu\n", orbits, count, nsol);
    fprintf(dst, "\n#runtime : %.6f secondes (mur)\n", temps);

    fclose( dst );
    free(sol);
    free(fct64);
    free(zip);

    return 0;
}
