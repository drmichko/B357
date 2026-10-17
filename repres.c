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


uchar *loadZipname(char *dirname, int target, basis_t *b)
{
    FILE *srcFile;
    size_t lus;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;
    size = size / 8;
    uchar *res;
    sprintf(fn, "%s/%d-%d-%d-%d.zip", dirname, target, b->s, b->t, b->m);
    srcFile = fopen(fn, "r");
    if ( ! srcFile ) {
	perror( fn );
	exit( 1 );
    } 
    res = calloc(size, sizeof(uchar));
    if (!res) {
        perror("allocation table");
        exit(1);
    }

    if (!srcFile) {
        perror(fn);
        exit(1);
    }
    lus = fread(res, 1, size, srcFile);
    if (lus != size) {
        printf("load problem:-(\n");
        printf("lus=%ld spaceSize=%ld\n", (long) lus, (long) size);
        exit(1);
    }
    fclose(srcFile);
    return res;
}

int value(FILE *src, uint64_t *v)
{
    char  *line = NULL;
    size_t cap  = 0;
    while (getline(&line, &cap, src) != -1)
        for (const char *p = line; *p; p++)
            if (sscanf(p, "#orbits=%" SCNu64, v) == 1) {
                free(line);
                return 1;
            }

    free(line);
    return 0;
}

static vector find_end(const uchar *zip, vector size, uint64_t orbits)
{
    if (orbits == 0) return 0;
    uint64_t count = 0;
    vector nbytes = size / 8;
    for (vector b = 0; b < nbytes; b++) {
        int pc = __builtin_popcount(zip[b]);
        if (count + pc >= orbits) {
            for (int k = 0; k < 8; k++)
                if ((zip[b] >> k) & 1)
                    if (++count == orbits) return b * 8 + k + 1;
        }
        count += pc;
    }
    return size;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *) a, y = *(const uint64_t *) b;
    return (x > y) - (x < y);
}

int main(int argc, char *argv[])
{
    int opt;
    int dimen = 7;
    int target = 0;
    int goal = 38;
    int nbproc =  omp_get_max_threads();
    int optwrite = 0;
    char *dirname = "data";
    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:g:ip:")) != -1) {
        switch (opt) {
        case 'w': optwrite++;break;
        case 'd': dirname = strdup( optarg ); break;
        case 'c': break;
        case 'p': nbproc = atoi( optarg );break;
        case 'm': dimen = atoi(optarg); break;
        case 'f': break;
        case 'g': goal = atoi(optarg); break;
        case 't': target = atoi(optarg); break;
        default:  exit(0);
        }
    }

    if (!target) return 0;

    char name[128];
    sprintf(name, "rho/rho-%d.out", target);

    FILE *src = fopen(name, "r");
    if (!src) { perror(name); exit(1); }
    initboole(dimen);
    initagldim(dimen);

    basis_t base = monomialBasis(3, 3, 7);

    boole h = loadBoole(src);
    uint64_t orbits;
    if (!value(src, &orbits)) exit(1);
    fclose(src);

    if (ffsize != N) {
        fprintf(stderr, "ffsize=%d, attendu %d\n", (int) ffsize, N);
        exit(1);
    }
    printf("\nUsing no %d  (threads: %d)\n", target, nbproc );
    fflush(stdout);

    uchar *zip = loadZipname(dirname, target, &base);
    vector vend = find_end(zip, base.size, orbits);


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

    #pragma omp parallel for schedule(dynamic, 4096) reduction(+:count)
    for (vector v = 0; v < vend; v++) {
        if (!((zip[v / 8] >> (v % 8)) & 1)) continue;
        count++;

        /* f + h directement sous forme de 2 mots de 64 bits */
        uint64_t tt[2] = { h64[0], h64[1] };
        for (int j = 0; j < nbl; j++)
            if (BITVECTOR(j, v)) {
                tt[0] ^= fct64[j][0];
                tt[1] ^= fct64[j][1];
            }

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

    double temps = omp_get_wtime() - t0;
    sprintf(name, "results/repres-%d.txt", target);
    FILE *dst = stdout;
    if (  optwrite ) 
    	dst = fopen( name, "a" );;
    /* --- Affichage séquentiel, dans l'ordre des v (comme la version série) --- */
    fprintf(dst, "\nUsing no %d  (threads: %d)\n", target, nbproc );
    fprintf(dst, "\n#%zu members have NL_2=38 :\n", nsol);
    qsort(sol, nsol, sizeof *sol, cmp_u64);
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
