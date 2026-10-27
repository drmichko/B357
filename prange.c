/*
 * nlq7_omp.c : test NQ(f) >= T pour des fonctions booléennes à 7 variables,
 * version parallèle OpenMP, avec pré-filtre de type Prange.
 * Table de vérité : f(x) = (tt[x >> 6] >> (x & 63)) & 1.
 *
 * Compilation : gcc -O3 -march=native -fopenmp -o nlq7_omp nlq7_omp.c [+ vos autres .c]
 *               (il faut au minimum -mpopcnt -mbmi2 pour popcnt et pext)
 * Exécution   : OMP_NUM_THREADS=16 ./nlq7_omp -t <target> -g <goal>
 *               variable d'environnement NPS = nombre d'ensembles Prange
 *               (défaut 1024, 0 = filtre désactivé)
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
#include <assert.h>
#include <immintrin.h>
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

/* ====================================================================== */
/*  Test exact : NQ(f) >= T ?  (Gray déroulé par 4, sortie anticipée)    */
/*  Entrée : c0,c1 = table de vérité de f (2 mots de 64 bits).            */
/*  Échec dès qu'une g de degré <= 2 vérifie d(f,g) < T ou > N - T.       */
/* ====================================================================== */
static inline int nq_at_least_m7(uint64_t c0, uint64_t c1, int T)
{
    const unsigned lo   = (unsigned)T;
    const unsigned span = (unsigned)(N - 2 * T);
#define CHK() if ((unsigned)(__builtin_popcountll(c0) + __builtin_popcountll(c1)) - lo > span) return 0
#define FLIP(k) do { c0 ^= mono[k][0]; c1 ^= mono[k][1]; } while (0)
    uint32_t m = 1;                          /* indice de bloc de 4 états */
    for (;;) {
        CHK(); FLIP(0);
        CHK(); FLIP(1);
        CHK(); FLIP(0);
        CHK();
        if (m >> (NMONO - 2)) return 1;      /* 2^28 états parcourus */
        FLIP(2 + __builtin_ctz(m));
        m++;
    }
#undef CHK
#undef FLIP
}

/* ====================================================================== */
/*  Filtre Prange : ensembles d'information de RM(2,7)                    */
/*  Un témoin q de RM(2,7) avec wt(tt+q) < T prouve NQ < T (rejet sûr).   */
/* ====================================================================== */
typedef struct {
    uint64_t m0, m1;          /* masque des 29 pivots dans tt[0], tt[1] */
    int      n0;              /* popcount(m0) */
    uint64_t row[29][2];      /* ligne RREF de chaque pivot, par position croissante */
} pset;

static pset *PS  = NULL;
static int   nps = 0;

static uint64_t pr_state = 0x9E3779B97F4A7C15ULL;
static uint64_t pr_rng(void)                 /* splitmix64 */
{
    uint64_t z = (pr_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

#define GETB(Mx, r, c) (((Mx)[r][(c) >> 6] >> ((c) & 63)) & 1)

/* Matrice génératrice de RM(2,7) : 1, x_i, x_i x_j (1+7+21 = 29 lignes) */
static void pr_build_G(uint64_t G[29][2])
{
    memset(G, 0, 29 * sizeof G[0]);
    int r = 0;
    for (int a = 0; a < 128; a++) G[r][a >> 6] |= 1ULL << (a & 63);
    r++;
    for (int i = 0; i < 7; i++, r++)
        for (int a = 0; a < 128; a++)
            if ((a >> i) & 1) G[r][a >> 6] |= 1ULL << (a & 63);
    for (int i = 0; i < 7; i++)
        for (int j = i + 1; j < 7; j++, r++)
            for (int a = 0; a < 128; a++)
                if (((a >> i) & (a >> j)) & 1) G[r][a >> 6] |= 1ULL << (a & 63);
}

/* Permutation aléatoire des colonnes + Gauss-Jordan -> forme systématique */
static void pr_make_pset(pset *s, const uint64_t G[29][2])
{
    uint64_t Mx[29][2];
    memcpy(Mx, G, sizeof Mx);

    int perm[128];
    for (int i = 0; i < 128; i++) perm[i] = i;
    for (int i = 127; i > 0; i--) {                 /* Fisher-Yates */
        int j = (int)(pr_rng() % (uint64_t)(i + 1));
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }

    int piv[29], rank = 0;
    for (int k = 0; k < 128 && rank < 29; k++) {
        int c = perm[k];
        int r = rank;
        while (r < 29 && !GETB(Mx, r, c)) r++;
        if (r == 29) continue;                      /* colonne dépendante */
        if (r != rank) {
            uint64_t t0 = Mx[r][0], t1 = Mx[r][1];
            Mx[r][0] = Mx[rank][0]; Mx[r][1] = Mx[rank][1];
            Mx[rank][0] = t0;       Mx[rank][1] = t1;
        }
        for (int i = 0; i < 29; i++)
            if (i != rank && GETB(Mx, i, c)) {
                Mx[i][0] ^= Mx[rank][0];
                Mx[i][1] ^= Mx[rank][1];
            }
        piv[rank++] = c;
    }
    assert(rank == 29);

    /* pext renvoie les bits par position croissante : on trie les lignes par pivot */
    for (int i = 1; i < 29; i++) {
        int pc = piv[i]; uint64_t a0 = Mx[i][0], a1 = Mx[i][1];
        int j = i - 1;
        while (j >= 0 && piv[j] > pc) {
            piv[j + 1] = piv[j];
            Mx[j + 1][0] = Mx[j][0]; Mx[j + 1][1] = Mx[j][1];
            j--;
        }
        piv[j + 1] = pc; Mx[j + 1][0] = a0; Mx[j + 1][1] = a1;
    }

    s->m0 = s->m1 = 0;
    for (int i = 0; i < 29; i++) {
        if (piv[i] < 64) s->m0 |= 1ULL << piv[i];
        else             s->m1 |= 1ULL << (piv[i] - 64);
        s->row[i][0] = Mx[i][0];
        s->row[i][1] = Mx[i][1];
    }
    s->n0 = __builtin_popcountll(s->m0);
}

/* Retourne 1 si un témoin q prouve NQ(tt) < T (rejet sûr), 0 sinon. */
static inline int prange_reject(uint64_t t0, uint64_t t1, int T)
{
    for (int r = 0; r < nps; r++) {
        const pset *s = &PS[r];
        uint64_t b = _pext_u64(t0, s->m0) | (_pext_u64(t1, s->m1) << s->n0);
        uint64_t c0 = t0, c1 = t1;                  /* on accumule tt + q */
        while (b) {
            int p = __builtin_ctzll(b);
            b &= b - 1;
            c0 ^= s->row[p][0];
            c1 ^= s->row[p][1];
        }
        if (__builtin_popcountll(c0) + __builtin_popcountll(c1) < T) return 1;
    }
    return 0;
}

/* ====================================================================== */
/*  Tables par tranches de 8 bits : tt = h + sum_j v_j fct64[j]           */
/* ====================================================================== */
#define NSL_MAX 8
static uint64_t TAB[NSL_MAX][256][2];
static int nsl;

static void build_tab(int nbl, const uint64_t (*fct64)[2])
{
    assert(nbl <= 8 * NSL_MAX);
    nsl = (nbl + 7) / 8;
    for (int k = 0; k < nsl; k++)
        for (int x = 0; x < 256; x++) {
            uint64_t a0 = 0, a1 = 0;
            for (int t = 0; t < 8; t++) {
                int j = 8 * k + t;
                if (j < nbl && ((x >> t) & 1)) { a0 ^= fct64[j][0]; a1 ^= fct64[j][1]; }
            }
            TAB[k][x][0] = a0; TAB[k][x][1] = a1;
        }
}

static inline void tt_of(uint64_t v, const uint64_t h64[2], uint64_t *t0, uint64_t *t1)
{
    uint64_t a = h64[0], b = h64[1];
    for (int k = 0; k < nsl; k++) {
        const uint64_t *e = TAB[k][(v >> (8 * k)) & 255];
        a ^= e[0]; b ^= e[1];
    }
    *t0 = a; *t1 = b;
}

static inline uint64_t zip_word(const uint8_t *zip, uint64_t b, uint64_t nbytes)
{
    uint64_t w = 0;                               /* octets absents => 0 */
    uint64_t off = 8 * b;
    size_t n = (nbytes - off < 8) ? (size_t)(nbytes - off) : 8;
    memcpy(&w, zip + off, n);                     /* little-endian */
    return w;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
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

    /* --- Tables par tranches + ensembles Prange (précalcul, hors région parallèle) --- */
    build_tab(nbl, (const uint64_t (*)[2]) fct64);

    {   /* contrôle : tables par tranches == calcul direct avec BITVECTOR */
        uint64_t nb = 8 * (last + 1);
        for (int s = 0; s < 1000; s++) {
            uint64_t vv = pr_rng() % nb;
            uint64_t r0 = h64[0], r1 = h64[1], a0, a1;
            for (int j = 0; j < nbl; j++)
                if (BITVECTOR(j, (vector) vv)) { r0 ^= fct64[j][0]; r1 ^= fct64[j][1]; }
            tt_of(vv, h64, &a0, &a1);
            if (a0 != r0 || a1 != r1) {
                fprintf(stderr, "erreur : BITVECTOR(j,v) != (v>>j)&1, adapter tt_of\n");
                exit(1);
            }
        }
    }

    nps = 1024;
    { const char *e = getenv("NPS"); if (e) nps = atoi(e); }
    if (nps > 0) {
        uint64_t G[29][2];
        PS = calloc((size_t) nps, sizeof *PS);
        if (!PS) { perror("PS"); exit(1); }
        pr_state = 12345;                          /* graine fixe : reproductible */
        pr_build_G(G);
        for (int r = 0; r < nps; r++) pr_make_pset(&PS[r], G);
    }

    /* --- Solutions : tampons locaux par thread, fusion sous critical --- */
    size_t cap = 0, nsol = 0;
    uint64_t *sol = NULL;
    int err = 0;

    uint64_t count = 0, nexact = 0;
    double t0 = omp_get_wtime();

    const uint64_t nbytes =  last + 1;
    const uint64_t nwords = (nbytes + 7) / 8;
    const int g = goal;

#pragma omp parallel reduction(+:count, nexact)
    {
        uint64_t *lsol = NULL;
        size_t lnsol = 0, lcap = 0;

#pragma omp for schedule(dynamic, 64) nowait
        for (uint64_t b = 0; b < nwords; b++) {
            uint64_t w = zip_word((const uint8_t *)zip, b, nbytes);

            while (w) {
                uint64_t v = 64 * b + (uint64_t) __builtin_ctzll(w);
                w &= w - 1;
                count++;

                uint64_t a0, a1;
                tt_of(v, h64, &a0, &a1);

                if (prange_reject(a0, a1, g)) continue;
                nexact++;
                if (!nq_at_least_m7(a0, a1, g)) continue;

                if (lnsol == lcap) {
                    lcap = lcap ? 2 * lcap : 256;
                    uint64_t *p = realloc(lsol, lcap * sizeof *p);
                    if (!p) { err = 1; break; }
                    lsol = p;
                }
                lsol[lnsol++] = v;
            }
        }

        if (lnsol) {
#pragma omp critical(sol)
            {
                if (!err && nsol + lnsol > cap) {
                    size_t nc = cap ? cap : 1024;
                    while (nsol + lnsol > nc) nc *= 2;
                    uint64_t *p = realloc(sol, nc * sizeof *p);
                    if (!p) err = 1; else { sol = p; cap = nc; }
                }
                if (!err) { memcpy(sol + nsol, lsol, lnsol * sizeof *lsol); nsol += lnsol; }
            }
        }
        free(lsol);
    }
    if (err) { perror("realloc"); exit(1); }
    if (nsol) qsort(sol, nsol, sizeof *sol, cmp_u64);   /* sortie reproductible */

    printf("candidat:%zu\n", nsol );
    double temps = omp_get_wtime() - t0;
    char name[128];
    sprintf(name, "results/res-%d.txt", target);
    FILE *dst = stdout;
    if ( optw ) dst= fopen( name, "a" );
    if ( !dst ) { perror(name); dst = stdout; }
    /* --- Affichage séquentiel, dans l'ordre des v (comme la version série) --- */
    fprintf(dst, "\nUsing no %d  (threads: %d)\n", target, nbproc );
    fprintf(dst, "\n#%zu members have NL_2 >=%d :\n", nsol, goal);
    for (size_t i = 0; i < nsol; i++) {
        boole f = vectortoboole((vector) sol[i], &base);
        for (int x = 0; x < ffsize; x++) f[x] ^= h[x];
        panf( dst, f);
        free(f);
    }
    fflush(dst);

    fprintf(dst, "\n#orbit=%" PRIu64 " count=%" PRIu64 "  %zu\n", orbits, count, nsol);
    fprintf(dst, "#prange: nps=%d  tests exacts=%" PRIu64 " (%.4f%% des candidats)\n",
            nps, nexact, count ? 100.0 * (double) nexact / (double) count : 0.0);
    fprintf(dst, "\n#runtime : %.6f secondes (mur)\n", temps);

    if (dst != stdout) fclose( dst );
    free(sol);
    free(PS);
    free(fct64);
    free(zip);

    return 0;
}
