// nq7_all.c : fichier ENTIER de fonctions de B(7) (ANF a..g, une par ligne, préfixe "anf=" accepté), OpenMP.
// Pour chaque f : NQ(f) = distance à RM(2,7), leaders (mots f+q de poids NQ, q dans RM(2,7)), couverture exacte.
// Méthode : q = Q + affine, Q parcourt les 2^21 quadratiques homogènes ; pour chaque Q, transformée de Walsh de f+Q :
//   min_affine wt(f+Q+l) = 64 - max|W|/2 ; leaders de la classe f+Q+RM(1,7) = les fréquences a où |W(a)| est maximal.
// Compilation : gcc -O3 -march=native -fopenmp -o nq7_all nq7_all.c
// Usage       : ./nq7_all entree.txt sortie.txt
// Sortie (une ligne par f, dans l'ordre) :  indice deg NQ nb_classes nb_leaders nb_classes_multi points_couverts
//   nb_classes       = nombre de Q (classes modulo RM(1,7)) contenant un leader
//   nb_leaders       = nombre total de leaders (mots de poids NQ dans f+RM(2,7))
//   nb_classes_multi = nombre de classes contenant PLUS d'un leader (0 si "au plus un leader modulo RM(1,7)" est vrai)
//   points_couverts  = taille de l'union des supports des leaders (128 : recouvrement total)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef int16_t i16;
static void anf2tt(const char *p, uint8_t *tt, int *deg) {
  uint8_t c[128]; memset(c, 0, 128); *deg = 0;
  if (!strncmp(p, "anf=", 4)) p += 4;
  while (*p && *p != '\n' && *p != '\r') {
    int m = 0, d = 0;
    while (*p && *p != '+' && *p != '\n' && *p != '\r') { if (*p >= 'a' && *p <= 'g') { m |= 1 << (*p - 'a'); d++; } p++; }
    c[m] ^= 1; if (d > *deg) *deg = d;
    if (*p == '+') p++;
  }
  for (int i = 0; i < 7; i++) for (int x = 0; x < 128; x++) if (x >> i & 1) c[x] ^= c[x ^ (1 << i)];
  memcpy(tt, c, 128);
}
static inline void wht(i16 *a) {
  for (int h = 1; h < 128; h <<= 1)
    for (int i = 0; i < 128; i += h << 1)
      for (int j = i; j < i + h; j++) { i16 u = a[j], v = a[j + h]; a[j] = u + v; a[j + h] = u - v; }
}
int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s entree.txt sortie.txt\n", argv[0]); return 1; }
  FILE *F = fopen(argv[1], "rb"); if (!F) { perror("entree"); return 1; }
  fseek(F, 0, SEEK_END); long sz = ftell(F); fseek(F, 0, SEEK_SET);
  char *buf = malloc(sz + 2);
  if (fread(buf, 1, sz, F) != (size_t)sz) { fprintf(stderr, "lecture incomplète\n"); return 1; }
  buf[sz] = '\n'; buf[sz + 1] = 0; fclose(F);
  long n = 0, cap = 1 << 20; char **ln = malloc(cap * sizeof(char *));
  for (char *p = buf; p < buf + sz + 1;) {
    char *e = strchr(p, '\n'); if (!e) break;
    if (e > p && !(e == p + 1 && *p == '\r')) { if (n == cap) { cap *= 2; ln = realloc(ln, cap * sizeof(char *)); } ln[n++] = p; }
    *e = 0; p = e + 1;
  }
  static i16 qs[21][128]; int k = 0;
  for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++) { for (int x = 0; x < 128; x++) qs[k][x] = ((x >> i & 1) & (x >> j & 1)) ? -1 : 1; k++; }
  char **out = calloc(n, sizeof(char *)); long fait = 0;
  #pragma omp parallel for schedule(dynamic, 4)
  for (long idx = 0; idx < n; idx++) {
    uint8_t tt[128]; int deg; anf2tt(ln[idx], tt, &deg);
    i16 g[128], w[128]; for (int x = 0; x < 128; x++) g[x] = tt[x] ? -1 : 1;
    int best = -1; uint64_t cov[2] = {0, 0}; long nl = 0, ncl = 0, nmulti = 0;
    for (uint32_t t = 0; t < (1u << 21); t++) {
      if (t) { int b = __builtin_ctz(t); for (int x = 0; x < 128; x++) g[x] *= qs[b][x]; }
      memcpy(w, g, sizeof w); wht(w);
      int m = 0; for (int a = 0; a < 128; a++) { int v = w[a] < 0 ? -w[a] : w[a]; if (v > m) m = v; }
      if (m < best) continue;
      if (m > best) { best = m; cov[0] = cov[1] = 0; nl = ncl = nmulti = 0; }
      int cnt = 0;
      for (int a = 0; a < 128; a++) {
        int v = w[a] < 0 ? -w[a] : w[a]; if (v != m) continue;
        cnt++; int target = w[a] > 0 ? -1 : 1;
        for (int x = 0; x < 128; x++) { int s = g[x] * (__builtin_parity(a & x) ? -1 : 1); if (s == target) cov[x >> 6] |= 1ULL << (x & 63); }
      }
      nl += cnt; ncl++; if (cnt > 1) nmulti++;
    }
    char s[160];
    int pts = __builtin_popcountll(cov[0]) + __builtin_popcountll(cov[1]);
    snprintf(s, sizeof s, "%ld %d %d %ld %ld %ld %d\n", idx, deg, 64 - best / 2, ncl, nl, nmulti, pts);
    out[idx] = strdup(s);
    long d;
    #pragma omp atomic capture
    d = ++fait;
    if (d % 1000 == 0) { fprintf(stderr, "\r%ld / %ld", d, n); fflush(stderr); }
  }
  FILE *G = fopen(argv[2], "w"); if (!G) { perror("sortie"); return 1; }
  for (long i = 0; i < n; i++) fputs(out[i], G);
  fclose(G); fprintf(stderr, "\nterminé : %ld fonctions\n", n);
  return 0;
}
