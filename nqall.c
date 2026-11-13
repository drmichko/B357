// nqall.c (= nq7_all.c + distribution de NL(f+Q)) : fichier ENTIER de fonctions de B(7) (ANF a..g, une par ligne, préfixe "anf=" accepté), OpenMP.
// Pour chaque f : NQ(f) = distance à RM(2,7), leaders (mots f+q de poids NQ, q dans RM(2,7)), couverture exacte.
// Méthode : q = Q + affine, Q parcourt les 2^21 quadratiques homogènes ; pour chaque Q, transformée de Walsh de f+Q :
//   min_affine wt(f+Q+l) = 64 - max|W|/2 ; leaders de la classe f+Q+RM(1,7) = les fréquences a où |W(a)| est maximal.
// Compilation : gcc -O3 -march=native -fopenmp -o nq7_all nq7_all.c
// Usage       : ./nqall entree.txt sortie.txt [distributions.txt]
// sortie.txt : 3 lignes par f, dans l'ordre du fichier :
//   anf=...
//   NQ= <distance à RM(2,7)> COV= <nombre de points (sur 128) recouverts par l'union des supports des leaders>
//   nl1:n1 nl2:n2 ...   (distribution de NL(f+Q) sur les 2^21 Q homogènes ; NL croissant, effectifs non nuls)
// distributions.txt (optionnel) : détail complet par fonction (format ci-dessous).
// Format du fichier optionnel, une ligne par f :
//   indice deg NQ nb_classes nb_leaders nb_classes_multi points_couverts  nl1:n1 nl2:n2 ...
//   et en plus la distribution de NL(f+Q) = 64 - sup_l|W|/2 sur les 2^21 Q homogènes (effectifs non nuls, NL croissant).
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
  char **out = calloc(n, sizeof(char *)), **out2 = calloc(n, sizeof(char *)); long fait = 0;
  #pragma omp parallel for schedule(dynamic, 4)
  for (long idx = 0; idx < n; idx++) {
    uint8_t tt[128]; int deg; anf2tt(ln[idx], tt, &deg);
    i16 g[128], w[128]; for (int x = 0; x < 128; x++) g[x] = tt[x] ? -1 : 1;
    uint32_t hist[65]; memset(hist, 0, sizeof hist);
    int best = -1; uint64_t cov[2] = {0, 0}; long nl = 0, ncl = 0, nmulti = 0;
    for (uint32_t t = 0; t < (1u << 21); t++) {
      if (t) { int b = __builtin_ctz(t); for (int x = 0; x < 128; x++) g[x] *= qs[b][x]; }
      memcpy(w, g, sizeof w); wht(w);
      int m = 0; for (int a = 0; a < 128; a++) { int v = w[a] < 0 ? -w[a] : w[a]; if (v > m) m = v; }
      hist[64 - m / 2]++;
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
    char s[65 * 16 + 160];
    int pts = __builtin_popcountll(cov[0]) + __builtin_popcountll(cov[1]);
    int len = snprintf(s, sizeof s, "%ld %d %d %ld %ld %ld %d", idx, deg, 64 - best / 2, ncl, nl, nmulti, pts);
    for (int v = 0; v <= 64; v++) if (hist[v]) len += sprintf(s + len, " %d:%u", v, hist[v]);
    s[len++] = '\n'; s[len] = 0;
    out2[idx] = strdup(s);
    { const char *a = ln[idx]; size_t L = strlen(a); while (L && (a[L-1]=='\r')) L--;
      const char *pre = strncmp(a, "anf=", 4) ? "anf=" : "";
      char *o = malloc(L + 65 * 16 + 96); int k2 = sprintf(o, "%s%.*s\nNQ= %d COV= %d\n", pre, (int)L, a, 64 - best / 2, pts);
      for (int v = 0, first = 1; v <= 64; v++) if (hist[v]) { k2 += sprintf(o + k2, first ? "%d:%u" : " %d:%u", v, hist[v]); first = 0; }
      o[k2++] = '\n'; o[k2] = 0; out[idx] = o; }
    long d;
    #pragma omp atomic capture
    d = ++fait;
    if (d % 1000 == 0) { fprintf(stderr, "\r%ld / %ld", d, n); fflush(stderr); }
  }
  FILE *G = fopen(argv[2], "w"); if (!G) { perror("sortie"); return 1; }
  for (long i = 0; i < n; i++) fputs(out[i], G);
  fclose(G);
  if (argc > 3) { FILE *H = fopen(argv[3], "w"); if (!H) { perror("distributions"); return 1; }
    for (long i = 0; i < n; i++) fputs(out2[i], H); fclose(H); } fprintf(stderr, "\nterminé : %ld fonctions\n", n);
  return 0;
}
