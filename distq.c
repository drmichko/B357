// nldist_all.c : fichier ENTIER de fonctions de B(7) (ANF a..g, une par ligne, préfixe "anf=" accepté), OpenMP.
// Pour chaque f : distribution de NL(f+q), q parcourant les 2^21 formes quadratiques homogènes en 7 variables,
// NL(g) = 64 - max_a |W_g(a)| / 2 (distance de g à RM(1,7)).
// Compilation : gcc -O3 -march=native -fopenmp -o nldist_all nldist_all.c
// Usage       : ./nldist_all entree.txt sortie.txt
// Sortie      : une ligne par f, dans l'ordre :  indice  nl1:n1 nl2:n2 ...  (valeurs de NL d'effectif non nul, croissantes)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef int16_t i16;
static void anf2tt(const char *p, uint8_t *tt) {
  uint8_t c[128]; memset(c, 0, 128);
  if (!strncmp(p, "anf=", 4)) p += 4;
  while (*p && *p != '\n' && *p != '\r') {
    int m = 0;
    while (*p && *p != '+' && *p != '\n' && *p != '\r') { if (*p >= 'a' && *p <= 'g') m |= 1 << (*p - 'a'); p++; }
    c[m] ^= 1;
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
    uint8_t tt[128]; anf2tt(ln[idx], tt);
    i16 g[128], w[128]; for (int x = 0; x < 128; x++) g[x] = tt[x] ? -1 : 1;
    uint32_t h[65]; memset(h, 0, sizeof h);
    for (uint32_t t = 0; t < (1u << 21); t++) {
      if (t) { int b = __builtin_ctz(t); for (int x = 0; x < 128; x++) g[x] *= qs[b][x]; }
      memcpy(w, g, sizeof w); wht(w);
      int m = 0; for (int a = 0; a < 128; a++) { int v = w[a] < 0 ? -w[a] : w[a]; if (v > m) m = v; }
      h[64 - m / 2]++;
    }
    char s[65 * 16 + 32]; int len = sprintf(s, "%ld", idx);
    for (int v = 0; v <= 64; v++) if (h[v]) len += sprintf(s + len, " %d:%u", v, h[v]);
    s[len++] = '\n'; s[len] = 0;
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
