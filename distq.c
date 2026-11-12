// distq_all.c : traite un fichier ENTIER de fonctions de B(7) (ANF en lettres a..g, une par ligne, préfixe "anf=" accepté).
// Pour chaque f : distribution de wt(f+q), q parcourant les 2^21 formes quadratiques homogènes en 7 variables.
// Parallélisé avec OpenMP (tous les cœurs disponibles, ou OMP_NUM_THREADS=n).
// Compilation : gcc -O3 -march=native -fopenmp -o distq_all distq_all.c
// Usage       : ./distq_all entree.txt sortie.txt
// Sortie      : une ligne par f, dans l'ordre du fichier :  indice  w1:n1 w2:n2 ...  (effectifs non nuls)
// Variables   : a=x1 (bit 0) ... g=x7 (bit 6) ; table de vérité indexée par x = sum x_i 2^(i-1).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static void anf2tt(const char *p, uint64_t tt[2]) {
  uint8_t c[128]; memset(c, 0, 128);
  if (!strncmp(p, "anf=", 4)) p += 4;
  while (*p && *p != '\n' && *p != '\r') {
    int m = 0;
    while (*p && *p != '+' && *p != '\n' && *p != '\r') {
      if (*p >= 'a' && *p <= 'g') m |= 1 << (*p - 'a');
      p++;
    }
    c[m] ^= 1;
    if (*p == '+') p++;
  }
  for (int i = 0; i < 7; i++)
    for (int x = 0; x < 128; x++) if (x >> i & 1) c[x] ^= c[x ^ (1 << i)];
  tt[0] = tt[1] = 0;
  for (int x = 0; x < 128; x++) if (c[x]) tt[x >> 6] |= 1ULL << (x & 63);
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s entree.txt sortie.txt\n", argv[0]); return 1; }
  FILE *F = fopen(argv[1], "rb");
  if (!F) { perror("entree"); return 1; }
  fseek(F, 0, SEEK_END); long sz = ftell(F); fseek(F, 0, SEEK_SET);
  char *buf = malloc(sz + 2);
  if (fread(buf, 1, sz, F) != (size_t)sz) { fprintf(stderr, "lecture incomplète\n"); return 1; }
  buf[sz] = '\n'; buf[sz + 1] = 0; fclose(F);

  long n = 0, cap = 1 << 20; char **ln = malloc(cap * sizeof(char *));   // début de chaque ligne non vide
  for (char *p = buf; p < buf + sz + 1;) {
    char *e = strchr(p, '\n'); if (!e) break;
    if (e > p && !(e == p + 1 && *p == '\r')) { if (n == cap) { cap *= 2; ln = realloc(ln, cap * sizeof(char *)); } ln[n++] = p; }
    *e = 0; p = e + 1;
  }

  uint64_t mq[21][2]; int k = 0;
  for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++) {
    mq[k][0] = mq[k][1] = 0;
    for (int x = 0; x < 128; x++) if ((x >> i & 1) && (x >> j & 1)) mq[k][x >> 6] |= 1ULL << (x & 63);
    k++;
  }

  char **out = calloc(n, sizeof(char *));
  long fait = 0;
  #pragma omp parallel for schedule(dynamic, 16)
  for (long i = 0; i < n; i++) {
    uint64_t tt[2]; anf2tt(ln[i], tt);
    uint64_t lo = tt[0], hi = tt[1];
    uint32_t h[129]; memset(h, 0, sizeof h);
    h[__builtin_popcountll(lo) + __builtin_popcountll(hi)]++;
    for (uint32_t t = 1; t < (1u << 21); t++) {
      int b = __builtin_ctz(t);
      lo ^= mq[b][0]; hi ^= mq[b][1];
      h[__builtin_popcountll(lo) + __builtin_popcountll(hi)]++;
    }
    char s[129 * 16 + 32]; int len = sprintf(s, "%ld", i);
    for (int w = 0; w <= 128; w++) if (h[w]) len += sprintf(s + len, " %d:%u", w, h[w]);
    s[len++] = '\n'; s[len] = 0;
    out[i] = strdup(s);
    long d;
    #pragma omp atomic capture
    d = ++fait;
    if (d % 5000 == 0) { fprintf(stderr, "\r%ld / %ld", d, n); fflush(stderr); }
  }
  FILE *G = fopen(argv[2], "w");
  if (!G) { perror("sortie"); return 1; }
  for (long i = 0; i < n; i++) fputs(out[i], G);
  fclose(G);
  fprintf(stderr, "\nterminé : %ld fonctions\n", n);
  return 0;
}
