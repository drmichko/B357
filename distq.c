// distq.c : pour chaque f de B(7) (ANF en lettres a..g, une par ligne, préfixe "anf=" accepté),
// distribution de wt(f+q) pour q parcourant les 2^21 formes quadratiques homogènes en 7 variables.
// Usage : ./distq fichier [debut=0] [nombre=tout] > sortie.txt
// Sortie, une ligne par f :  indice  w1:n1 w2:n2 ...   (seulement les poids de effectif non nul, poids croissants)
// Variables : a=x1 (bit 0) ... g=x7 (bit 6) ; la table de vérité est indexée par x = sum x_i 2^(i-1).
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
    c[m] ^= 1;                       // monôme 1 (m=0) = constante
    if (*p == '+') p++;
  }
  for (int i = 0; i < 7; i++)         // transformée de Möbius
    for (int x = 0; x < 128; x++) if (x >> i & 1) c[x] ^= c[x ^ (1 << i)];
  tt[0] = tt[1] = 0;
  for (int x = 0; x < 128; x++) if (c[x]) tt[x >> 6] |= 1ULL << (x & 63);
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s fichier [debut] [nombre]\n", argv[0]); return 1; }
  FILE *F = fopen(argv[1], "r");
  if (!F) { perror("open"); return 1; }
  long start = argc > 2 ? atol(argv[2]) : 0;
  long cnt = argc > 3 ? atol(argv[3]) : (1L << 60);

  uint64_t mq[21][2]; int k = 0;      // tables de vérité des monômes x_i x_j
  for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++) {
    mq[k][0] = mq[k][1] = 0;
    for (int x = 0; x < 128; x++) if ((x >> i & 1) && (x >> j & 1)) mq[k][x >> 6] |= 1ULL << (x & 63);
    k++;
  }
  static char line[1 << 16];
  long idx = 0, done = 0;
  while (fgets(line, sizeof line, F)) {
    if (idx < start) { idx++; continue; }
    if (done >= cnt) break;
    uint64_t lo, hi, tt[2]; anf2tt(line, tt); lo = tt[0]; hi = tt[1];
    uint32_t h[129]; memset(h, 0, sizeof h);
    h[__builtin_popcountll(lo) + __builtin_popcountll(hi)]++;      // q = 0
    for (uint32_t t = 1; t < (1u << 21); t++) {                    // code de Gray
      int b = __builtin_ctz(t);
      lo ^= mq[b][0]; hi ^= mq[b][1];
      h[__builtin_popcountll(lo) + __builtin_popcountll(hi)]++;
    }
    printf("%ld", idx);
    for (int w = 0; w <= 128; w++) if (h[w]) printf(" %d:%u", w, h[w]);
    printf("\n");
    idx++; done++;
  }
  fclose(F);
  return 0;
}
