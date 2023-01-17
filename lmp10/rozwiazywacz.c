#include "gaus/matrix.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void rozwiaz(matrix_t *macierz)
{

//rozwiazania zapisze jednoczesnie do pliku zeby mogl je pobierac kod z pliku hermit.c
FILE *out = fopen("wspolczynniki","w");
    printf ("Rozwiązania: \n");
    matrix_t *m = macierz;
    int sym = 0;
    if (m != NULL) {
      matrix_t *c = NULL;
      int *row_per = malloc (m->rn * sizeof *row_per);
      printf ("\nMacierz:\n");
      write_matrix (m, stdout);
      c = pivot_ge_matrix (m, row_per);
        
     if (c != NULL) {
        int i;
       // printf ("\nPo elim. Gaussa:\n");
       // write_matrix (c, stdout);
       // printf ("Permutacja:");
        for (i = 0; i < c->rn; i++)
        //  printf (" %d", row_per[i]);
       // printf ("\n");
        if (bs_matrix (c) == 0) {
          int j;
          int *iper = pivot_get_inv_per (c, row_per);
         // printf ("Permutacja odwrotna:");
          for (i = 0; i < c->rn; i++)
         //   printf (" %d", iper[i]);
         // printf ("\n");
         // printf ("\nPo podstawieniu wstecz:\n");
         // write_matrix (c, stdout);
        //  printf ("Rozwiązania: \n");
          for (j = 0; j < c->cn - c->rn; j++) {
            printf ("Nr %d:\n", j + 1);
            for (i = 0; i < c->rn; i++) {
              int oi = sym ? iper[i] : i;
              printf ("\t%f", *(c->e + oi * c->cn + j + c->rn));
              fprintf(out, "%f\n", *(c->e + oi * c->cn + j + c->rn));
            }
            printf ("\n");
          }
        }
      }
      else
        printf ("nie lula!\n");
    }
    fclose(out);
 
}
