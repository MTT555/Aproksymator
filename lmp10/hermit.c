//jako baze funkcji wybieramy wielomian Hermite'a, ponizszy kod odpowiada za obliczenie wartosci takiego wielomianu
#include "makespl.h"
#include "piv_ge_solver.h"
#include "rozwiazywacz.h"
#include "gaus/matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>




double hermit(int n,double x){

//obliczamy rekurencyjnie ze wzoru
if (n == 0) return 1;
else if (n == 1) return 2*x;
else return 2*x*hermit(n-1, x)-2*(n-1)*hermit(n-2,x);

}
//kolejne pochodne, pierwsza druga i trzecia
double hermit1(int n, double x){

if (n == 0) return 0;
else if (n == 1) return 2;
else return 2*hermit(n-1, x) + 2*x*hermit1(n-1,x)-2*(n-1)*hermit1(n-2,x);


}

double hermit2(int n, double x){

if (n == 0) return 0;
else if (n == 1) return 0;
else return 4*hermit1(n-1, x) + 2*x*hermit2(n-1,x)-2*(n-1)*hermit2(n-2,x);


}

double hermit3(int n, double x){

if (n == 0) return 0;
else if (n == 1) return 0;
else return 6*hermit2(n-1, x) + 2*x*hermit3(n-1,x)-2*(n-1)*hermit3(n-2,x);


}
//wartosc funkcji
double f(double *wspolczynniki, double xx, double nb){
double sumuj = 0;
for (int i = 0; i < nb; i++){
sumuj += wspolczynniki[i]*hermit(i,xx);
}
return sumuj;
}
// pierwsza pochodna
double f1(double *wspolczynniki, double xx, double nb){
double sumuj = 0;
for (int i = 0; i < nb; i++){
sumuj += wspolczynniki[i]*hermit1(i,xx);
}
return sumuj;

}
//druga pochodna
double f2(double *wspolczynniki, double xx, double nb){
double sumuj = 0;
for (int i = 0; i < nb; i++){
sumuj += wspolczynniki[i]*hermit2(i,xx);
}
return sumuj;

}
//trzecia pochodna
double f3(double *wspolczynniki, double xx, double nb){
double sumuj = 0;
for (int i = 0; i < nb; i++){
sumuj += wspolczynniki[i]*hermit3(i,xx);
}
return sumuj;

}
//kod odpowiedzialny za aproksymacje sredniokwadratowa z baza wielomianow Hermita

void make_spl(points_t *pts, spline_t *spl){
matrix_t *macierz = NULL; //macierz przechowujaca uklad rownan


// ile chcemy miec funkcji bazowych?
int nb = pts->n - 3 > 10 ? 10 : pts->n - 3;

  char *nbEnv= getenv( "APPROX_BASE_SIZE" ); // pobieramy ze zmiennej srodowiskowej

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	macierz = make_matrix(nb, nb + 1); // tworzymy macierz


printf("nb: %d\n",nb);
printf("Punktow jest: %d\n", pts->n);
//pobieram sobie wektory z punktami i ilosc tych punktow
double *x = pts->x;
double *y = pts->y;
double n = pts->n;

// a moze uloze sobie taka pomocnicza macierz zeby sie latwiej liczylo
double *pom[nb];
for (int i = 0; i < nb; i++) pom[i] = malloc(sizeof(double) * (nb+1));
// wypelniam zerami bo trzeba bedzie potem sumowac
for (int i = 0; i < nb; i++){
for (int j = 0; j < nb+1; j++){
 pom[i][j] = 0;
}
}

//ukladam uklad rownan - najpierw czesc macierzy poza kolumna najbardziej po prawej

double a = 0, b = 0; // a,b - kolejno indeksy funkcji bazowych

double *wskaznik = pts->x;
for (int i = 0; i < nb; i++){
for (int j = 0; j < nb; j++){
 for (int k = 0; k < n; k++){
pom[i][j] += hermit(a,*wskaznik)*hermit(b,*wskaznik);
wskaznik++;
 }
wskaznik = pts->x;
a++;
}
b++;
a = 0;
wskaznik = pts->x;
}

a = 0;
wskaznik = pts->x;
double *wskaznik_y = pts->y;
// teraz kolumna po prawej
for (int i = 0; i < nb; i++){
for (int k = 0; k < n; k++){
pom[i][nb] += hermit(a, *wskaznik)*(*wskaznik_y);
wskaznik++;
wskaznik_y++;

}
a++;
wskaznik = pts->x;
wskaznik_y = pts->y;
}

// dodajemy dane do macierzy matrix_t

for (int i = 0; i < nb; i++){
for (int j = 0; j < nb+1; j++){
put_entry_matrix(macierz, i,j,pom[i][j]);

}
}
write_matrix(macierz,stdout);
rozwiaz(macierz);

//pobieram wspolczynniki z pliku

FILE *in = fopen("wspolczynniki","r");
if (in == NULL) printf("Blad ladowania pliku ze wspolczynnikami\n");
double pom1;
double *wspolczynniki = malloc(sizeof(double)*nb);
int iterator = 0;
while (fscanf(in, "%lg", &pom1) != EOF){
wspolczynniki[iterator] = pom1;
iterator++;
}
fclose(in);
//pobieram pierwszy i ostatni punkt kolejno do zmiennych a i b
a = x[0];
b = x[pts->n-1];
if (alloc_spl(spl, nb) == 0) {
		for (int i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] =  f(wspolczynniki, xx, nb);
			spl->f1[i] = f1(wspolczynniki, xx, nb);
			spl->f2[i] = f2(wspolczynniki, xx, nb);
			spl->f3[i] = f3(wspolczynniki, xx, nb);
			
		}
	}
	for (int i = 0; i < nb; i++)free(pom[i]);
	free(wspolczynniki);
	free_matrix(macierz);

}
