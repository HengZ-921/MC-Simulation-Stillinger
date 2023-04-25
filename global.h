#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define A 7.049556277
#define B 0.6022245584
#define alpha 1.8
#define lambda 21.0
#define gama 1.2
#define sigma 2.1
#define sPi 4.1887902048 // 4Pi/3

typedef struct NBL {
	int NAtom, *num, **atom;
	double cL, Rk, rho, RnSqr, **coord, **v;
} NBL;

double** readxyz ( int* numAtom, double* cell );
void free2D ( int N, double** mtx );
void cleanAll ( int N, NBL nbLi, double** xyz );

NBL initNBL ( int N, double L );
void updateNBL ( double** xyz, NBL nbLi );

double wrap ( double pos, double cell );
void minVec ( double* I, double* J, double* vec, double L );
double dSqr ( double* vec );
double potSingle ( int I, double* Ixyz, NBL li, double** xyz );

int trialMove ( double* trial, double dr, double T, NBL li, double** xyz, double* dtt );
void refreshNBL ( double* dtt, double** xyz, NBL li );
void eqSweep ( double *dr, double* T, double** xyz, NBL li );
void prodSweep ( double dr, double T, double** xyz, NBL li );
