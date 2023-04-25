#include "global.h"

NBL initNBL ( int N, double L ) {
	int numNB;
	double Rc, Rk, Rn;
	NBL li;
	scanf( "%*s %lf\n", &Rc );
	scanf( "%*s %lf\n", &Rk );
	Rn = Rc + Rk;
	li.NAtom = N;
	li.cL = L;
	li.Rk = Rk;
	li.RnSqr = Rn * Rn;
	li.rho = N * pow( L, -3.0 );
	numNB = 2 * li.rho * sPi * pow( Rn, 3.0 );
	li.num = calloc( N, sizeof(int) );
	li.atom = malloc( N * sizeof(int*) );
	li.coord = malloc( N * sizeof(double*) );
	li.v = malloc( 3 * sizeof(double*) );
	for ( int i = 0; i < N; i++ ) {
		li.atom[i] = malloc( numNB * sizeof(int) );
		li.coord[i] = malloc( 3 * sizeof(double) );
	}
	for ( int i = 0; i < 3; i++ ) {
		li.v[i] = calloc( 3, sizeof(double) );
	}
	return li;
}

void updateNBL ( double** xyz, NBL li ) {
	for ( int i = 0; i < li.NAtom; i++ ) { li.num[i] = 0; }
	for ( int i = 0; i < li.NAtom; i++ ) {
		for ( int k = 0; k < 3; k++ ) {
			li.coord[i][k] = xyz[i][k];
		}
		for ( int j = i + 1; j < li.NAtom; j++ ) {
			minVec( xyz[i], xyz[j], li.v[1], li.cL );
			if ( dSqr( li.v[1] ) < li.RnSqr ) {
				li.atom[i][li.num[i]] = j;
				li.atom[j][li.num[j]] = i;
				li.num[i] += 1;
				li.num[j] += 1;
			}
		}
	}
	return;
}

double wrap ( double pos, double cell ) {
	while ( pos < 0 ) { pos += cell; }
	while ( pos >= cell ) { pos -= cell; }
	return pos;
}

void minVec ( double* I, double* J, double* vec, double L ) {
	double dr, hf = L * 0.5;
	for ( int i = 0; i < 3; i++ ) {
		dr = J[i] - I[i];
		if ( dr < -hf ) { vec[i] = dr + L; }
		else { vec[i] = ( dr > hf ) ? dr - L : dr; }
	}
	return;
}

double dSqr ( double* vec ) {
	double sum = 0.0;
	for ( int i = 0; i < 3; i++ ) { sum += vec[i] * vec[i]; }
	return sum;
}
