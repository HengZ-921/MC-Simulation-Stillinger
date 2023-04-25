#include "global.h"

double div_a ( double r ) { return 1.0 / ( r - alpha ) ; }

double f_2BD ( double r ) {
	return A * ( B * pow( r, -4.0 ) - 1.0 ) * exp( div_a( r ) );
}

double h_exp ( double r1, double r2 ) {
	return lambda * exp( gama * ( div_a( r1 ) + div_a ( r2 ) ) );
}

double cosTheta ( double* v1, double* v2, double r1, double r2 ) {
	double sum = 0.0;
	for ( int i = 0; i < 3; i++ ) { sum += v1[i] * v2[i]; }
	sum = sum / r1 / r2 + 1.0 / 3.0;
	return sum * sum;
}

double h_3BD ( double* v1, double* v2, double r1, double a2 ) {
	double h = 0.0, r2 = dSqr( v2 );
	if ( r2 < a2 ) {
		r2 = sqrt( r2 );
		h = h_exp( r1, r2 ) * cosTheta( v1, v2, r1, r2 );
	}
	return h;
}

double potSingle ( int I, double* Ixyz, NBL li, double** xyz ) {
	double r1, E2 = 0.0, E3 = 0.0, a2 = alpha * alpha;
	int J, K;
	for ( int j = 0; j < li.num[I]; j++ ) {
		J = li.atom[I][j];
		minVec( Ixyz, xyz[J], li.v[1], li.cL );
		r1 = dSqr( li.v[1] );
		if ( r1 < a2 ) {
			r1 = sqrt( r1 );
			E2 += f_2BD( r1 );
			for ( int k = j + 1; k < li.num[I]; k++ ) {
				K = li.atom[I][k];
				minVec( Ixyz, xyz[K], li.v[2], li.cL );
				E3 += h_3BD( li.v[1], li.v[2], r1, a2 );
			}
			minVec( xyz[J], Ixyz, li.v[1], li.cL );
			for ( int k = 0; k < li.num[J]; k++ ) {
				K = li.atom[J][k];
				minVec( xyz[J], xyz[K], li.v[2], li.cL );
				E3 += h_3BD( li.v[1], li.v[2], r1, a2 );
			}
			E3 -= h_exp( r1, r1 ) * cosTheta( li.v[1], li.v[1], r1, r1 );
		}
	}
//	printf( "%lf\t%lf\t", E2, E3 );
	return E2 + E3;
}
