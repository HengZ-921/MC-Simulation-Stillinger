#include "global.h"

void updateRDF ( double bin, double** xyz, int* hist, NBL li ) {
	int N;
	for ( int i = 0; i < li.NAtom; i++ ) {
		for ( int j = i + 1; j < li.NAtom; j++ ) {
			minVec( xyz[i], xyz[j], li.v[1], li.cL );
			N = sqrt( dSqr( li.v[1] ) ) / bin;
			hist[N] += 1;
		}
	}
	return;
}

double cubeDiff ( int n ) {
	double diff = 3.0 * n * n + 3.0 * n + 1.0;
	return 1.0 / diff;
}

void outputRDF ( int Nacc, int N, double bin, int* hist, NBL li ) {
	double c = pow( bin, -3.0) * 2.0 / Nacc / li.NAtom / sPi / li.rho;
	double xval, yval;
	for (int i = 0; i < N; i++ ) {
		xval = bin * ( i + 0.5 );
		yval = hist[i] * c * cubeDiff( i );
		printf( "%lf\t%lf\n", xval, yval );
	}
	return;
}

int* accumRDF ( int stp, int freq, double bin, int* Nacc, double dr, double T, NBL li, double** xyz ) {
	int count = 0,  Nbin = li.cL / bin;
	int* hist = calloc( Nbin, sizeof( int ) );
	double dtt = 0.0;
	while ( count <= stp ) {
		if ( count == *Nacc * freq ) {
			updateRDF( bin, xyz, hist, li );
			*Nacc += 1;
		//	printf("%d\t%d\n", count, *Nacc);
		}
		count += trialMove( li.v[0], dr, T, li, xyz, &dtt );
	//	printf("%d\t %lf\n", count, dtt);
		refreshNBL( &dtt, xyz, li );
	}
	return hist;
}

void prodSweep ( double dr, double T, double** xyz, NBL li ) {
	int *hist, stp, freq, Nacc = 0;
	double bin, lim;
	scanf( "%*s %d\n", &stp );
	scanf( "%*s %lf\n", &bin );
	scanf( "%*s %lf\n", &lim );
	scanf( "%*s %d\n", &freq );
	updateNBL( xyz, li );
	hist = accumRDF( stp, freq, bin, &Nacc, dr, T, li, xyz );
	int Nbin = lim / bin;
	//printf("%d\t%d\n", Nacc, Nbin);
	outputRDF( Nacc, Nbin, bin, hist, li );
	free( hist );
	return;
}
