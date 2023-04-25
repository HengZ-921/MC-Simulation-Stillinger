#include "global.h"

double ranf () { return rand() / (double) RAND_MAX * 2 - 1.0; }

double rdlg () { return -log( rand() / (double) RAND_MAX ); }

int trialMove ( double* trial, double dr, double T, NBL li, double** xyz, double* dtt ) {
	int accept = 0, K = rand() % li.NAtom;
	double dV;
	for ( int i = 0; i < 3; i++ ) { trial[i] = wrap( xyz[K][i] + dr * ranf(), li.cL ); }
	dV = potSingle( K, trial, li, xyz ) - potSingle( K, xyz[K], li, xyz );
	if ( dV <= T * rdlg() ) {
		accept = 1;
		for ( int i = 0; i < 3; i++ ) { xyz[K][i] = trial[i]; }
	}
	minVec( xyz[K], li.coord[K], trial, li.cL );
	*dtt += sqrt( dSqr( trial ) );
	return accept;
}

void refreshNBL ( double* dtt, double** xyz, NBL li ) {
	if ( *dtt > li.Rk * 0.5 ) {
		updateNBL( xyz, li );
		*dtt = 0.0;
	}
	return;	
}

void sweepSteps ( int stp, double dr, double T, NBL li, double** xyz ) {
	int count = 0;
	double dtt = 0.0;
	while ( count < stp ) {
		count += trialMove( li.v[0], dr, T, li, xyz, &dtt );
	//	printf("%d\t %lf\n", count, dtt);
		refreshNBL( &dtt, xyz, li );
	}
}

void eqSweep ( double *dr, double* T, double** xyz, NBL li ) {
	int stp;
	scanf( "%*s %lf\n", dr );
	scanf( "%*s %lf\n", T );
	scanf( "%*s %d\n", &stp );
	srand( time(0) );
	updateNBL( xyz, li );
	sweepSteps( stp, *dr, *T, li, xyz );
	return;
}
