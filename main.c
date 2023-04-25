#include "global.h"

int main () {
	int numAtom;
	double cell, drMx, temp;
	double** Coord = readxyz( &numAtom, &cell );
	NBL nbLi = initNBL( numAtom, cell );
//	updateNBL( Coord, nbLi );
//	for ( int i = 0; i < numAtom; i++ ) { printf( "E = %lf\n", potSingle( i, Coord[i], nbLi, Coord ) ); }
	eqSweep( &drMx, &temp, Coord, nbLi );
	prodSweep( drMx, temp, Coord, nbLi );

	cleanAll( numAtom, nbLi, Coord );
	return 0;
}
//	for ( int i = 0; i < numAtom; i++ ) { printf("%lf\t%lf\t%lf\n", Coord[i][0], Coord[i][1], Coord[i][2]); }
//	for ( int i = 0; i < 2; i++ ) { printf( "%d\n", nbLi.num[i]); } //for ( int j = 0; j < nbLi.num[i]; j++ ) { printf("%d ", nbLi.atom[i][j]); } 
