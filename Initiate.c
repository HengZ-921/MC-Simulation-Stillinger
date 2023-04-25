#include "global.h"

double** readxyz ( int* numAtom, double* cell ) {
	/* read the coordinates from input file and wrap them in the box. */
	char fileName[100];
	scanf( "%s\n", fileName );
	scanf( "%*s %lf\n", cell );
	FILE* fp = fopen( fileName, "r" );
	if ( fp == NULL ) { 
		printf("No file!\n");
		return 0;
	}
	fscanf( fp, "%d\n", numAtom);
	double x, y, z;
	double **Coord = malloc( *numAtom * sizeof(double*) );
	fscanf( fp, "%*s\n" );
	for ( int i = 0; i < *numAtom; i++ ) {
		Coord[i] = malloc( 3 * sizeof(double) );
		fscanf( fp, "%*s %lf %lf %lf\n", &x, &y, &z );
		Coord[i][0] = wrap( x / sigma, *cell );
		Coord[i][1] = wrap( y / sigma, *cell );
		Coord[i][2] = wrap( z / sigma, *cell );
	}
	fclose( fp );
	return Coord;
}

void free2D ( int N, double** mtx ) {
	for ( int i = 0; i < N; i++ ) { free(mtx[i]); }
	free(mtx);
	return;
}

void cleanAll ( int N, NBL nbLi, double** xyz ) {
	for ( int i = 0; i < N; i++ ) { free( nbLi.atom[i] ); }
	free( nbLi.atom );
	free2D( N, nbLi.coord );
	free2D( 3, nbLi.v );
	free2D( N, xyz );
	free(nbLi.num);
	return;
}
