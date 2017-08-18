#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <lapacke/lapacke.h>
#include <time.h>
#include <sys/time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

/* Lattice size and dimensions */
#define ND 4
/* Size of the lattice dimensions */
static int Ldim[ND] = { 2, 2, 2, 2 };

#define NDIRS (2*ND)
int VOLUME;

int * block_lattice_sites; //does not include positive boundary
int ** move_pairs;
int n_move_pairs;
int this_block;
#define BLOCKSIZE 4
#define N_BLOCK_UPDATES 1
int BLOCK_VOLUME;

#define PERIODIC
//#define ANTIPERIODIC
//#define MASS_IN_MATRIX

/* Flavors */
#define N_FLAVOR 2
#define OCCUPIED 1
#define NONE -100   //A meta value for sites that don't exist

/* Choose here the method for calculatign the derivative */
//#define SMALL_U
//#define LARGE_U
//#define FULL_DETERMINANT
#define FLUCTUATION_DETERMINANT


/* In Staggered.c */
int site_vector_to_index( int vector[] );
void site_index_to_vector(int index, int vector[] );
void print_config();

/* In fermion_matrix.c */
double fM_index( int x1[ND], int x2[ND] );
void calc_Dinv( );
double determinant();

#ifndef ANTIPERIODIC_MASS
double det_add_monomers(int x1, int x2, int do_flavor[N_FLAVOR]);
double det_move_monomers(int x1, int x2, int do_flavor[N_FLAVOR]);
double det_remove_monomers(int x1, int x2, int do_flavor[N_FLAVOR]);
#else
double det_add_monomers(int x1, int do_flavor[N_FLAVOR]);
double det_remove_monomers(int x1, int do_flavor[N_FLAVOR]);
#endif
void update_current_determinant( int do_flavor[N_FLAVOR] );
void update_background( int monomer );







