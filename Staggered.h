#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <lapacke/lapacke.h>
#include <stdbool.h>
//#include <lapacke.h>
#include <time.h>
#include <sys/time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

/* Lattice size and dimensions */
#define ND 4
/* Size of the lattice dimensions */
static int Ldim[ND] = { 12, 12, 12, 12 };

#define NDIRS (2*ND)
int VOLUME;

// Choose update method
#define WORM_UPDATE
//#define LOCAL_UPDATE

//Choose spatial boundary conditions (temporal alway antiperiodic)
//#define PERIODIC
#define ANTIPERIODIC

//Uncomment to include site mass in the determinant
//(does not work with fluctuation determinant)
//#define MASS_IN_MATRIX

/* Flavors */
#define N_FLAVOR 2

/* Choose here the method for calculating the derivative */
//#define FULL_DETERMINANT
#define FLUCTUATION_DETERMINANT

/* Toggle direct measurement for site and link vev */
//#ifdef MEASUREVEV

/* Save memory to fit large lattices
 * Reduces memory use significantly and
 * slightly increases compute times
 */
#define SAVE_MEMORY

/* field values */
#define UNOCCUPIED 0
#define OCCUPIED 1
#define NONE -100   //A meta value for sites that don't exist


/* In Staggered.c */
int site_vector_to_index( int vector[] );
void site_index_to_vector(int index, int vector[] );
void print_config();

/* In fermion_matrix.c */
double fM_index( int x1[ND], int x2[ND] );
void calc_Dinv( );
double determinant();


double det_add_monomers(int x1, int x2, int do_flavor[N_FLAVOR]);
double det_move_monomers(int x1, int x2, int do_flavor[N_FLAVOR]);
double det_remove_monomers(int x1, int x2, int do_flavor[N_FLAVOR]);
#ifdef MASS_IN_MATRIX
double det_add_monomers1(int x1, int do_flavor[N_FLAVOR]);
double det_remove_monomers1(int x1, int do_flavor[N_FLAVOR]);
#endif
void update_current_determinant( int do_flavor[N_FLAVOR] );
void update_background( int monomer );

void measure_vev( double * linkvev, double * sitevev );

/* update functions */
double worm_update( int *additions, int *removals, int *m_additions, int *m_removals, int *switches );
void local_update( int *additions, int *removals, int *moves,int *m_additions, int *m_removals, int *m_moves, int *switches );



