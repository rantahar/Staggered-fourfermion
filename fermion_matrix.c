/****************************************************************
 * Fermion matrix and determinant
 *
 ****************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "Staggered.h"

extern bool **occupation_field;
extern int n_occupied[N_FLAVOR];
extern int **neighbour;
extern int **eta;
extern int **ksi;   //Staggered ksi matrix for link mass
extern int max_fluctuations;
extern double linkmass;
extern double m;
double *D;
double *Dinv;


/* Map indexes of even and odd sites to global indexes */
int * even_to_global_index;
int * odd_to_global_index;
int * global_to_even_index;
int * global_to_odd_index;


/* Remember the current determinant to avoid recalculating 
 * update after accepting a new configuration
 */
double current_determinant[N_FLAVOR];
double new_determinant[N_FLAVOR];
void update_current_determinant( int do_flavor[N_FLAVOR] ){
  for(int n=0; n<N_FLAVOR;n++) if(do_flavor[n]) current_determinant[n] = new_determinant[n];
}

void * wrapmalloc( int size, size_t bitsize ){
  void * pointer=NULL;
  pointer = calloc(size,bitsize);
  if( pointer == NULL ){
    printf("Could not allocate memory. Size %d,%lu.\n",size,bitsize);
  }
  return pointer;
}



/* Definitions of the fermion matrix */
#ifdef MASS_IN_MATRIX
#ifdef ANTIPERIODIC
void init_fermion_matrix( ) {
 static int init = 1;
 if(init==1){
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = wrapmalloc(n*n,sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    D[x1+n*x1] += m;
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int v1[ND];
      int bc_sign = 1;
      site_index_to_vector(x1,v1);
      if( v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      D[x1+n*x2] += 0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x1,v1);
      bc_sign = 1;
      if( v1[nu] == 0 ) bc_sign = -1;
    }
  }
  init = 0;
 }
}


void Dv( double * y, double * x ){
  /* Set up the fermion matrix */
  int n = VOLUME;
  int v1[ND];
  
  for( int x1=0;x1<n;x1++) y[x1] = m*x[x1];
  for( int x1=0;x1<n;x1++){
    site_index_to_vector(x1,v1);
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int bc_sign = 1;
      if( v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      y[x1] += bc_sign * ( 0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      bc_sign = 1;
      if( v1[nu] == 0 ) bc_sign = -1;
      y[x1] += bc_sign * ( -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;
    }
  }
}

#endif

#ifdef PERIODIC //Periodic boundaries in spatial directions
void init_fermion_matrix( )
{
 static int init = 1;
 if(init==1){
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = wrapmalloc(n*n,sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    D[x1+n*x1] += m;
    int nu=0;
    {
      //First the time direction, antiperiodic
      int x2 = neighbour[nu][x1];
      int v1[ND];
      int bc_sign = 1;
      site_index_to_vector(x1,v1);
      if( v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      D[x1+n*x2] += 0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x1,v1);
      bc_sign = 1;
      if( v1[nu] == 0 ) bc_sign = -1;
      D[x1+n*x2] += -0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;
    }
    //Spatial directions, periodic
    for( int nu=1;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      D[x1+n*x2] += 0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      D[x1+n*x2] += -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;
    }
  }
  init = 0;
 }
}

void Dv( double * y, double * x ){
  /* Set up the fermion matrix */
  int n = VOLUME;
  int v1[ND];
  
  for( int x1=0;x1<n;x1++) y[x1] = m*x[x1];
  for( int x1=0;x1<n;x1++){
    site_index_to_vector(x1,v1);
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int bc_sign = 1;
      if( nu==0 && v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      y[x1] += bc_sign * ( 0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      bc_sign = 1;
      if( nu==0 && v1[nu] == 0 ) bc_sign = -1;
      y[x1] += bc_sign * ( -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;
    }
  }
}

#endif

#else //MASS_IN_MATRIX

#ifdef ANTIPERIODIC //Antiperiodic boundaries
void init_fermion_matrix( )
{
 static int init = 1;
 if(init==1){
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = wrapmalloc(n*n,sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int v1[ND];
      int bc_sign = 1;
      site_index_to_vector(x1,v1);
      if( v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      D[x1+n*x2] += 0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x1,v1);
      bc_sign = 1;
      if( v1[nu] == 0 ) bc_sign = -1;
      D[x1+n*x2] += -0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;
    }
  }
  init = 0;
 }
}

void Dv( double * y, double * x ){
  /* Set up the fermion matrix */
  int n = VOLUME;
  int v1[ND];
  
  for( int x1=0;x1<n;x1++) y[x1] = 0;
  for( int x1=0;x1<n;x1++){
    site_index_to_vector(x1,v1);
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int bc_sign = 1;
      if( v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      y[x1] += bc_sign * ( 0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      bc_sign = 1;
      if( v1[nu] == 0 ) bc_sign = -1;
      y[x1] += bc_sign * ( -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;
    }
  }
}


#endif

#ifdef PERIODIC //Periodic boundaries in spatial directions
void init_fermion_matrix( )
{
 static int init = 1;
 if(init==1){
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = wrapmalloc(n*n,sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    int nu=0;
    {
      //First the time direction, antiperiodic
      int x2 = neighbour[nu][x1];
      int v1[ND];
      int bc_sign = 1;
      site_index_to_vector(x1,v1);
      if( v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      D[x1+n*x2] += 0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x1,v1);
      bc_sign = 1;
      if( v1[nu] == 0 ) bc_sign = -1;
      D[x1+n*x2] += -0.5 * bc_sign * eta[nu][x1] + bc_sign * linkmass * ksi[nu][x1] ;
    }
    //Spatial directions, periodic
    for( int nu=1;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      D[x1+n*x2] += 0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      D[x1+n*x2] += -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;
    }
  }
  init = 0;
 }
}

void Dv( double * y, double * x ){
  /* Set up the fermion matrix */
  int n = VOLUME;
  int v1[ND];
  
  for( int x1=0;x1<n;x1++) y[x1] = 0;
  for( int x1=0;x1<n;x1++){
    site_index_to_vector(x1,v1);
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int bc_sign = 1;
      if( nu==0 && v1[nu] == (Ldim[nu]-1) ) bc_sign = -1;
      y[x1] += bc_sign * ( 0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      bc_sign = 1;
      if( nu==0 && v1[nu] == 0 ) bc_sign = -1;
      y[x1] += bc_sign * ( -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1]) * x[x2] ;
    }
  }
}

#endif
#endif //MASS_IN_MATRIX


/* Set up indexes at even and odd sites */
void init_evenodd_index( )
{
  even_to_global_index = wrapmalloc( VOLUME/2,sizeof(int) );
  odd_to_global_index  = wrapmalloc( VOLUME/2,sizeof(int) );
  global_to_even_index = wrapmalloc( VOLUME,sizeof(int) );
  global_to_odd_index = wrapmalloc( VOLUME,sizeof(int) );
  int ie=0, io=0;
  int v[ND];
  for(int x=0; x<VOLUME; x++) {
    site_index_to_vector(x,v);
    int n = 0;
    for( int nu=0; nu<ND; nu++ ) n+=v[nu];
    if( n%2 == 0 ){
      even_to_global_index[ie] = x;
      global_to_even_index[x] = ie;
      global_to_odd_index[x] = -1;
      ie++;
    } else {
      odd_to_global_index[io] = x;
      global_to_odd_index[x] = io;
      global_to_even_index[x] = -1;
      io++;
    }
  }
}





/* The index of a site on a staggered block */
int block_index( int vector[ND] ){
  int index = vector[ND-1]%2;
  for (int dir=ND-1; dir--;) index = index*2 + vector[dir]%2;
  return index;
}


/* Return the inverse of the Dirac operator for two given sites */
double Dinvf( int io, int ie ){
    int boundary_sign = 1;
    int ve[ND],vo[ND];
    //printf(" %d %d\n",io,ie);

    //Get sites as vectors
    io = odd_to_global_index[io];
    site_index_to_vector(io,vo);
    ie = even_to_global_index[ie];
    site_index_to_vector(ie,ve);

    //printf(" (%d,%d,%d,%d) (%d,%d,%d,%d)\n",vo[0],vo[1],vo[2],vo[3],ve[0],ve[1],ve[2],ve[3]);

    //Translate so that the even site is on the first block
    for( int nu = 0; nu<ND; nu++ ) {
       int evenblock = ve[nu]/2;
       vo[nu] -= evenblock*2 ;
       ve[nu] -= evenblock*2 ;

       //Correct for antiperiodic boundaries
       if( vo[nu] < 0 ){
          boundary_sign *=-1;
          vo[nu] += Ldim[nu];
       }
    }

    //printf(" (%d,%d,%d,%d) (%d,%d,%d,%d)\n", vo[0],vo[1],vo[2],vo[3], ve[0],ve[1],ve[2],ve[3] );

    //Back to indexes
    io = site_vector_to_index(vo);
    io = global_to_odd_index[io];
    ie = block_index(ve);

    //printf(" %d %d %g\n\n",io,ie,boundary_sign*Dinv[ io + ie*VOLUME/2 ]);

    return boundary_sign*Dinv[ io + ie*VOLUME/2 ];
}




#ifdef FULL_DETERMINANT
/* Just calculate the determinant of the full fermion matrix */
double determinant( int f ){

  /* Initialize everything on the first time */
  static int init = 1;
  if(init) {
    init = 0;
    current_determinant[0] = current_determinant[1] = 1;
  }
  init_fermion_matrix( );
 
  /* The size of the matrix */
  int n=(VOLUME-n_occupied[f]);

  /* Determinant of a rank 0 matrix is 1 */
  if( n==0 ) {
    new_determinant[f] = 1;
    return( new_determinant[f]/current_determinant[f] );
  }

  /* Build a list of free sites */
  int v[ND];
  int *freelist;
  freelist = wrapmalloc( n,sizeof(int) );
  int ind=0;
  for(int x=0; x<VOLUME; x++) if( occupation_field[f][x]==0 ) {
    freelist[ind] = x;
    ind++;
  }
  
  /* Construct the matrix of free sites*/
  double *M = wrapmalloc(n*n,sizeof(double));
  for(int i=0; i<n*n; i++) M[i] = 0;
  for(int i1=0; i1<n; i1++) for(int i2=0; i2<n; i2++){
    M[i1+n*i2] = D[freelist[i1]+VOLUME*freelist[i2]] ;
  }

  /* Use lapack to calculate the LU decompisition */
  int info;
  int *ipiv;
  ipiv = wrapmalloc( n,sizeof(int) );
  LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );

  /* Calculate the determinant */
  double det = 1;
  for(int i=0; i<n; i++) {
    det *= M[i*n+i];
    if( ipiv[i] != i+1 ) det*=-1; 
  }
  
#ifndef MASS_IN_MATRIX
  /* Check for negative determinants (in most cases this is debugging) */
  /* Note that individual negative determinant are natural with MASS_IN_MATRIX,
   * both flavors need to be take into account to get a positive determinant*/
  if((det/current_determinant[f])<-1e-30) {
    static int n=0;
    printf("NEGATIVE determinant %d %g %g\n",n, det/current_determinant[f] ,det);
    n++;
  }
#endif

  free(M);
  free(ipiv);
  free(freelist);
  
  new_determinant[f] = det;
  return( new_determinant[f]/current_determinant[f] );
}
#endif





#ifdef FLUCTUATION_DETERMINANT
/* A more efficient method using a background configuration
 * and calculating the change in the determinant compared to
 * that background.
 * Uses an even to odd matrix and does not work with MASS_IN_MATRIX
 */

#ifdef MASS_IN_MATRIX
fluctuation determinant only written for an even-odd decomposed matrix
disable mass_in_matrix
#endif

/* The background configuration */
int n_bc_monomers[N_FLAVOR];
int **bc_occupation_field;
int ** bc_evenlist, ** bc_oddlist;

/* Matrices needed for the fluctuation matrix */
double ***BAP, ***PAC, ***BGC, **Ginv;
unsigned char **newsite,***newcombination;

/* Invert the matrix of propagators between occupied sites
 * Assign the current configuration as the background
 */
void update_background( int f )
{
  /* Initialize everything on the first call */
  static int init = 1;
  if(init==1){
    Ginv = wrapmalloc(N_FLAVOR,sizeof(double*));
    BAP = wrapmalloc(N_FLAVOR,sizeof(double*));
    PAC = wrapmalloc(N_FLAVOR,sizeof(double*));
    BGC = wrapmalloc(N_FLAVOR,sizeof(double*));
    newsite = wrapmalloc(N_FLAVOR,sizeof(unsigned char*));
    newcombination = wrapmalloc(N_FLAVOR,sizeof(int**));
    bc_occupation_field = wrapmalloc(N_FLAVOR,sizeof(int*));
    bc_evenlist = wrapmalloc(N_FLAVOR,sizeof(int*));
    bc_oddlist = wrapmalloc(N_FLAVOR,sizeof(int*));

    for( int n=0; n<N_FLAVOR; n++){
      Ginv[n] = NULL;
      BAP[n] = wrapmalloc(VOLUME/2,sizeof(double*));
      PAC[n] = wrapmalloc(VOLUME/2,sizeof(double*));
      for(int i=0; i<VOLUME/2; i++){
        BAP[n][i] = NULL;
        PAC[n][i] = NULL;
      }
   
      newsite[n] = wrapmalloc( VOLUME,sizeof(bool) );
#ifndef SAVE_MEMORY
      BGC[n] = wrapmalloc(VOLUME/2,sizeof(double*));
      for(int i=0; i<VOLUME/2; i++) BGC[n][i] = wrapmalloc(VOLUME/2,sizeof(double));
      newcombination[n] = wrapmalloc( VOLUME,sizeof(bool*) );
      for(int i=0; i<VOLUME; i++)
        newcombination[n][i] = wrapmalloc( VOLUME,sizeof(bool) );
#else
      BGC[n] = wrapmalloc(max_fluctuations,sizeof(double*));
      for(int i=0; i<max_fluctuations; i++) BGC[n][i] = wrapmalloc(max_fluctuations,sizeof(double));
#endif

      bc_occupation_field[n] = wrapmalloc( VOLUME,sizeof(int) );
      bc_evenlist[n] = wrapmalloc( VOLUME/2,sizeof(int) );
      bc_oddlist[n]  = wrapmalloc( VOLUME/2,sizeof(int) );
    }
    init_Dinv_cg();
  }

  /* The current configuration becomes the background */
  n_bc_monomers[f] = n_occupied[f];
  for (int x=0; x<VOLUME; x++) {
    bc_occupation_field[f][x] = occupation_field[f][x];
  }
  
  /* Size of the background matrix */
  int n_bc = n_occupied[f]/2;

  /* Invert the background matrix
   * Lapack fails for rank zero, so skip   */
  if( n_bc > 0 ){

    if( Ginv[f] != NULL) free(Ginv[f]);
    Ginv[f] = wrapmalloc(n_bc*n_bc,sizeof(double));

    int *ipiv = wrapmalloc( n_bc,sizeof(int) );
    int info;
  
    /* Find occupied and unoccupied sites, construct lists */
    int v[ND];
    int ioe=0, ioo=0, iue=0, iuo=0;
    for(int x=0; x<VOLUME; x++) {
      site_index_to_vector(x,v);
      int n = 0;
      for( int nu=0; nu<ND; nu++ ) n+=v[nu];
      if( n%2 == 0 ){
        if( occupation_field[f][x]==OCCUPIED ){
          bc_evenlist[f][ioe] = global_to_even_index[x];
          ioe++;
        }
      } else {
        if( occupation_field[f][x]==OCCUPIED ){
          bc_oddlist[f][ioo] = global_to_odd_index[x];
          ioo++;
        }
      }
    }
    if(ioe!=ioo || ioo!=n_bc){
      printf("Number on occupied sites doesn't match, i=%d, j=%d, n_bc=%d, n_occupied=%d\n",ioe,ioo,n_bc,n_occupied[f]);
      print_config();
      exit(1);
    }

    /* Construct the inverse of the matrix propagators over occupied sites */
    for(int ie=0; ie<n_bc; ie++) for(int io=0; io<n_bc; io++){
      Ginv[f][ie+n_bc*io] = Dinvf(bc_oddlist[f][io],bc_evenlist[f][ie]);
    }
    LAPACK_dgetrf( &n_bc, &n_bc, Ginv[f], &n_bc, ipiv, &info );
    if( info != 0 ) {
      printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
      exit(-1);
    }

    /* The inversion */
    int lwork=n_bc*n_bc;
    double *work;
    work = wrapmalloc( lwork,sizeof(double) );
    if (NULL == work) {
      printf("failed to allocate work matrix in update_background (size n_bc=%d) \n",n_bc);
      exit(-1);
    }
    LAPACK_dgetri(&n_bc, Ginv[f], &n_bc, ipiv, work, &lwork, &info);

    free(work);
    free(ipiv);
  }

  /* Update the saved determinant */
  current_determinant[f] = 1;

  /* Mark fluctuation matrix at all sites as not calculated */
  for( int a=0; a<VOLUME; a++ ) newsite[f][a] = 1;
#ifndef SAVE_MEMORY
  for( int a=0; a<VOLUME; a++ ) for( int b=0; b<VOLUME; b++ ) newcombination[f][a][b] = 1;
#endif
  for(int i=0; i<VOLUME/2; i++){
    if( PAC[f][i] != NULL ){
      free(PAC[f][i]);
      PAC[f][i] = NULL;
    }
    if( BAP[f][i] != NULL ){
      free(BAP[f][i]);
      BAP[f][i] = NULL; 
    }
  }


  if(init==1) init = 0;

  #ifdef DEBUG
  printf("New background for fermion type %d with %d monomers\n", f, n_bc_monomers[f]);
  #endif
}


/* Find site in background array using binary search */
int find_in_bc( int x, int * list, int n ){
  
  int i=n/2, u=n, l=0;
  do{
    if( list[i] == x ) return i;
    if( list[i] > x ) {
      u = i;
      i = l+(u-l)/2;
    } else {
      l = i;
      i = l+(u-l)/2 ;
    }
  } while(1);

  /*
  int i;
  for( i=0; i<n; i++) if( list[i] == x ) break ;
  return i;
  */
}


/* Calculate the determinant of the current configuration */
int n_fluctuations[N_FLAVOR];
int *added_evenlist,*added_oddlist,*removed_evenlist,*removed_oddlist;
double determinant( int f ){
  int n_bc = n_bc_monomers[f]/2;
  int n_added_even=0, n_removed_even=0, n_added_odd=0, n_removed_odd=0;
  static int ndet = 0;

  /* Construct a list of changes to the background */
  int v[ND];
  static int init = 1;
  if(init == 1){
    added_evenlist = wrapmalloc( max_fluctuations,sizeof(int) );
    added_oddlist  = wrapmalloc( max_fluctuations,sizeof(int) );
    removed_evenlist = wrapmalloc( max_fluctuations,sizeof(int) );
    removed_oddlist  = wrapmalloc( max_fluctuations,sizeof(int) );
    init = 0;
  }
  for(int x=0; x<VOLUME; x++) if(occupation_field[f][x] != bc_occupation_field[f][x]) {
    site_index_to_vector(x,v);
    int n = 0;
    for( int nu=0; nu<ND; nu++ ) n+=v[nu];
    if( n%2 == 0 ){
      if( occupation_field[f][x] == OCCUPIED ){
         added_evenlist[n_added_even] = global_to_even_index[x];
         n_added_even++;
      } else  {
         /* index in the background matrix required for removed sites */
         removed_evenlist[n_removed_even] = find_in_bc( global_to_even_index[x], bc_evenlist[f], n_bc_monomers[f]/2 );
         n_removed_even++;
      }
    } else {
      if( occupation_field[f][x]==OCCUPIED ){
         added_oddlist[n_added_odd] = global_to_odd_index[x];
         n_added_odd++;
      } else  {
         /* index in the background matrix required for removed sites */
         removed_oddlist[n_removed_odd] = find_in_bc( global_to_odd_index[x], bc_oddlist[f], n_bc_monomers[f]/2 );
         n_removed_odd++;
      }
    }
  }
  if( n_added_even-n_removed_even != n_added_odd-n_removed_odd ) {
    printf("Mismatch, %d %d  %d %d\n",n_added_even,n_added_odd,n_removed_even,n_removed_odd);
    exit(1);
  }

  /* The size of the fluctuation matrix */
  int n_fl = n_added_even + n_removed_odd;
  n_fluctuations[f] = n_fl;

  if(n_fl==0){  /* No fluctuations, trivial case */
    new_determinant[f] = 1;
#ifdef DEBUG
    printf("Det %d  %g %g (n_fl=0)\n", ndet, 1.0, new_determinant[f]/current_determinant[f] );
    ndet++;
#endif
    return new_determinant[f]/current_determinant[f];
  }


  /* Construct new parts of the fluctuation matrix */
  for( int b=0; b<n_added_even; b++) if( newsite[f][added_evenlist[b]] == 1 ){
    double D[n_bc];
    for( int l=0; l<n_bc; l++) D[l] = Dinvf(bc_oddlist[f][l],added_evenlist[b]); 
    BAP[f][added_evenlist[b]] = wrapmalloc(VOLUME/2,sizeof(double));
    for( int a=0; a<n_bc; a++)
    {
      double sum=0;
      double G[n_bc];
      for( int l=0; l<n_bc; l++) G[l] = Ginv[f][l+a*n_bc];
      for( int l=0; l<n_bc; l++) sum += D[l]*G[l];
      BAP[f][added_evenlist[b]][a] = sum;
    }
    newsite[f][added_evenlist[b]] = 0;
  }

  for( int a=0; a<n_added_odd; a++) if( newsite[f][VOLUME/2+added_oddlist[a]] == 1 ){
    double D[n_bc];
    for( int k=0; k<n_bc; k++) D[k] = Dinvf(added_oddlist[a],bc_evenlist[f][k]);
    PAC[f][added_oddlist[a]] = wrapmalloc(VOLUME/2,sizeof(double));
    for( int b=0; b<n_bc; b++)  { 
      double sum = 0;
      double G[n_bc];
      for( int k=0; k<n_bc; k++) G[k] = Ginv[f][b+k*n_bc];
      for( int k=0; k<n_bc; k++) sum += D[k]*G[k];
      PAC[f][added_oddlist[a]][b] = -sum;
    }
    newsite[f][VOLUME/2+added_oddlist[a]] = 0;
  }


#ifndef SAVE_MEMORY
  for( int b=0; b<n_added_even; b++) for( int a=0; a<n_added_odd; a++)
    if( newcombination[f][added_evenlist[b]][VOLUME/2+added_oddlist[a]] == 1 )
  {
    double sum = 0;
    double D[n_bc], P[n_bc];
    for( int l=0; l<n_bc; l++) D[l] = Dinvf(bc_oddlist[f][l],added_evenlist[b]);
    for( int l=0; l<n_bc; l++) P[l] = PAC[f][added_oddlist[a]][l];
    for( int l=0; l<n_bc; l++) sum += D[l]*P[l];
    BGC[f][added_evenlist[b]][added_oddlist[a]] =
      Dinvf(added_oddlist[a],added_evenlist[b]) + sum;
    newcombination[f][added_evenlist[b]][VOLUME/2+added_oddlist[a]] = 0;
  }
#else
  for( int b=0; b<n_added_even; b++) for( int a=0; a<n_added_odd; a++)
  {
    double sum = 0;
    double D[n_bc], P[n_bc];
    for( int l=0; l<n_bc; l++) D[l] = Dinvf(bc_oddlist[f][l],added_evenlist[b]);
    for( int l=0; l<n_bc; l++) P[l] = PAC[f][added_oddlist[a]][l];
    for( int l=0; l<n_bc; l++) sum += D[l]*P[l];
    BGC[f][b][a] =
      Dinvf(added_oddlist[a],added_evenlist[b]) + sum;
  }
#endif



  /* Pick entries for the fluctuation matrix */
  double *F;
  F = wrapmalloc(n_fl*n_fl,sizeof(double));
  for( int b=0; b<n_removed_odd; b++) for( int a=0; a<n_removed_even; a++)
    F[a*n_fl+b] = Ginv[f][removed_oddlist[b]+removed_evenlist[a]*n_bc];
  for( int a=0; a<n_added_odd; a++) for( int b=0; b<n_removed_odd; b++)
    F[(n_removed_even+a)*n_fl+b] = PAC[f][added_oddlist[a]][removed_oddlist[b]];
  for( int a=0; a<n_removed_even; a++) for( int b=0; b<n_added_even; b++)
    F[a*n_fl+(n_removed_odd+b)] = BAP[f][added_evenlist[b]][removed_evenlist[a]];
#ifndef SAVE_MEMORY
  for( int b=0; b<n_added_even; b++) for( int a=0; a<n_added_odd; a++)
    F[(n_removed_even+a)*n_fl+(n_removed_odd+b)] = BGC[f][added_evenlist[b]][added_oddlist[a]];
#else
  for( int b=0; b<n_added_even; b++) for( int a=0; a<n_added_odd; a++)
    F[(n_removed_even+a)*n_fl+(n_removed_odd+b)] = BGC[f][b][a];
#endif

  /* determinant from LU */
  double det=1;
  int ipiv[n_fl], info;
  LAPACK_dgetrf( &n_fl, &n_fl, F, &n_fl, ipiv, &info );
  for(int i=0; i<n_fl; i++) {
    det *= F[i*n_fl+i];
  }

  free(F);


  new_determinant[f] = det*det;

  return( new_determinant[f]/current_determinant[f] );
}


#endif //FLUCTUATION_DETERMINANT




/* Calculate weight for specific operations simply by changing the field
   and calling determinant().
   Unfortunately the background needs to be updated here when necessary. */
int init = 1;
double det_add_monomers(int x1, int x2, int do_flavor[N_FLAVOR]){

  #ifdef FLUCTUATION_DETERMINANT
  /* Initiate fluctuations */
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif
  
  /* Calculate the determinant after the change */
  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    /* Just change the field */
    occupation_field[f][x1] = OCCUPIED;
    occupation_field[f][x2] = OCCUPIED;
    n_occupied[f] += 2;
    /* And get the determinant */
    det *= determinant(f);
    occupation_field[f][x1] = UNOCCUPIED;
    occupation_field[f][x2] = UNOCCUPIED;
    n_occupied[f] -= 2;
  }

  return det;
}

double det_move_monomers(int x1, int x2, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  /* Initiate fluctuations */
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif

  /* Calculate the determinant after the change */
  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    /* Just change the field */
    occupation_field[f][x1] = UNOCCUPIED;
    occupation_field[f][x2] = OCCUPIED;
    det *= determinant( f );
    /* And get the determinant */
    occupation_field[f][x1] = OCCUPIED;
    occupation_field[f][x2] = UNOCCUPIED;
  }

  return det;
}


double det_remove_monomers(int x1, int x2, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  /* Initiate fluctuations */
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif

  /* Calculate the determinant after the change */  
  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    /* Just change the field */
    occupation_field[f][x1] = UNOCCUPIED;
    occupation_field[f][x2] = UNOCCUPIED;
    n_occupied[f] -= 2;
    det *= determinant(f);
    /* And get the determinant */
    occupation_field[f][x1] = OCCUPIED;
    occupation_field[f][x2] = OCCUPIED;
    n_occupied[f] += 2;
  }

  return det;
}

#ifdef MASS_IN_MATRIX

double det_add_monomers1(int x1, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  /* Initiate fluctuations */
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif
  
  /* Calculate the determinant after the change */  
  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]) {
    /* Just change the field */
    occupation_field[f][x1] = OCCUPIED;
    n_occupied[f] += 1;
    /* And get the determinant */
    det *= fabs(determinant(f));
    occupation_field[f][x1] = UNOCCUPIED;
    n_occupied[f] -= 1;
  }

  return det;
}

double det_remove_monomers1(int x1, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  /* Initiate fluctuations */
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif

  /* Calculate the determinant after the change */    
  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    /* Just change the field */
    occupation_field[f][x1] = UNOCCUPIED;
    n_occupied[f] -= 1;
    /* And get the determinant */
    det *= fabs(determinant(f));
    occupation_field[f][x1] = OCCUPIED;
    n_occupied[f] += 1;
  }

  return det;
}
#endif




#ifdef MEASUREVEV
/* Calculate the link vev */
void measure_vev( double * linkvev, double * sitevev ){
  init_fermion_matrix( );

  int n = VOLUME ;
  
  int info;
  int *ipiv;
  ipiv = wrapmalloc( n,sizeof(int) );
  double *Dinv = wrapmalloc( n*n,sizeof(double) );
  for(int i=0; i<n*n; i++) Dinv[i] = 0;

  /* Construct the even to odd matrix
   * using only unoccupied sites
   */
  for(int i1=0; i1<n; i1++){
    if( occupation_field[0][i1] == UNOCCUPIED ){
      for(int i2=0; i2<n; i2++){
        if( occupation_field[0][i2] == UNOCCUPIED ) {
          Dinv[i1+n*i2] = D[ i1+VOLUME*i2 ] ;
        }
      }
    } else{
      Dinv[i1+n*i1] = 1;
    }
  }

  /* Invert using lapack routines */
  int lwork=n*n;
  double *work;
  work = wrapmalloc( lwork,sizeof(double) );
  LAPACK_dgetrf( &n, &n, Dinv, &n, ipiv, &info );
  if( info != 0 ) {
    printf("init_Dinv: sgetrf returned an error %d! \n", info);
    exit(-1);
  }

  LAPACK_dgetri(&n, Dinv, &n, ipiv, work, &lwork, &info); 
  if( info != 0 ) {
    printf("init_Dinv: sgetri returned an error %d! \n", info);
    exit(-1);
  }
  
  for(int i1=0; i1<n; i1++){
    int mu=0;
    int i2 = neighbour[mu][i1];
    int v1[ND];
    int bc_sign = 1;
    site_index_to_vector(i1,v1);
    if( v1[mu] == (Ldim[mu]-1) ) bc_sign = -1;
    *linkvev += bc_sign*ksi[mu][i1] * Dinv[i1+n*i2];

    
    int omu = ( mu + ND ) % NDIRS;
    i2 = neighbour[omu][i1];
    bc_sign = 1;
    if( v1[mu] == 0 ) bc_sign = -1;
    *linkvev += bc_sign*ksi[mu][i1] * Dinv[i1+n*i2];
    
    for( int mu=1;mu<ND;mu++){
      int i2 = neighbour[mu][i1];
      *linkvev += ksi[mu][i1] * Dinv[i1+n*i2];

      int omu = ( mu + ND ) % NDIRS;
      i2 = neighbour[omu][i1];
      *linkvev += ksi[mu][i1] * Dinv[i1+n*i2];
    }

    if( occupation_field[0][i1] == UNOCCUPIED )
      *sitevev += Dinv[i1+n*i1];
    //if(i1 == 10) printf(" site %g\n", Dinv[i1+n*i1]); 
    
  }

  /*linkvev = 0;
  for(int i1=0; i1<n; i1++){
    for(int mu=0; mu<ND; mu++) {
      int i2 = neighbour[mu][i1];
      linkvev += ksi[mu][i2]*Dinv[i1+n*i2];
    }
  }*/
  
  
  free(work);
  free(ipiv);
  free(Dinv);

}
#endif



















