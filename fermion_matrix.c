/****************************************************************
 * Simulate 2 flavors of staggered fermions using the fermion bag  
 * algorithm (arXiv:0910.5736).
 *
 ****************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "Staggered.h"

extern int **occupation_field;
extern int n_occupied[N_FLAVOR];
extern int **neighbour;
extern int **eta;
extern int **ksi;   //Staggered ksi matrix for link mass
extern int max_fluctuations;
extern double linkmass;
extern double m;
double *D;
double * Dinv;


/* Map indexes of even and odd sites to global indexes */
int * even_to_global_index;
int * odd_to_global_index;
int * global_to_even_index;
int * global_to_odd_index;


//Remember the current determinant to avoid recalculating
double current_determinant[N_FLAVOR];
double new_determinant[N_FLAVOR];
/* Update the current determinant after accepting a configuration */
void update_current_determinant( int do_flavor[N_FLAVOR] ){
  for(int n=0; n<N_FLAVOR;n++) if(do_flavor[n]) current_determinant[n] = new_determinant[n];
}


#ifdef MASS_IN_MATRIX //Antiperiodic boundaries
#ifdef ANTIPERIODIC
void init_fermion_matrix( )
{
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = malloc(n*n*sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    D[x1+n*x1] += m;
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int v1[ND],v2[ND];
      site_index_to_vector(x1,v1);
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu][x2] + linkmass * ksi[nu][x2] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu][x2] - linkmass * ksi[nu][x2] ;
      }

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu][x2] - linkmass * ksi[nu][x2] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu][x2] + linkmass * ksi[nu][x2] ;
      }
    }
  }
}
#endif

#ifdef PERIODIC //Periodic boundaries in spatial directions
void init_fermion_matrix( )
{
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = malloc(n*n*sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    D[x1+n*x1] += m;
    int nu=0;
    {
      int x2 = neighbour[nu][x1];
      int v1[ND],v2[ND];
      site_index_to_vector(x1,v1);
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu%ND][x1] + linkmass * ksi[nu][x1] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu%ND][x1] - linkmass * ksi[nu][x1] ;
      }

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu][x1] - linkmass * ksi[nu][x1] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;
      }
    }
    for( int nu=1;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int v1[ND],v2[ND];
      site_index_to_vector(x1,v1);
      site_index_to_vector(x2,v2);
      D[x1+n*x2] += 0.5 * eta[nu%ND][x1] + linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x2,v2);
      D[x1+n*x2] += -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;
    }
  }
}
#endif

#else //MASS_IN_MATRIX

#ifdef ANTIPERIODIC //Antiperiodic boundaries
void init_fermion_matrix( )
{
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = malloc(n*n*sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    for( int nu=0;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int v1[ND],v2[ND];
      site_index_to_vector(x1,v1);
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu][x2] + linkmass * ksi[nu][x2] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu][x2] - linkmass * ksi[nu][x2] ;
      }

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu][x2] - linkmass * ksi[nu][x2] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu][x2] + linkmass * ksi[nu][x2] ;
      }
      
    }
  }

  /*for( int i=0;i<n;i++) {
    for( int j=0;j<n;j++) {
      printf(" %4.2f ",D[i+n*j]);
    }
    printf("\n");
  }*/
}
#endif

#ifdef PERIODIC //Periodic boundaries in spatial directions
void init_fermion_matrix( )
{
  /* Set up the fermion matrix */
  int n = VOLUME;
  D = malloc(n*n*sizeof(double));
  for( int i=0;i<n*n;i++)  D[i] = 0;
  
  for( int x1=0;x1<n;x1++){
    int nu=0;
    {
      int x2 = neighbour[nu][x1];
      int v1[ND],v2[ND];
      site_index_to_vector(x1,v1);
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu%ND][x1] + linkmass * ksi[nu][x1] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu%ND][x1] - linkmass * ksi[nu][x1] ;
      }

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x2,v2);
      if( v2[nu] > v1[nu] ) {
        D[x1+n*x2] += 0.5 * eta[nu][x1] - linkmass * ksi[nu][x1] ;
      } else {
        D[x1+n*x2] += -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;
      }
    }
    for( int nu=1;nu<ND;nu++){
      int x2 = neighbour[nu][x1];
      int v1[ND],v2[ND];
      site_index_to_vector(x1,v1);
      site_index_to_vector(x2,v2);
      D[x1+n*x2] += 0.5 * eta[nu%ND][x1] + linkmass * ksi[nu][x1] ;

      int onu = ( nu + ND ) % NDIRS;
      x2 = neighbour[onu][x1];
      site_index_to_vector(x2,v2);
      D[x1+n*x2] += -0.5 * eta[nu][x1] + linkmass * ksi[nu][x1] ;
    }
  }

}
#endif
#endif //MASS_IN_MATRIX


void init_evenodd_index( )
{
  /* Set up indexes at even and odd sites */
  even_to_global_index = malloc( VOLUME/2*sizeof(int) );
  odd_to_global_index  = malloc( VOLUME/2*sizeof(int) );
  global_to_even_index = malloc( VOLUME*sizeof(int) );
  global_to_odd_index = malloc( VOLUME*sizeof(int) );
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





void write_Dinv(){
  FILE * Dinv_file;
  char filename[100] = "Dinv";
  int n = VOLUME/2;

  Dinv_file = fopen(filename,"wb");
  if (Dinv_file){
    fwrite(Dinv, n*n, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {
    printf("Could not write Dinv file\n");
    exit(1);
  }
}


void init_Dinv(){
  
  init_evenodd_index();

  FILE * Dinv_file;
  char filename[100] = "Dinv";
  int n = VOLUME/2;
  Dinv = malloc(n*n*sizeof(double));

  Dinv_file = fopen(filename,"rb");
  if (Dinv_file){
    printf("Reading Dinv\n");
    fread(Dinv, n*n, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {

    printf("No Dinv file\n");
    init_fermion_matrix( );
    
    int info;
    int *ipiv;
    ipiv = malloc( n*sizeof(int) );
    for(int i=0; i<n*n; i++) Dinv[i] = 0;

    /* Construct the even to odd matrix */
    for(int ie=0; ie<n; ie++) for(int io=0; io<n; io++){
      Dinv[ie+n*io] = D[even_to_global_index[ie]+VOLUME*odd_to_global_index[io]] ;
    }

    int lwork=n*n;
    double *work;
    work = malloc( lwork*sizeof(double) );
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
  
    free(work);
    free(ipiv);
    free(D);
    write_Dinv();
  }
  
}




#ifndef FLUCTUATION_DETERMINANT


#ifdef FULL_DETERMINANT
/* Just calculate the determinant of the full fermion matrix */
double determinant( int f ){
  static int init = 1;
  if(init) {
    init_fermion_matrix( );
    init = 0;
    current_determinant[0] = current_determinant[1] = 1;
  }
  int n=(VOLUME-n_occupied[f]);
  if( n==0 ) {
    new_determinant[f] = 1;
    return( new_determinant[f]/current_determinant[f] );
  }

  //printf("Full determinant, n=%d\n",n);

  int info;
  int *ipiv;
  ipiv = malloc( n*sizeof(int) );
  double det = 1;
  double *M = malloc(n*n*sizeof(double));
  for(int i=0; i<n*n; i++) M[i] = 0;

  /* Construct the even to odd matrix of unoccupied sites */
  int v[ND];
  int *freelist;
  freelist = malloc( n*sizeof(int) );
  int ind=0;
  for(int x=0; x<VOLUME; x++) if( occupation_field[f][x]==0 ) {
    freelist[ind] = x;
    ind++;
  }

  for(int i1=0; i1<n; i1++) for(int i2=0; i2<n; i2++){
    M[i1+n*i2] = D[freelist[i1]+VOLUME*freelist[i2]] ;
  }

  /*for(int i1=0; i1<n; i1++) {
    for(int i2=0; i2<n; i2++) {
      printf(" %4.2f", M[i1+n*i2]);
    }
    printf("\n");
  }
  printf("\n");
  */


  LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
  //if(info) printf("Zero at %d\n",info); 
  for(int i=0; i<n; i++) {
    det *= M[i*n+i];
    if( ipiv[i] != i+1 ) det*=-1; 
    //printf(" %g %g %d\n",det, M[i*n+i], ipiv[i]!=i+1 );
  }
    
  if((det/current_determinant[f])<-1e-30) {
    static int n=0;
    printf("NEGATIVE determinant %d %g %g\n",n, det/current_determinant[f] ,det);
    n++;
  }

  free(M);
  free(ipiv);
  free(freelist);
  new_determinant[f] = det;
  //printf("det %g\n",det);
  return( new_determinant[f]/current_determinant[f] );
}
#endif





#else // defined FLUCTUATION_DETERMINANT
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
int **newsite, ***new_combination;

/* Invert the matrix of propagators between occupied sites
 * Assign the current configuration as the background
 */
void update_background( int f )
{
  static int init = 1;
  if(init==1){
    Ginv = malloc(N_FLAVOR*sizeof(double*));
    BAP = malloc(N_FLAVOR*sizeof(double*));
    PAC = malloc(N_FLAVOR*sizeof(double*));
    BGC = malloc(N_FLAVOR*sizeof(double*));
    //new_combination = malloc(N_FLAVOR*sizeof(int*));
    newsite = malloc(N_FLAVOR*sizeof(int*));
    bc_occupation_field = malloc(N_FLAVOR*sizeof(int*));
    bc_evenlist = malloc(N_FLAVOR*sizeof(int*));
    bc_oddlist = malloc(N_FLAVOR*sizeof(int*));

    for( int n=0; n<N_FLAVOR; n++){
      Ginv[n] = malloc(VOLUME*VOLUME/4*sizeof(double));
      BAP[n] = malloc(VOLUME/2*sizeof(double*));
      PAC[n] = malloc(VOLUME/2*sizeof(double*));
      BGC[n] = malloc(VOLUME/2*sizeof(double*));
      //new_combination[n] = malloc( sizeof(int*)*VOLUME/2 );
      for(int i=0; i<VOLUME/2; i++){
        BAP[n][i] = malloc(VOLUME/2*sizeof(double));
        PAC[n][i] = malloc(VOLUME/2*sizeof(double));
        BGC[n][i] = malloc(VOLUME/2*sizeof(double));
        //new_combination[n][i] = malloc( sizeof(int)*VOLUME/2 );
     }
   
      newsite[n] = malloc( sizeof(int)*VOLUME );
  
      bc_occupation_field[n] = malloc( VOLUME*sizeof(int) );
      bc_evenlist[n] = malloc( VOLUME/2*sizeof(int) );
      bc_oddlist[n]  = malloc( VOLUME/2*sizeof(int) );
    }
    init_Dinv( );
  }

  /* The current configuration becomes the background */
  n_bc_monomers[f] = n_occupied[f];
  for (int x=0; x<VOLUME; x++) {
    bc_occupation_field[f][x] = occupation_field[f][x];
  }
  
  int n_bc = n_occupied[f]/2;

  /*If no occupied sites, determinant of rank 0 is 1 */
  if( n_bc > 0 ){
    int *ipiv = malloc( n_bc*sizeof(int) );
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
          //printf("Occupied site %d (odd %d), bc index %d\n",x,bc_oddlist[ioo],ioo);
          ioo++;
        }
      }
    }
    if(ioe!=ioo || ioo!=n_bc){
      printf("Number on occupied sites doesn't match, i=%d, j=%d, n_bc=%d, n_occupied=%d\n",ioe,ioo,n_bc,n_occupied[f]);
      print_config();
      exit(1);
    }

    /* Construct the inverse of the occupied to occupied propagator */
    for(int ie=0; ie<n_bc; ie++) for(int io=0; io<n_bc; io++){
      Ginv[f][ie+n_bc*io] = Dinv[bc_oddlist[f][io]+VOLUME/2*bc_evenlist[f][ie]];
    }
    LAPACK_dgetrf( &n_bc, &n_bc, Ginv[f], &n_bc, ipiv, &info );
    if( info != 0 ) {
      printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
      exit(-1);
    }
    /* The inversion */
    int lwork=n_bc*n_bc;
    double *work;
    work = malloc( lwork*sizeof(double) );
    if (NULL == work) {
      printf("failed to allocate work matrix in update_background (size n_bc=%d) \n",n_bc);
      exit(-1);
    }
    LAPACK_dgetri(&n_bc, Ginv[f], &n_bc, ipiv, work, &lwork, &info);

    free(work);
    free(ipiv);
  }


  current_determinant[f] = 1;
  /* Mark fluctuation matrix at all sites as not calculated */
  for( int a=0; a<VOLUME; a++ ) newsite[f][a] = 1;
  //for( int a=0; a<VOLUME/2; a++ ) for( int b=0; b<VOLUME/2; b++ ) 
  //  new_combination[f][a][b] = 1;

  if(init==1) init = 0;
  #ifdef DEBUG
  printf("New background for fermion type %d with %d monomers\n", f, n_bc_monomers[f]);
  #endif
}


int n_fluctuations[N_FLAVOR];
/* Determinants after adding or removing links */
double determinant( int f ){
  int n_bc = n_bc_monomers[f]/2;
  int n_added_even=0, n_removed_even=0, n_added_odd=0, n_removed_odd=0;
  static int ndet = 0;

  /* Construct a list of changes to the background */
  int v[ND];
  int *added_evenlist,*added_oddlist,*removed_evenlist,*removed_oddlist;
  added_evenlist = malloc( max_fluctuations*sizeof(int) );
  added_oddlist  = malloc( max_fluctuations*sizeof(int) );
  removed_evenlist = malloc( max_fluctuations*sizeof(int) );
  removed_oddlist  = malloc( max_fluctuations*sizeof(int) );
  for(int x=0; x<VOLUME; x++) if(occupation_field[f][x] != bc_occupation_field[f][x]) {
    site_index_to_vector(x,v);
    int n = 0;
    for( int nu=0; nu<ND; nu++ ) n+=v[nu];
    if( n%2 == 0 ){
      if( occupation_field[f][x] == OCCUPIED ){
         added_evenlist[n_added_even] = global_to_even_index[x];
         //printf("Even site number %d, %d, added\n",x, added_evenlist[n_added_even]);
         n_added_even++;
      } else  {
         /* index in the background matrix required for removed sites */
         int i;
         for( i=0; i<n_bc_monomers[f]; i++) if( bc_evenlist[f][i] == global_to_even_index[x] ) break ;
         removed_evenlist[n_removed_even] = i;
         //printf("Even site number %d, %d, removed\n",i, bc_evenlist[i]);
         n_removed_even++;
      }
    } else {
      if( occupation_field[f][x]==OCCUPIED ){
         added_oddlist[n_added_odd] = global_to_odd_index[x];
         //printf("Odd site number %d, %d, added\n",x, added_oddlist[n_added_odd]);
         n_added_odd++;
      } else  {
         int i;
         for( i=0; i<n_bc_monomers[f]; i++) if( bc_oddlist[f][i] == global_to_odd_index[x] ) break ;
         removed_oddlist[n_removed_odd] = i;
         //printf("Odd site number %d, %d, removed\n",i, bc_oddlist[i]);
         n_removed_odd++;
      }
    }
  }
  if( n_added_even-n_removed_even != n_added_odd-n_removed_odd ) {
    printf("Mismatch, %d %d  %d %d\n",n_added_even,n_added_odd,n_removed_even,n_removed_odd);
    exit(1);
  }

  int n_fl = n_added_even + n_removed_odd;
  n_fluctuations[f] = n_fl;

  //printf("Calculating fluctuation with %d additions and removals\n",n_fl);
  if(n_fl==0){  /* No fluctuations, trivial case */
    new_determinant[f] = 1;
#ifdef DEBUG
    printf("Det %d  %g %g (n_fl=0)\n", ndet, 1.0, new_determinant[f]/current_determinant[f] );
    ndet++;
#endif
    free(removed_evenlist);
    free(removed_oddlist);
    free(added_evenlist);
    free(added_oddlist);
    return new_determinant[f]/current_determinant[f];
  }

  /* Construct new parts of the fluctuation matrix */
  for( int b=0; b<n_added_even; b++) if( newsite[f][added_evenlist[b]] == 1 ){
    for( int a=0; a<n_bc; a++)
    {
      double sum=0;
      for( int l=0; l<n_bc; l++) sum += Dinv[bc_oddlist[f][l]+added_evenlist[b]*VOLUME/2]*Ginv[f][l+a*n_bc];
      BAP[f][a][added_evenlist[b]] = sum;
    }
    newsite[f][added_evenlist[b]] = 0;
  }

  for( int a=0; a<n_added_odd; a++) if( newsite[f][VOLUME/2+added_oddlist[a]] == 1 ){
    for( int b=0; b<n_bc; b++)  { 
      double sum = 0;
      for( int k=0; k<n_bc; k++) sum += Ginv[f][b+k*n_bc]*Dinv[added_oddlist[a]+bc_evenlist[f][k]*VOLUME/2];
      PAC[f][added_oddlist[a]][b] = -sum;
    }
    newsite[f][VOLUME/2+added_oddlist[a]] = 0;
  }


  for( int b=0; b<n_added_even; b++) for( int a=0; a<n_added_odd; a++) 
  //if( new_combination[f][added_evenlist[b]][added_oddlist[a]] == 1 )
  {
    double sum = 0;
    for( int l=0; l<n_bc; l++)
      sum += Dinv[bc_oddlist[f][l]+added_evenlist[b]*VOLUME/2]*PAC[f][added_oddlist[a]][l];
    BGC[f][b][a] = Dinv[added_oddlist[a]+added_evenlist[b]*VOLUME/2] + sum;

    //new_combination[f][added_evenlist[b]][added_oddlist[a]] = 0;
  }


  /* Pick entries for the fluctuation matrix */
  double *F;
  F = malloc(n_fl*n_fl*sizeof(double));
  for( int b=0; b<n_removed_odd; b++) for( int a=0; a<n_removed_even; a++)
    F[a*n_fl+b] = Ginv[f][removed_oddlist[b]+removed_evenlist[a]*n_bc];
  for( int a=0; a<n_added_odd; a++) for( int b=0; b<n_removed_odd; b++)
    F[(n_removed_even+a)*n_fl+b] = PAC[f][added_oddlist[a]][removed_oddlist[b]];
  for( int a=0; a<n_removed_even; a++) for( int b=0; b<n_added_even; b++)
    F[a*n_fl+(n_removed_odd+b)] = BAP[f][removed_evenlist[a]][added_evenlist[b]];
  for( int b=0; b<n_added_even; b++) for( int a=0; a<n_added_odd; a++)
    F[(n_removed_even+a)*n_fl+(n_removed_odd+b)] = BGC[f][b][a];

  /* determinant from LU */
  double det=1;
  int ipiv[n_fl], info;
  LAPACK_dgetrf( &n_fl, &n_fl, F, &n_fl, ipiv, &info );
  double logdet = 0;
  for(int i=0; i<n_fl; i++) {
    logdet += log(fabs(F[i*n_fl+i]));
    //if( ipiv[i] != i+1 ) det*=-1; 
  }
  det *= exp(logdet);


  new_determinant[f] = det*det;
#ifdef DEBUG
  printf("Det %d %g  %g %g  %g\n", ndet, det, new_determinant[f], current_determinant[f], new_determinant[f]/current_determinant[f] );
  ndet++;
  if( new_determinant[f]/current_determinant[f] > 1e50 ){
    for(int b=0; b<n_bc_monomers[f]/2; b++) {
      for(int a=0; a<n_bc_monomers[f]/2; a++) {
        printf(" %5.2f", Ginv[f][b+n*a]);
      }
      printf("\n");
    }
    printf("\n");

    for(int ie=0; ie<n_fl; ie++) {
      for(int io=0; io<n_fl; io++) {
        printf(" %8.5f", F[ie+n_fl*io]);
      }
      printf("\n");
    }
    printf("\n");
    exit(1);
  }
#endif

  free(F);
  free(removed_evenlist);
  free(removed_oddlist);
  free(added_evenlist);
  free(added_oddlist);

  return( new_determinant[f]/current_determinant[f] );
}


#endif //FLUCTUATION_DETERMINANT



/* Calculate weight for specific operations simply by changing the field
   and calling determinant().
   Unfortunately the background needs to be updated here when necessary. */
#ifndef ANTIPERIODIC_MASS
double det_add_monomers(int x1, int x2, int do_flavor[N_FLAVOR]){

  #ifdef FLUCTUATION_DETERMINANT
  static int init = 1;
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif
  
  //current_determinant[0] = current_determinant[1] = 1;
  double det = 1;//./(determinant( 0 )*determinant( 1 ));
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    occupation_field[f][x1] = OCCUPIED;
    occupation_field[f][x2] = OCCUPIED;
    n_occupied[f] += 2;
  //}
    det *= determinant(f);// determinant( 0 )*determinant( 1 );
  //for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    occupation_field[f][x1] = 0;
    occupation_field[f][x2] = 0;
    n_occupied[f] -= 2;
  }

  return det;
}

double det_move_monomers(int x1, int x2, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  static int init = 1;
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif

  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    occupation_field[f][x1] = 0;
    occupation_field[f][x2] = OCCUPIED;
    det *= determinant( f );
    occupation_field[f][x1] = OCCUPIED;
    occupation_field[f][x2] = 0;
  }

  return det;
}


double det_remove_monomers(int x1, int x2, int do_flavor[N_FLAVOR]){

  #ifdef FLUCTUATION_DETERMINANT
  static int init = 1;
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif

  //current_determinant[0] = current_determinant[1] = 1;
  double det = 1;//./(determinant( 0 )*determinant( 1 ));
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    occupation_field[f][x1] = 0;
    occupation_field[f][x2] = 0;
    n_occupied[f] -= 2;
  //}
    det *= determinant(f);//determinant( 0 )*determinant( 1 );
  //for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    occupation_field[f][x1] = OCCUPIED;
    occupation_field[f][x2] = OCCUPIED;
    n_occupied[f] += 2;
  }

  return det;
}

#else

double det_add_monomers(int x1, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  static int init = 1;
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif
  
  //current_determinant[0] = current_determinant[1] = 1;
  double det = 1;
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]) {
    occupation_field[f][x1] = OCCUPIED;
    n_occupied[f] += 1;
    det *= fabs(determinant(f));
    occupation_field[f][x1] = 0;
    n_occupied[f] -= 1;
  }

  return det;
}

double det_remove_monomers(int x1, int do_flavor[N_FLAVOR]){
  #ifdef FLUCTUATION_DETERMINANT
  static int init = 1;
  if( init ) {
    for( int f=0; f<N_FLAVOR; f++ ) n_fluctuations[f] = 0;
    init = 0;
  }
  for( int f=0; f<N_FLAVOR; f++ ) if( n_fluctuations[f] > max_fluctuations - 2 ){
    update_background( f );
  }
  #endif

  //current_determinant[0] = current_determinant[1] = 1;
  double det = 1;//./(determinant( 0 )*determinant( 1 ));
  for(int f=0; f<N_FLAVOR; f++) if(do_flavor[f]){
    occupation_field[f][x1] = 0;
    n_occupied[f] -= 1;
    det *= fabs(determinant(f));
    occupation_field[f][x1] = OCCUPIED;
    n_occupied[f] += 1;
  }

  return det;
}
#endif
