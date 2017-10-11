/****************************************************************
 * Invert the fermion matrix using the conjugate gradient algorithm
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
extern double *D;
extern double *Dinv;


/* Map indexes of even and odd sites to global indexes */
extern int * even_to_global_index;
extern int * odd_to_global_index;
extern int * global_to_even_index;
extern int * global_to_odd_index;


void init_evenodd_index( );
void write_Dinv();
int block_index( int v[ND] );



/* Save the propagator matrix to avoid recalculation when possible */
void write_Dinv( int block ){
  FILE * Dinv_file;
  char filename[100] = "Dinv";
  int n = VOLUME/2, nsites;
  if(block==1) {
    int ns = 1;
    for( int nu=0; nu<ND; nu++) ns*=2;
    nsites = n*ns;
  } else {
    nsites = n*n;
  }

  Dinv_file = fopen(filename,"wb");
  if (Dinv_file){
    fwrite(Dinv, nsites, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {
    printf("Could not write Dinv file\n");
    exit(1);
  }
}






#define CG_ACCURACY 1e-10
#define CG_MAX_ITER 100

void cg( double *source ){

  double rr, pMp, a;
  double *r = malloc(VOLUME*sizeof(double));
  double *p = malloc(VOLUME*sizeof(double));
  double *Mp = malloc(VOLUME*sizeof(double));
  double *MMp = malloc(VOLUME*sizeof(double));

  Dv( r, source );

  for(int x=0; x<VOLUME; x++) p[x] = r[x];
  for(int x=0; x<VOLUME; x++) source[x] = 0;

  double rr_old = 0;
  for(int x=0; x<VOLUME; x++) rr_old += r[x]*r[x];
  double rr_init = rr_old;
  if( rr_old < CG_ACCURACY ){
    return ;
  }
  //printf("CG: %d  %g\n",0,rr_old);

  int k;
  for( k=1; k < CG_MAX_ITER; k++ )
  {
    Dv( Mp, p );
    Dv( MMp, Mp );
    pMp = 0;
    for(int x=0; x<VOLUME; x++) pMp += p[x]*MMp[x];
    a = rr_old / pMp ;
    for(int x=0; x<VOLUME; x++) source[x] = source[x] + p[x]*a;
    for(int x=0; x<VOLUME; x++) r[x] = r[x] - MMp[x]*a;
    
    rr=0;
    for(int x=0; x<VOLUME; x++) rr += r[x]*r[x];
    //printf("CG: %d  %g %g %g %g\n",k,rr,rr_old,pMp,a);
    if( rr < CG_ACCURACY )
      break;

    double b = rr / rr_old ;
    for(int x=0; x<VOLUME; x++) p[x] = r[x] + p[x]*b;

    rr_old = rr;
  }

  free(r);
  free(p);
  free(Mp);
  free(MMp);
  printf("CG: %d %g %g\n",k,rr/rr_init,rr_init);
  
}





/* Initialize the propagator matrix, from the file if possible */
/* Use conjugate gradient otherwise */
void init_Dinv_cg(){
  
  init_evenodd_index();

  // Try loading from the file
  FILE * Dinv_file;
  char filename[100] = "Dinv";
  int n = VOLUME/2;

  //Number of sites in a staggered block
  int ns = 1;
  for( int nu=0; nu<ND; nu++) ns*=2;

  Dinv = malloc(n*ns*sizeof(double));

  Dinv_file = fopen(filename,"rb");
  /*if (Dinv_file){
    printf("Reading Dinv\n");
    fread(Dinv, n*ns, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {*/

    /* File not found or loading failed */
    printf("No Dinv file\n");

    double * source = malloc(VOLUME*sizeof(double));
    
    // Loop over necessary sources
    for( int s=0; s<ns; s++){

      int ie = s;
      int v[ND], sum=0;
      for( int nu=0; nu<ND; nu++){
        v[nu] = ie%2; ie/=2; sum+=v[nu];
      }
      if( sum%2 == 0 ){
        printf("(%d,%d,%d,%d)\n",v[0],v[1],v[2],v[3]);
        int ig = site_vector_to_index(v);
        ie = block_index(v);
   
    //for(int ie=0; ie<VOLUME/2; ie++){
        for(int x=0; x<VOLUME; x++) source[x] = 0;
        //int ig = even_to_global_index[ie];
        source[ig] = 1;

        /* Now we have a source with index ig set to 1
         * Find the inverse */
        cg( source );
        int v[ND];
        site_index_to_vector(ig,v);
        v[0] += 1;
        ig = site_vector_to_index(v);

        for(int x=0; x<VOLUME; x++){
          int io = global_to_odd_index[x];
          if( io != -1 ) Dinv[ io + ie*VOLUME/2 ] = source[x];
        }
      }
    }
    
    free(source);

    //Save to file
    //write_Dinv( 1 );
  //}

  printf("Propagator matrix initialized\n");

}















#if 0

/* Initialize the propagator matrix, from the file if possible */
/* Deprecated: Requires the construction of the full even-odd matrix.
 * Instead we invert on a staggered hypercube and use translations to
 * find the other propagators. */

double Df( int x1, int x2 );

void init_Dinv_deprecated(){
  
  init_evenodd_index();

  // Try loading from the file
  FILE * Dinv_file;
  char filename[100] = "Dinv";
  int n = VOLUME/2;
  Dinv = wrapmalloc(n*n*sizeof(double));

  Dinv_file = fopen(filename,"rb");
  if (Dinv_file){
    printf("Reading Dinv\n");
    fread(Dinv, n*n, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {

    /* File not found or loading failed */
    printf("No Dinv file\n");
    //init_fermion_matrix( );
    //Save memory by constructing the Dirac operator on the run
    
    int info;
    int *ipiv;
    ipiv = wrapmalloc( n*sizeof(int) );
    for(int i=0; i<n*n; i++) Dinv[i] = 0;

    /* Construct the even to odd matrix */
    for(int ie=0; ie<n; ie++) for(int io=0; io<n; io++){
      Dinv[ie+n*io] = Df(even_to_global_index[ie],odd_to_global_index[io]);
    }

    /* Invert using lapack routines */
    int lwork=n*n;
    double *work;
    work = wrapmalloc( lwork*sizeof(double) );
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
    //Save to file
    //write_Dinv(0);
  }

  printf("Propagator matrix initialized\n");

}

#endif








