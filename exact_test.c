/****************************************************************
 * An exact calculation by enumerating all configurations
 *
 ****************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "Staggered.h"

#ifdef FLUCTUATION_DETERMINANT
fluctuation determinant not available for exact test
use full determinant
#endif

/* storage */
int    **eta;   //Staggered eta matrix
int    **ksi;   //Staggered ksi matrix for link mass
int    *occupied_sites;   //List of occupied sites
double U;
double m;
double linkmass;

/* Neighbour index arrays
 */
int **neighbour;


/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_occupied[N_FLAVOR];
int **occupation_field;

int n_fourfermion_monomer;
int n_mass_monomer[N_FLAVOR];
int *fourfermion_monomer;
int **mass_monomer;

static inline int opp_dir(int dir){
  return ( dir + ND ) % NDIRS;
}


int site_vector_to_index( int vector[] ){

  int index = vector[ND-1];
  for (int dir=ND-1; dir--;) index = index*Ldim[dir] + vector[dir];

  return index;
}

void site_index_to_vector(int index, int vector[] ){
  int i = index;

  for (int dir=0; dir<ND; dir++){
    vector[dir] = i%Ldim[dir];
    i = i/Ldim[dir];
  }

}




double Z=0;
double ff_density=0;
double m_density=0;
double W = 1;
extern double current_determinant[N_FLAVOR];
void all_configs( int x0 ){
    int max=VOLUME+1, min=0;
    
    //Do measurements on the current config
    current_determinant[0] = current_determinant[1] = 1;


      double det = determinant(0)*determinant(1);
      Z += W*det;
      ff_density += W*det*n_fourfermion_monomer;
      m_density += W*det*n_mass_monomer[0];
    
    //printf("measurements %g %d %d %d\n",det,n_fourfermion_monomer,n_mass_monomer[0],n_mass_monomer[1]);

    for(int x=x0; x<VOLUME; x++){
        double W_old = W;
        
        //Make it a four fermion monomer
        fourfermion_monomer[x] = 1;
        occupation_field[0][x] = 1;
        occupation_field[1][x] = 1;
        n_fourfermion_monomer++;
        n_occupied[0]++;
        n_occupied[1]++;
        W*=U;
        all_configs( x+1 );
        W=W_old;
        fourfermion_monomer[x] = 0;
        occupation_field[0][x] = 0;
        occupation_field[1][x] = 0;
        n_fourfermion_monomer--;
        n_occupied[0]--;
        n_occupied[1]--;

#ifndef MASS_IN_MATRIX     
        //Make it a mass monomer
        mass_monomer[0][x] = 1;
        occupation_field[0][x] = 1;
        n_mass_monomer[0]++;
        n_occupied[0]++;
        W*=m;
        all_configs( x+1 );
        
        mass_monomer[1][x] = 1;
        occupation_field[1][x] = 1;
        n_mass_monomer[1]++;
        n_occupied[1]++;
        W*=m;
        all_configs( x+1 );
        
        mass_monomer[0][x] = 0;
        occupation_field[0][x] = 0;
        n_mass_monomer[0]--;
        n_occupied[0]--;
        W/=m;
        all_configs( x+1 );
        
        mass_monomer[1][x] = 0;
        occupation_field[1][x] = 0;
        n_mass_monomer[1]--;
        n_occupied[1]--;
        W = W_old;
#endif
    }
}




void print_config()
{
  int s[ND];
  for ( s[0]=0; s[0]<Ldim[0]; s[0]++) {
    for ( s[1]=0; s[1]<Ldim[1]; s[1]++){
#if ND==4
      for ( s[2]=0; s[2]<Ldim[2]; s[2]++) {
        for ( s[3]=0; s[3]<Ldim[3]; s[3]++){
#endif
#if ND==3
      for ( s[2]=0; s[2]<Ldim[2]; s[2]++) {
#endif
          int empty = 1;
          int index = site_vector_to_index(s);
          if(fourfermion_monomer[index]==1)  { empty = 0; printf(" o "); }
          if(mass_monomer[0][index]==1)  { empty = 0;
            if(mass_monomer[1][index]==1) printf(" + ");
            else printf(" - "); } 
          else if(mass_monomer[1][index]==1)  { empty = 0; printf(" | "); }
          if(empty==1 && occupation_field[0][index] == 1 ) { empty = 0;
               if( occupation_field[1][index]==1) printf(" x ");
               else printf(" \\ "); }
          if(empty==1 && occupation_field[1][index] == 1 ) { empty = 0; printf(" / "); }
          if(empty==1) { printf(" . "); }
#if ND==4
        }
        printf(" ");
      }
      printf(" \n");
#endif
#if ND==3
      }
      printf(" ");
#endif
    }
    printf(" \n");
  }
  printf(" \n");
  //usleep(1000000);
}





/* Main function
 */
int main(int argc, char* argv[])
{
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  int i,n_updates_per_measure,n_mass_updates_per_measure,n_measure,n_average;
  long seed;

  /* Read in the input */
  printf(" Discarded : ");
  scanf("%d",&n_updates_per_measure);

  printf(" Discarded : ");
  scanf("%d",&n_measure);

  printf(" Discarded : ");
  scanf("%d",&n_mass_updates_per_measure);
  
  printf(" Discarded : ");
  scanf("%d",&n_average);

  int max_fluctuations;
  printf(" Discarded : ");
  scanf("%d",&max_fluctuations);

  printf(" U : ");
  scanf("%lf",&U);
 
  printf(" m : ");
  scanf("%lf",&m);

  printf(" linkmass : ");
  scanf("%lf",&linkmass);

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" %dD staggered fermion, ( ",ND);
  for( int nu=0; nu<ND-1; nu++){
    printf("%d , ", Ldim[nu] );
  }
  printf("%d ) lattice\n", Ldim[ND-1]);
  printf(" U %f \n", U);
  printf(" m %f %f \n", m, linkmass);
  fflush(stdout);

  VOLUME = 1;
  for( int nu=0; nu<ND; nu++) VOLUME*=Ldim[nu];

  occupation_field = malloc( N_FLAVOR*sizeof(int*) );
  for( int i=0; i<N_FLAVOR; i++ ){
    occupation_field[i] = malloc( VOLUME*sizeof(int) );
    for (int x=0; x<VOLUME; x++) occupation_field[i][x] = 0;
    n_occupied[i] = 0;
  }

  mass_monomer = malloc( N_FLAVOR*sizeof(int*) );
  for( int i=0; i<N_FLAVOR; i++ ){ 
    mass_monomer[i] = malloc( VOLUME*sizeof(int) );
    for (int x=0; x<VOLUME; x++) mass_monomer[i][x] = 0;
    n_mass_monomer[i] = 0;
  }
  fourfermion_monomer = malloc( VOLUME*sizeof(int) );
  for (int x=0; x<VOLUME; x++) fourfermion_monomer[x] = 0;
  n_fourfermion_monomer = 0;

  /* The staggered eta matrix */
  eta = malloc( ND*sizeof(int *) );
  for( int nu=0; nu<ND; nu++) eta[nu] = malloc( VOLUME*sizeof(int) );
  for (int x=0; x<VOLUME; x++) {
    int vector[ND];
    site_index_to_vector(x,vector);
    for( int nu=0; nu<ND; nu++) {
      int eta_exponent = 0;
      for( int mu=0; mu<nu; mu++){
         eta_exponent += vector[mu];
      }
      if( eta_exponent%2 == 0 ){
        eta[nu][x] =  1;
      } else {
        eta[nu][x] = -1;
      }
    }
  }

  /* The staggered eps*ksi matrix for link mass */
  ksi = malloc( ND*sizeof(int *) );
  for (int nu=0; nu<ND; nu++) ksi[nu] = malloc( VOLUME*sizeof(int) );
  for (int x=0; x<VOLUME; x++) {
    int vector[ND];
    site_index_to_vector(x,vector);
    int eps_exponent = 0;
    for( int nu=0; nu<ND; nu++) eps_exponent += vector[nu];
    for( int nu=0; nu<ND; nu++) {
      int ksi_exponent = 0;
      for( int mu=nu+1; mu<ND; mu++)
         ksi_exponent += vector[mu];
      if( (eps_exponent+ksi_exponent)%2 == 0 ){
        ksi[nu][x] =  1;
      } else {
        ksi[nu][x] = -1;
      }
    }
  } 


  /* The nb array */
  neighbour = malloc( NDIRS*sizeof(int *) );
  for( int nu=0; nu<NDIRS; nu++) neighbour[nu] = malloc( VOLUME*sizeof(int) );
  for (int x=0; x<VOLUME; x++) {
     int vector[ND];
     site_index_to_vector(x,vector);
     for( int nu=0; nu<ND; nu++) {
        int nb[ND];
        for( int mu=0; mu<ND; mu++) nb[mu] = vector[mu];
        nb[nu] = (nb[nu]+1)%Ldim[nu];
        neighbour[nu][x] = site_vector_to_index(nb);

        for( int mu=0; mu<ND; mu++) nb[mu] = vector[mu];
        nb[nu] = (nb[nu]+Ldim[nu]-1)%Ldim[nu];
        neighbour[opp_dir(nu)][x] = site_vector_to_index(nb);
     }
  } 


  
  /* calculate propagators */
#ifdef PROPAGATOR_MATRIX
  calc_Dinv( );
#endif


  struct timeval start, end;
  gettimeofday(&start,NULL);

  /* Background and fluctuation matrix */
#ifdef FLUCTUATION_DETERMINANT
  for(int m=0; m<N_FLAVOR; m++) update_background( m );
#endif

  all_configs( 0 );

  printf("MONOMERDENSITY %g %g\n", ff_density/(Z*VOLUME), m_density/(Z*VOLUME) );

  
  printf(" ** simulation done\n");

  return(1);
}



























