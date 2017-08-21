/****************************************************************
 * Simulate 2 flavors of staggered fermions with a four fermion
 * interaction using the fermion bag algorithm (arXiv:0910.5736).
 * Measure the bilinear condensate using a mass term as an
 * external source , similarly to arxiv 1609.08541.
 ****************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "Staggered.h"

/* storage */
int    **eta;   //Staggered eta matrix
int    **ksi;   //Staggered ksi matrix for link mass
int    *occupied_sites;   //List of occupied sites
double U;
double m;
double linkmass;

/* Maximum number of fluctuations from the background configuration */
int max_fluctuations;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_occupied[N_FLAVOR];
int **occupation_field;

int n_fourfermion_monomer;
int n_mass_monomer[N_FLAVOR];
int *fourfermion_monomer;
int **mass_monomer;

/* Neighbour index arrays
 */
int **neighbour;

/* Opposite of a direction */
static inline int opp_dir(int dir){
  return ( dir + ND ) % NDIRS;
}

void add_monomer( int x, int x2);
void remove_monomer( int x, int x2 );




/* Write the current configuration into a file */
void write_config(){
  FILE * config_file;
  char filename[100];
  sprintf(filename, "config_checkpoint_U%.6gm%.6g",U,m);

  int * buffer = malloc((1+2*N_FLAVOR)*VOLUME*sizeof(int));
  for (int i=0; i<VOLUME; i++) {
    buffer[(1+2*N_FLAVOR)*i] =  fourfermion_monomer[i];
    for (int f=0; f<N_FLAVOR; f++)
      buffer[(1+2*N_FLAVOR)*i+1+f] =  mass_monomer[f][i];
    for (int f=0; f<N_FLAVOR; f++)
      buffer[(1+2*N_FLAVOR)*i+3+f] =  occupation_field[f][i];
  }
  
  config_file = fopen(filename,"wb");
  if (config_file){
    fwrite(buffer, (1+2*N_FLAVOR)*VOLUME, sizeof(int), config_file);
    fclose(config_file);
  } else {
    printf("Could not write configuration\n");
    exit(1);
  }
  free(buffer);
}

/* Read the configuration file */
void read_config(){
  FILE * config_file;
  char filename[100];
  sprintf(filename, "config_checkpoint_U%.6gm%.6g",U,m);

  config_file = fopen(filename,"rb");
  int * buffer = malloc((1+2*N_FLAVOR)*VOLUME*sizeof(int));
  if (config_file){
    fread(buffer, (1+2*N_FLAVOR)*VOLUME, sizeof(int), config_file);
    fclose(config_file);
  
    for (int i=0; i<VOLUME; i++) {
      fourfermion_monomer[i] = buffer[(1+2*N_FLAVOR)*i];
      for (int f=0; f<N_FLAVOR; f++)
        mass_monomer[f][i] = buffer[(1+2*N_FLAVOR)*i+1+f];
      for (int f=0; f<N_FLAVOR; f++)
        occupation_field[f][i] = buffer[(1+2*N_FLAVOR)*i+3+f];
    }
  } else {
    printf("No configuration file\n");
    printf("Starting from an empty configuration\n");
  }

  free(buffer);
}







/* Map ND-vector to an index */
int site_vector_to_index( int vector[] ){

  int index = vector[ND-1];
  for (int dir=ND-1; dir--;) index = index*Ldim[dir] + vector[dir];

  return index;
}

/* Map site index to an ND-vector */
void site_index_to_vector(int index, int vector[] ){
  int i = index;

  for (int dir=0; dir<ND; dir++){
    vector[dir] = i%Ldim[dir];
    i = i/Ldim[dir];
  }

}







/* Turn a monomer on at a link */
void occupy_site(int x1, int x2, int m){
  if ( occupation_field[m][x1] == UNOCCUPIED && occupation_field[m][x2] == UNOCCUPIED ){
    occupation_field[m][x1] = OCCUPIED;
    occupation_field[m][x2] = OCCUPIED;
#ifdef DEBUG
    printf("Turned on flavor %d at %d and %d\n",m,x1,x2);
#endif
  } else {
    printf("Site already occupied\n");
    exit(1);
  }
}

/* Turn a link off */
void free_site(int x1, int x2, int m){
  if ( occupation_field[m][x1] != UNOCCUPIED && occupation_field[m][x2] != UNOCCUPIED ){
    occupation_field[m][x1] = UNOCCUPIED;
    occupation_field[m][x2] = UNOCCUPIED;
#ifdef DEBUG
    printf("Turned off flavor type %d at %d and %d\n",m,x1,x2);
#endif
  } else {
    printf("Flavor %d at sites %d and %d already unoccupied ( %d %d )\n",m,x1,x2,occupation_field[m][x1],occupation_field[m][x2]);
    exit(1);
  }
}






#ifndef MASS_IN_MATRIX
/* Measure the susceptibilities using a worm update method
 * aa is measures the same flavor susceptibility and 
 * ab measures the mixed flavor susceptibility
 *
 * The worm begins by creating two sources with the probability 
 * from a partition function with two insertions. One of the sources 
 * moves on the lattice, according to the same probability density,
 * changing the configuration on the way. The number of configurations
 * counts the weight of the modified partition function against the
 * unmodified one, that is, the measurement. 
 */
int n_susc1_measurements=1000;
void susceptibility_aa(){
  int step=0;
  /* Choose flavor */
  int f = (int) (mersenne()*N_FLAVOR);
  int flavorlist[N_FLAVOR] = {0,0};
  flavorlist[f] = 1;

  for( int mea=0; mea<n_susc1_measurements; mea++ ){
    int x1,x2;
    x1 = (int) (mersenne()*VOLUME);
    int nu = (int) (mersenne()*NDIRS);
    x2 = neighbour[nu][x1];
    
    /* Sites needs to be unoccupied */
    if ( (occupation_field[f][x1] == UNOCCUPIED) && (occupation_field[f][x2] == UNOCCUPIED) ) {

      int notdone = 0;
      /* Try to turn on two monomers */
      double p = det_add_monomers( x1, x2, flavorlist );
      if( mersenne() < p ){
        occupation_field[f][x1]=1;
        occupation_field[f][x2]=1;
        update_current_determinant(flavorlist);
        notdone = 1;
      }

      /* Move the monomer */
      while(notdone) { 
        step++;
        int dir, newx;

        dir = (int) (mersenne()*NDIRS);
        newx = neighbour[dir][x2];
        if( mersenne() < 0.5 ) { 
          if( newx==x1 ){
            double p = det_remove_monomers( x1, x2, flavorlist );
            if( mersenne() < p ){
              occupation_field[f][x1]=UNOCCUPIED;
              occupation_field[f][x2]=UNOCCUPIED;
              update_current_determinant(flavorlist);
              notdone = 0;
            }
          } else if( mass_monomer[f][newx] == 1 ) {
            mass_monomer[f][newx] = UNOCCUPIED;
            mass_monomer[f][x2] = OCCUPIED;
            x2 = newx;
          }
        } else {
          
          dir = (int) (mersenne()*NDIRS);
          newx = neighbour[dir][newx];
          if(occupation_field[f][newx] == UNOCCUPIED) {

            double p = det_move_monomers( x2, newx, flavorlist );
            if( mersenne() < p ) {
              occupation_field[f][x2]=UNOCCUPIED;
              occupation_field[f][newx]=OCCUPIED;
              update_current_determinant(flavorlist);
              x2=newx;
            }
          }
        }
      }
    }
  }
  printf("SUSC1 %g\n",(double)step/(double)(4*n_susc1_measurements));
}

int n_susc2_measurements=100;
void susceptibility_ab(){
  int step=0;

  if( n_fourfermion_monomer > 0 ) for( int mea=0; mea<n_susc2_measurements; mea++ ){
    int x1;
    do { /* Site needs to be fully occupied */
      x1 = (int) (mersenne()*VOLUME);
    } while ( fourfermion_monomer[x1] == UNOCCUPIED );
    /* Turn the monomer into a pair of sources */
    fourfermion_monomer[x1] = UNOCCUPIED;
    //Occupation field remains 1
    int x2=x1;
    int notdone = 1;

    /* Choose flavor to move */
    int f = (int) (mersenne()*N_FLAVOR);
      
    /* Move the monomer */
    while(notdone) { //print_config();
      if( x1!=x2 ) step++;
      int dir, newx;

        dir = (int) (mersenne()*NDIRS);
        newx = neighbour[dir][x2];
        if( mass_monomer[f][newx] == OCCUPIED ) {
          mass_monomer[f][newx] = UNOCCUPIED;
          mass_monomer[f][x2] = OCCUPIED;
          x2 = newx;
        } else {
        
        dir = (int) (mersenne()*NDIRS);
        newx = neighbour[dir][newx];
        if(occupation_field[f][newx] == UNOCCUPIED) {
          int flavorlist[N_FLAVOR] = {0,0};
          flavorlist[f] = 1;
          double p = det_move_monomers( x2, newx, flavorlist );
          if( mersenne() < p ) {
            occupation_field[f][x2]=UNOCCUPIED;
            occupation_field[f][newx]=OCCUPIED;
            update_current_determinant(flavorlist);
            x2=newx;
          }
        }
      }

      if( x1==x2 ){
        notdone=0;
        fourfermion_monomer[x1] = OCCUPIED;
      }
     
    }
  }
  printf("SUSC2 %g\n",(double)step*n_fourfermion_monomer/(double)(2*U*n_susc2_measurements*VOLUME));
}
#endif





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
          if(fourfermion_monomer[index]==OCCUPIED)  { empty = 0; printf(" o "); }
          if(mass_monomer[0][index]==OCCUPIED)  { empty = 0;
            if(mass_monomer[1][index]==OCCUPIED) printf(" + ");
            else printf(" - "); } 
          else if(mass_monomer[1][index]==OCCUPIED)  { empty = 0; printf(" | "); }
          if(empty==1 && occupation_field[0][index] == OCCUPIED ) { empty = 0;
               if( occupation_field[1][index]==OCCUPIED) printf(" x ");
               else printf(" \\ "); }
          if(empty==1 && occupation_field[1][index] == OCCUPIED ) { empty = 0; printf(" / "); }
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

  int i,n_measure,n_average;
  long seed;

  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_measure);
  
  printf(" Averaged over : ");
  scanf("%d",&n_average);

  printf(" Size of the fluctuation matrix : ");
  scanf("%d",&max_fluctuations);

  printf(" U : ");
  scanf("%lf",&U);
 
  printf(" m : ");
  scanf("%lf",&m);

  printf(" linkmass : ");
  scanf("%lf",&linkmass);

  printf(" Random number : ");
  scanf("%ld",&seed);
  seed_mersenne( seed );

  /* "Warm up" the rng generator */
  for (i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" %dD staggered fermion, ( ",ND);
  for( int nu=0; nu<ND-1; nu++){
    printf("%d , ", Ldim[nu] );
  }
  printf("%d ) lattice\n", Ldim[ND-1]);
  printf(" %d updates\n", n_measure );
  printf(" averaged over %d\n", n_average );
  printf(" fluctuation matrix size %d\n", max_fluctuations );
  printf(" U %f \n", U);
  printf(" m %f %f \n", m, linkmass);
  printf(" Random seed %ld\n", seed );
  fflush(stdout);

  /* Some usefull definitions and necessary fields */
  VOLUME = 1;
  for( int nu=0; nu<ND; nu++) VOLUME*=Ldim[nu];

  /* field marking fields as occupied of unoccupied */
  occupation_field = malloc( N_FLAVOR*sizeof(int*) );
  for( int i=0; i<N_FLAVOR; i++ ){
    occupation_field[i] = malloc( VOLUME*sizeof(int) );
    for (int x=0; x<VOLUME; x++) occupation_field[i][x] = 0;
    n_occupied[i] = 0;
  }

  /* field for marking four fermion and mass monomers */
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


  /* The neighbour array */
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

  /* Setup starting configuration */
  read_config();
  for (int x=0; x<VOLUME; x++) if( fourfermion_monomer[x] == 1 ) 
    n_fourfermion_monomer++;
  for (int x=0; x<VOLUME; x++) for(int m=0; m<N_FLAVOR; m++) if( mass_monomer[m][x] == 1 ) 
    n_mass_monomer[m]++;
  for (int x=0; x<VOLUME; x++) for(int m=0; m<N_FLAVOR; m++) if( occupation_field[m][x] == 1 ) 
    n_occupied[m]++;


  /* Start timing */ 
  struct timeval start, end;
  gettimeofday(&start,NULL);

  /* Initialize background and fluctuation matrix */
#ifdef FLUCTUATION_DETERMINANT
  for(int m=0; m<N_FLAVOR; m++) update_background( m );
#endif

  /* The update/measure loop */
  for (i=1; i<n_measure+1; i++) {

    /* Zero measurements */
    double fourfermion_monomer_density=0;
    double mass_monomer_density[N_FLAVOR];
    double linkvev = 0;
    for(int m=0; m<N_FLAVOR; m++) mass_monomer_density[m] = 0;

    /* For keeping track of changes */
    int additions = 0, removals = 0, m_additions = 0, m_removals = 0, switches = 0;
    #ifdef LOCAL_UPDATE
    int moves = 0, m_moves = 0;
    #endif

    /* updates */
    for(int j=0; j<n_average; j++){
      /* Update */
      #ifdef WORM_UPDATE
      worm_update( &additions, &removals, &m_additions, &m_removals, &switches );
      #endif
      #ifdef LOCAL_UPDATE
      local_update( &additions, &removals, &moves, &m_additions, &m_removals, &m_moves, &switches );
      #endif 

      /* sum measurements */
      fourfermion_monomer_density += n_fourfermion_monomer/(double)VOLUME;
      for(int m=0; m<N_FLAVOR; m++) mass_monomer_density[m] += n_mass_monomer[m]/(double)VOLUME;

      linkvev += link_vev();
    }

    

    /* Time and report */
    gettimeofday(&end,NULL);
    unsigned long diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
    #ifdef LOCAL_UPDATE
    printf("\n%d, %d updates in %.3g seconds, %d additions, %d removals, %d moves \n", i, n_average, 1e-6*diff, additions, removals, moves );
    printf("%d, %d mass updates, %d additions, %d removals, %d moves, %d switches \n", i, n_average, m_additions, m_removals, m_moves, switches );
    #endif
    #ifdef WORM_UPDATE
    printf("\n%d, %d updates in %.3g seconds, %d additions, %d removals, %d mass additions, %d mass removals, %d switches\n", i, n_average, 1e-6*diff, additions, removals, m_additions, m_removals, switches );
    #endif
 
    /* Basic measurements, monomer densities */
    printf("MONOMERDENSITY %g ", fourfermion_monomer_density/(double)n_average);
    for(int m=0; m<N_FLAVOR; m++) printf("%g ", mass_monomer_density[m]/(double)n_average);
    printf("\n");
    printf("LINKVEV %g \n", linkvev/(double)(n_average*VOLUME*2));
    
    /* Susceptibilities */
    //susceptibility_ab();
    //susceptibility_aa();
    diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
    printf("Inversion in %.3g seconds\n", i, n_average, 1e-6*diff );
    /* write checkpoint */
    write_config();

    gettimeofday(&start,NULL);
  }

  printf(" ** simulation done\n");

  return(1);
}



























