/****************************************************************
 * Worm update 
 * Usually more efficient for adding monomer
 ****************************************************************/

#ifdef DEBUG
#include <fenv.h>
#endif

#include "Staggered.h"

#ifdef WORM_UPDATE

extern double U;
extern double m;
extern double linkmass;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
extern int n_occupied[N_FLAVOR];
extern bool **occupation_field;

extern int n_fourfermion_monomer;
extern int n_mass_monomer[N_FLAVOR];
extern bool *fourfermion_monomer;
extern bool **mass_monomer;

/* Neighbour index arrays
 */
extern int **neighbour;


/* Try to turn a four fermion monomer into two mass monomers
 * or a mass monomer into two fourfermion monomers
 * This does not require calculating a determinant
 */
#define N_SWITCH_ATTEMPS VOLUME
int switch_monomers()
{
  int success = 0, x1;
  for( int attempt=0; attempt < N_SWITCH_ATTEMPS; attempt++) { 
  x1 = mersenne()*VOLUME;
  if ( (occupation_field[0][x1] == 1) &&  (occupation_field[1][x1] == 1) ) {
   
    if( fourfermion_monomer[x1] == 1){
     double p = m*m/U ;
     if( mersenne() < p ) {
        fourfermion_monomer[x1]=0;
        mass_monomer[0][x1]=1;
        mass_monomer[1][x1]=1;
        n_fourfermion_monomer-=1;
        n_mass_monomer[0]+=1;
        n_mass_monomer[1]+=1;
        success = 1;
     }
    } else if( fourfermion_monomer[x1] == 0 ){
     double p = U/(m*m) ;
     if( mersenne() < p ) {
        fourfermion_monomer[x1]=1;
        mass_monomer[0][x1]=0;
        mass_monomer[1][x1]=0;
        n_fourfermion_monomer+=1;
        n_mass_monomer[0]-=1;
        n_mass_monomer[1]-=1;
        success = 1;
     }
    }
    
  }
  }
  return success;
}



double worm_update( int *additions, int *removals, int *m_additions, int *m_removals, int *switches ){
  
   /* Try flipping monomers, helps at large mass */
   *switches += switch_monomers();

    /* Starting site */
   int x = mersenne()*VOLUME;
   int f = mersenne()*2;
   int mu = mersenne()*NDIRS;
   int x2 = neighbour[mu][x];
   int flavorlist[N_FLAVOR] = {0,0};
   flavorlist[f] = 1;
   double notdone = 0;

   /* During the worm measure the number of configurations with an added source,
    * Provides a measurement of the bilinear expectation value */
   int nsteps = 0;

   /* Try starting a worm at a random site */
   if( (occupation_field[f][x] == 0) && (occupation_field[f][x2] == 0) ){
     /* Start by creating two monomers
      * one is the head and the other is a mass monomer 
      * Note the missing weight m, added at the end
      */
     double p = det_add_monomers( x, x2, flavorlist );
     if( mersenne() < p ){
        occupation_field[f][x] = 1;
        occupation_field[f][x2] = 1;
        mass_monomer[f][x2] = 1;
        update_current_determinant(flavorlist);
        notdone = 1;
        n_mass_monomer[f]+=1;
        n_occupied[f]+=2;
        *m_additions += 1;
     }
   } else if(mass_monomer[f][x] == 1) {
     /* Start by turning a mass monomer into  the head
      */
     if( mersenne() < 1./(m*m) ){
        mass_monomer[f][x] = 0;
        notdone = 1;
        n_mass_monomer[f]-=1;
        *m_removals += 1;
     }
  }

  while(notdone) {
    int step = mersenne()*5;
    int mu = mersenne()*NDIRS;
    int newx = neighbour[mu][x];

    nsteps++;

    switch(step){
    case 0:
      if( mersenne() < m*m ){
        /* End worm here, make endpoint a mass monomer */
        mass_monomer[f][x] = 1;
        n_mass_monomer[f] += 1;
        notdone = 0;
        *m_additions += 1;
      }
      break;

    case 1:
      if( mass_monomer[f][newx] == 1 ) {
        /* Worm annihilates with a neighbouring mass monomer */
        double p = det_remove_monomers( x, newx, flavorlist );
        if( mersenne() < p ){
          occupation_field[f][newx] = 0;
          occupation_field[f][x] = 0;
          mass_monomer[f][newx] = 0;
          update_current_determinant(flavorlist);
          notdone = 0;
          n_mass_monomer[f]-=1;
          n_occupied[f]-=2;
          *m_removals += 1;
        }
      }
      break;

    case 2:
      if( fourfermion_monomer[newx] == 1 ) {
        /* Remove the four fermion monomer,
         * head of worm moves and changes flavor
         */
        double p = det_remove_monomers( x, newx, flavorlist );
        if( mersenne() < p/U ){
           occupation_field[f][newx] = 0;
           occupation_field[f][x] = 0;
           fourfermion_monomer[newx] = 0;
           x = newx;
           update_current_determinant(flavorlist);
           n_occupied[f]-=2;
           flavorlist[f] = 0;
           f = !f;
           flavorlist[f] = 1;
           n_fourfermion_monomer-=1;
           *removals += 1;
        }
      }
      break;

    case 3:
      if( occupation_field[!f][x] == 0 && occupation_field[!f][newx] == 0   ) {
        /* Try to create a four fermion monomer */
        int f2 = !f; int flist2[2] = {0,0};
        flist2[f2]=1;
        double p = det_add_monomers( x, newx, flist2 );
        if( mersenne() < p*U ){
           occupation_field[f2][x] = 1;
           occupation_field[f2][newx] = 1;
           fourfermion_monomer[x] = 1;
           x = newx; 
           update_current_determinant(flist2);
           n_occupied[f2]+=2;
           flavorlist[f] = 0;
           f = !f;
           flavorlist[f] = 1;
           n_fourfermion_monomer+=1;
           *additions += 1;
          }
        } 
       break;

    case 4:
      mu = mersenne()*NDIRS;
      newx = neighbour[mu][newx];
      if( occupation_field[f][newx] == 0   ) {
        /* Just move the head of the worm
         * works even at U=0 */
        double p = det_move_monomers( x, newx, flavorlist );
        if( mersenne() < p ){
           occupation_field[f][x] = 0;
           occupation_field[f][newx] = 1;
           x = newx; 
           update_current_determinant(flavorlist);
          } 
        } 
       break;
    } //switch 
  } //while notdone
  
  /* Fix the weight
   * There are 5 possible steps, count each as 1/5 */
  return( (double)nsteps*m/5 ); 
}

#endif


