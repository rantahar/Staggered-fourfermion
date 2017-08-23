/****************************************************************
 * Local update 
 ****************************************************************/

#ifdef DEBUG
#include <fenv.h>
#endif

#include "Staggered.h"

#ifdef LOCAL_UPDATE

extern double U;
extern double m;
extern double linkmass;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
extern int n_occupied[N_FLAVOR];
extern int **occupation_field;

extern int n_fourfermion_monomer;
extern int n_mass_monomer[N_FLAVOR];
extern int *fourfermion_monomer;
extern int **mass_monomer;

/* Neighbour index arrays
 */
extern int **neighbour;



/* Count free sites that can be occupied
 */
int addable_sites(){
  int free_sites=0;
  for(int i=0; i<VOLUME ; i++){
    for(int nu=0; nu<ND; nu++){
      int j = neighbour[nu][i];
      if ( (fourfermion_monomer[i] == 0) && (fourfermion_monomer[j] == 0) &&
           (occupation_field[0][i] == 0) && (occupation_field[0][j] == 0) &&
           (occupation_field[1][i] == 0) && (occupation_field[1][j] == 0) ){
        free_sites++;
      }
    }
  }
  return free_sites;
}

/* Count sites where mass monomers can be added */
int addable_sites_mass( int f ){
  int free_sites=0;
  for(int i=0; i<VOLUME ; i++){  
    for(int nu=0; nu<ND; nu++){
      int j = neighbour[nu][i];
      if ( (mass_monomer[f][i] == 0) && (mass_monomer[f][j] == 0) &&
           (occupation_field[f][i] == 0) && (occupation_field[f][j] == 0) ){
        free_sites++;
      }
    }
  }
  return free_sites;
}

/* The number of addable sites after removing from two sites */
int addable_sites_after_removing( int x1, int x2){
  fourfermion_monomer[x1] = 0; fourfermion_monomer[x2] = 0;
  occupation_field[0][x1] = 0; occupation_field[0][x2] = 0;
  occupation_field[1][x1] = 0; occupation_field[1][x2] = 0;
  int free_sites=addable_sites();
  fourfermion_monomer[x1] = 1; fourfermion_monomer[x2] = 1;
  occupation_field[0][x1] = 1; occupation_field[0][x2] = 1;
  occupation_field[1][x1] = 1; occupation_field[1][x2] = 1;
  return free_sites;
}

/* The number of addable sites after removing from two sites */
int addable_sites_mass_after_removing( int x1, int x2, int f ){
  mass_monomer[f][x1] = 0; mass_monomer[f][x2] = 0;
  occupation_field[f][x1] = 0; occupation_field[f][x2] = 0;
  int free_sites=addable_sites_mass( f );
  mass_monomer[f][x1] = 1; mass_monomer[f][x2] = 1;
  occupation_field[f][x1] = 1; occupation_field[f][x2] = 1;
  return free_sites;
}


/* The number of sites where monomers can be removed */
int removable_sites(){
  int occupied_sites=0;
  for(int i=0; i<VOLUME ; i++){
    for(int nu=0; nu<ND; nu++){
      int j = neighbour[nu][i];
      if ( (fourfermion_monomer[i] == 1) && (fourfermion_monomer[j] == 1) &&
           (occupation_field[0][i] == 1) && (occupation_field[0][j] == 1) &&
           (occupation_field[1][i] == 1) && (occupation_field[1][j] == 1) ){
        occupied_sites++;
      }
    }
  }
  return occupied_sites;
}

/* The number of sites where mass monomers can be removed */
int removable_sites_mass( int f ){
  int occupied_sites=0;
  for(int i=0; i<VOLUME ; i++){
    for(int nu=0; nu<ND; nu++){
      int j = neighbour[nu][i];
      if ( (mass_monomer[f][i] == 1) && (mass_monomer[f][j] == 1) &&
           (occupation_field[f][i] == 1) && (occupation_field[f][j] == 1) ){
        occupied_sites++;
      }
    }
  }
  return occupied_sites;
}

/* The number of removable sites after adding two sites */
int removable_sites_after_adding( int x1, int x2){
  fourfermion_monomer[x1] = 1; fourfermion_monomer[x2] = 1;
  occupation_field[0][x1] = 1; occupation_field[0][x2] = 1;
  occupation_field[1][x1] = 1; occupation_field[1][x2] = 1;
  int occupied_sites=removable_sites();
  fourfermion_monomer[x1] = 0; fourfermion_monomer[x2] = 0;
  occupation_field[0][x1] = 0; occupation_field[0][x2] = 0;
  occupation_field[1][x1] = 0; occupation_field[1][x2] = 0;
  return occupied_sites;
}

/* The number of removable sites after adding two sites */
int removable_sites_mass_after_adding( int x1, int x2, int f){
  mass_monomer[f][x1] = 1; mass_monomer[f][x2] = 1;
  occupation_field[f][x1] = 1; occupation_field[f][x2] = 1;
  int occupied_sites=removable_sites_mass( f );
  mass_monomer[f][x1] = 0; mass_monomer[f][x2] = 0;
  occupation_field[f][x1] = 0; occupation_field[f][x2] = 0;
  return occupied_sites;
}







/* Turn a monomer on at a link */
void occupy_site(int x1, int x2, int m);

/* Turn a link off */
void free_site(int x1, int x2, int m);











#ifndef MASS_IN_MATRIX
/* Suggest adding monomers
 */
int add_fourfermion_monomer()
{
  int success = 0, x1,x2;
  int flavorlist[N_FLAVOR] = {1,1};
  int addable = addable_sites();
  if( addable > 0 ){

  do { 
   x1 = (int) (mersenne()*VOLUME);
   int nu = (int) (mersenne()*ND);
   x2 = neighbour[nu][x1];
  } while ( (fourfermion_monomer[x1] == 1) || (fourfermion_monomer[x2] == 1) ||
            (occupation_field[0][x1] == 1) || (occupation_field[0][x2] == 1) ||
            (occupation_field[1][x1] == 1) || (occupation_field[1][x2] == 1) ) ;

    double d = det_add_monomers( x1, x2, flavorlist );
    double n_factor = addable / (double)removable_sites_after_adding(x1,x2);
    double p = U*U*d * n_factor ;
    if( mersenne() < p ) {
      occupy_site( x1, x2, 0 );
      occupy_site( x1, x2, 1 );
      n_occupied[0]+=2;
      n_occupied[1]+=2;
      fourfermion_monomer[x1]=1;
      fourfermion_monomer[x2]=1;
      n_fourfermion_monomer+=2;
      success = 1;
      update_current_determinant(flavorlist);
    }
  

  }
  return success;
}

/* Suggest adding monomers
 */
int add_mass_monomer()
{
  int success = 0;
  int x1, x2;
  int f = (int) (mersenne()*N_FLAVOR);
  int flavorlist[N_FLAVOR] = {0,0};
  flavorlist[f] = 1;

  int addable = addable_sites_mass(f);
  if( addable > 0 ){

  do {
   x1 = (int) (mersenne()*VOLUME);
   int nu = (int) (mersenne()*ND);
   x2 = neighbour[nu][x1];
  } while ( (mass_monomer[f][x1] == 1)     || (mass_monomer[f][x2] == 1) ||
            (occupation_field[f][x1] == 1) || (occupation_field[f][x2] == 1) ) ;

    double d = det_add_monomers( x1, x2, flavorlist );
    double p = m*m*d * addable / (double)removable_sites_mass_after_adding(x1,x2,f);
    if( mersenne() < p ) {
      occupy_site( x1, x2, f );
      n_occupied[f]+=2;
      mass_monomer[f][x1] = 1;
      mass_monomer[f][x2] = 1;
      n_mass_monomer[f] += 2;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }

  return success;
}

/* Suggest removing monomers
 */
int remove_fourfermion_monomer()
{
  int success = 0;
  int x1, x2;
  int flavorlist[N_FLAVOR] = {1,1};
  int removable = removable_sites();
  if( removable > 0 ){
  
  do { /* Site needs to be fully occupied */
   x1 = (int) (mersenne()*VOLUME);
   int nu = (int) (mersenne()*ND);
   x2 = neighbour[nu][x1];
  } while ( (fourfermion_monomer[x1] == 0) || (fourfermion_monomer[x2] == 0) ||
            (occupation_field[0][x1] == 0) || (occupation_field[0][x2] == 0) ||
            (occupation_field[1][x1] == 0) || (occupation_field[1][x2] == 0) ) ;

    double d = det_remove_monomers( x1, x2, flavorlist  );
    double n_factor = removable / (double)addable_sites_after_removing(x1,x2);
    double p = d/(U*U) * n_factor ;
    if( mersenne() < p ){
      free_site( x1, x2, 0 );
      free_site( x1, x2, 1 );
      n_occupied[0]-=2;
      n_occupied[1]-=2;
      fourfermion_monomer[x1]=0;
      fourfermion_monomer[x2]=0;
      n_fourfermion_monomer-=2;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }
  return success;
}

/* Suggest removing monomers
 */
int remove_mass_monomer()
{
  int success = 0;
  int x1, x2;
  int f = (int) (mersenne()*N_FLAVOR);
  int flavorlist[N_FLAVOR] = {0,0};
  flavorlist[f] = 1;

  int removable = removable_sites_mass(f);
  if( removable > 0 ){

  do { /* Fermion f at the sites needs to be occupied */
   x1 = (int) (mersenne()*VOLUME);
   int nu = (int) (mersenne()*ND);
   x2 = neighbour[nu][x1];
  } while ( (mass_monomer[f][x1] == 0)     || (mass_monomer[f][x2] == 0) ||
            (occupation_field[f][x1] == 0) || (occupation_field[f][x2] == 0) ) ;
    double d = det_remove_monomers( x1, x2, flavorlist );
    double p = d/(m*m) * (double)removable / (double)addable_sites_mass_after_removing(x1,x2,f);
    if( mersenne() < p ) {
      free_site( x1, x2, f );
      n_occupied[f]-=2;
      mass_monomer[f][x1] = 0;
      mass_monomer[f][x2] = 0;
      n_mass_monomer[f] -= 2;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }

  return success;
}


#define N_MOVE_ATTEMPS 1

/* Suggest moving monomers
 */
int move_fourfermion_monomer()
{
  int moves = 0;
  for( int attempt=0; attempt < N_MOVE_ATTEMPS; attempt++) {
    int x1, x2;
    int flavorlist[N_FLAVOR] = {1,1};

    x1 = (int) (mersenne()*VOLUME);
    int nu = mersenne()*NDIRS;
    x2 = neighbour[nu][x1];
    nu = mersenne()*NDIRS;
    x2 = neighbour[nu][x2];

    if ( (fourfermion_monomer[x1] == 1) && (fourfermion_monomer[x2] == 0) &&
         (occupation_field[0][x1] == 1) && (occupation_field[0][x2] == 0) &&
         (occupation_field[1][x1] == 1) && (occupation_field[1][x2] == 0) ) 
    {
      double p = det_move_monomers( x1, x2, flavorlist );
      if( mersenne() < p ) {
        occupation_field[0][x1]=0;
        occupation_field[1][x1]=0;
        occupation_field[0][x2]=OCCUPIED;
        occupation_field[1][x2]=OCCUPIED;
        fourfermion_monomer[x1] = 0;
        fourfermion_monomer[x2] = 1;
        moves++;
        update_current_determinant(flavorlist);
      }
    }
  }
  return moves;
}

/* Suggest moving monomers
 */
int move_mass_monomer()
{
  int moves = 0;
  int f = (int) (mersenne()*N_FLAVOR);
  for( int attempt=0; attempt < N_MOVE_ATTEMPS; attempt++) {
    int x1, x2;
    int flavorlist[N_FLAVOR] = {0,0};
    flavorlist[f] = 1;

    x1 = (int) (mersenne()*VOLUME);
    int nu = mersenne()*NDIRS;
    x2 = neighbour[nu][x1];
    nu = mersenne()*NDIRS;
    x2 = neighbour[nu][x2];

    if ( (mass_monomer[f][x1] == 1) && (mass_monomer[f][x2] == 0) &&
         (occupation_field[f][x1] == 1) && (occupation_field[f][x2] == 0) )
    {
      double p = det_move_monomers( x1, x2, flavorlist );
      if( mersenne() < p ) {
        occupation_field[f][x1]=0;
        occupation_field[f][x2]=OCCUPIED;
        mass_monomer[f][x1] = 0;
        mass_monomer[f][x2] = 1;
        moves++;
        update_current_determinant(flavorlist);
      }
    }
    
  }
  return moves;
}

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


/* Do a number of random update attempts  */
#define N_updates 6
void local_update( int *additions, int *removals, int *moves,int *m_additions, int *m_removals, int *m_moves, int *switches )
{
  int changes=0;

  for( int i=0; i<N_updates; i++){
    int step = mersenne()*6; 
    switch(step){
      case 0:
        *additions += add_fourfermion_monomer();
        break;
      case 1:
        *removals += remove_fourfermion_monomer();
        break;
      case 2:
        *moves += move_fourfermion_monomer();
        break;
      case 3:
        *m_moves += move_mass_monomer();
        break;
      case 4:
        *m_additions += add_mass_monomer();
        break;
      case 5:
        *m_removals += remove_mass_monomer();
        break;
      default:
        *switches += switch_monomers();
       break;
    }
  }
}




















#else //MASS_IN_MATRIX

int add_fourfermion_monomer()
{
  int success = 0, x1,x2;
  int flavorlist[N_FLAVOR] = {1,1};

  x1 = (int) (mersenne()*VOLUME);

  if( (fourfermion_monomer[x1] == 0) && (occupation_field[0][x1] == 0) &&
      (occupation_field[1][x1] == 0) ){
    double d = det_add_monomers1( x1, flavorlist );
    double p = U*d;
    if( mersenne() < p ) {
      occupation_field[0][x1] = 1;
      occupation_field[1][x1] = 1;
      n_occupied[0]+=1;
      n_occupied[1]+=1;
      fourfermion_monomer[x1]=1;
      n_fourfermion_monomer+=1;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }
  return success;
}

int remove_fourfermion_monomer()
{
  int success = 0;
  int x1, x2;
  int flavorlist[N_FLAVOR] = {1,1};

  x1 = (int) (mersenne()*VOLUME);

  if( (fourfermion_monomer[x1] == 1) && (occupation_field[0][x1] == 1) &&
      (occupation_field[1][x1] == 1) ) {
    double d = det_remove_monomers1( x1, flavorlist  );
    double p = d/U;
    if( mersenne() < p ){
      occupation_field[0][x1] = 0;
      occupation_field[1][x1] = 0;
      n_occupied[0]-=1;
      n_occupied[1]-=1;
      fourfermion_monomer[x1]=0;
      n_fourfermion_monomer-=1;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }
  return success;
}


/* Suggest adding two monomers
 */
int add_fourfermion_monomer2()
{
  int success = 0, x1,x2;
  int flavorlist[N_FLAVOR] = {1,1};
 
  x1 = (int) (mersenne()*VOLUME);
  int nu = (int) (mersenne()*ND);
  x2 = neighbour[nu][x1];

  if( (fourfermion_monomer[x1] == 0) && (occupation_field[0][x1] == 0) &&
      (occupation_field[1][x1] == 0) && (fourfermion_monomer[x2] == 0) && 
      (occupation_field[0][x2] == 0) && (occupation_field[1][x2] == 0) ) {

    double d = det_add_monomers( x1, x2, flavorlist );
    double p = U*U*d ;
    if( mersenne() < p ) {
      occupy_site( x1, x2, 0 );
      occupy_site( x1, x2, 1 );
      n_occupied[0]+=2;
      n_occupied[1]+=2;
      fourfermion_monomer[x1]=1;
      fourfermion_monomer[x2]=1;
      n_fourfermion_monomer+=2;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }
  return success;
}


/* Suggest removing two monomers
 */
int remove_fourfermion_monomer2()
{
  int success = 0;
  int x1, x2;
  int flavorlist[N_FLAVOR] = {1,1};
  
  x1 = (int) (mersenne()*VOLUME);
  int nu = (int) (mersenne()*ND);
  x2 = neighbour[nu][x1];
  if( (fourfermion_monomer[x1] == 1) && (occupation_field[0][x1] == 1) &&
      (occupation_field[1][x1] == 1) && (fourfermion_monomer[x2] == 1) && 
      (occupation_field[0][x2] == 1) && (occupation_field[1][x2] == 1) ) {

    double d = det_remove_monomers( x1, x2, flavorlist  );
    double p = d/(U*U) ;
    if( mersenne() < p ){
      free_site( x1, x2, 0 );
      free_site( x1, x2, 1 );
      n_occupied[0]-=2;
      n_occupied[1]-=2;
      fourfermion_monomer[x1]=0;
      fourfermion_monomer[x2]=0;
      n_fourfermion_monomer-=2;
      success = 1;
      update_current_determinant(flavorlist);
    }
  }
  return success;
}




/* Do a number of random update attempts  */
#define N_updates 4
void local_update( int *additions, int *removals, int *moves,int *m_additions, int *m_removals, int *m_moves, int *switches )
{
  int changes=0;

  for( int i=0; i<N_updates; i++){
    int step = mersenne()*4; 
    switch(step){
      case 0:
        *additions += add_fourfermion_monomer();
        break;
      case 1:
        *removals += remove_fourfermion_monomer();
        break;
      case 2:
        *additions += 2*add_fourfermion_monomer2();
        break;
      case 3:
        *removals += 2*remove_fourfermion_monomer2();
        break;
    }
  }
}









#endif





#endif //LOCAL_UPDATE
