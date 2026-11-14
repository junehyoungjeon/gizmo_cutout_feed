#ifdef FOFRAD

/* variables for fofrad.c code */
/*
 */


#include "../allvars.h"

#define FOF_MINGROUPLEN 8	// minimum length of fofrad group
#define FOFRAD_HALO_RHOCRIT 500	// overdensity relative to mean for identifying halos
#define FOFRAD_MAXMASS 2e12     // maximum mass of single fofrad group

#define NDIM 3

typedef struct fof_part_data
{
  double Pos[NDIM];
  double Mass;
#ifdef SIGMASFR_WIND
  double Sfr;
  double r2;
#endif
  short int Type;
#ifdef OUTPUT_FOFRAD_PARTICLES
  MyIDType ID;
#endif

  unsigned int localIndex;
  unsigned int Task;
  unsigned int galID;
}
fof_part_data;

typedef struct fof_group_data
{
  //int jlo,jhi;
  double xlo[NDIM],xhi[NDIM];
  double Mbar, Mgas, Mstar;
#ifdef SIGMASFR_WIND
  double Sfr,rw2;
#endif
  double COM[NDIM];
}
fof_group_data;

int fof_compare_galID(),fof_compare_xdir(),fof_compare_ydir(),fof_compare_zdir(),fof_compare_ID();
double fofrad_periodic();
#ifdef SIGMASFR_WIND
int fof_compare_r2();
#endif


#endif
