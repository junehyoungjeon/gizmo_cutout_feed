/*
  TODO:
  * add to cpu.txt

  - what about groups that are dropped?  wasting memory!
  * forget about em!
*/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "../allvars.h"
#include "../proto.h"
#include "fofrad.h"

#ifdef FOFRAD

/*! \file fofrad.c
 *  Fast serial approximate FOF group finder, by Romeel Dave, last updated 2/19/2015
 *
 *  The algorithm uses quicksort to group particles that are within +/-LinkL of each other in
 *  each direction (x,y,z) independently.  This is fast, but it's not exactly FOF, since
 *  instead of using a spherical volume of radius LinkL for the search, it uses a cube
 *  of length 2*LinkL on a side.  
 *
 *  LinkL is in code distance units, *not* in units of mean interparticle spacing.
 */

#define nHthresh 0.13  // H atoms/cm^3
#define overthresh 1000  // overdensity threshold for including gas in galaxy (mostly for hi-z)

fof_group_data *Groups;
int NPfof,Ngroup,NumPart;

/* 
   links particles into FOF groups.  nHthresh in H atoms/cm^3 
   
   call with a positive number to output groups to a file
*/
int fof_rad(int snapnum)
{
  int i,j,idim,jprev,Nsub,Nstar;
  int NPtot;
  double xprev,Mgas,Mstar,Mdark,COM[NDIM];
  double densconv,massconv;
#ifdef SIGMASFR_WIND
  double Sfr,sfrwind;
#endif
  double celldistance(),rhovir_over_rhobar(),NFW_massfactor();
  //double localMinParticleMass  = 1.0e8;
  //double globalMinParticleMass = 1.0e8;

  double fofrad_periodic();
  fof_part_data *Pfof, *All_Pfof;

  double tbeg = 0.0;
  double tend = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  tbeg = my_second();

  /* calculate Linking-Length */
  int ndm = 0;
  long long ndmtot = 0;
  double LinkL, rhodm, massDM, massDMtot;
  massDM = 0.0;
  massDMtot = 0.0;
  for(i=0, ndm=0; i<NumPart; i++){
    if(P[i].Type == 1){
      ndm += 1;
      massDM += P[i].Mass;
    }
  }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&massDM, &massDMtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
  LinkL = FOFRAD * pow(massDMtot / ndmtot / rhodm, 1.0 / 3);
  /* ------------------------ */
  
  densconv = All.HubbleParam*All.HubbleParam * 
    All.cf_a3inv*HYDROGEN_MASSFRAC*All.UnitMass_in_g /
    (All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm*PROTONMASS);	// converts to H number density
  massconv = All.UnitMass_in_g/(1.989e33*All.HubbleParam); // converts to solar masses
  double nHcrit = 1.88e-29 * All.HubbleParam * All.HubbleParam * All.cf_a3inv * All.OmegaBaryon * HYDROGEN_MASSFRAC / PROTONMASS; // physical nH at critical density 
  if( ThisTask==0 ) printf("\nFOFRAD: z=%g LL(comoving kpc/h)=%g nHthresh=%g overdensthresh=%g\n",1./All.cf_atime-1,LinkL,nHthresh,overthresh*nHcrit);

  t0 = my_second();
  for( i=0,NPfof=0; i<NumPart; i++){
    P[i].galID = 0;

    /*get minimum particle mass 
    if( P[i].Type==0 || P[i].Type==4 )
      if( P[i].Mass < localMinParticleMass )
	localMinParticleMass = P[i].Mass;*/
	
    if( (P[i].Type==0 && SphP[i].Density*densconv > nHthresh ) || P[i].Type==4 )
      NPfof += 1;
  }
  Pfof = (fof_part_data *) mymalloc("FOFRAD_particles",NPfof*sizeof(fof_part_data));
  for( i=0,NPfof=0; i<NumPart; i++ ) { 		// select particles to group
    if( (P[i].Type==0 && SphP[i].Density*densconv > nHthresh && SphP[i].Density*densconv > overthresh*nHcrit ) || P[i].Type==4 ) {

      for( j=0; j<NDIM; j++ ) Pfof[NPfof].Pos[j] = P[i].Pos[j];
      Pfof[NPfof].localIndex = i;
      Pfof[NPfof].Mass  = P[i].Mass;
      Pfof[NPfof].Type  = P[i].Type;
#ifdef SIGMASFR_WIND
      Pfof[NPfof].Sfr  = SphP[i].Sfr;
#endif
      Pfof[NPfof].Task  = ThisTask;
      Pfof[NPfof].galID = 0;
#ifdef OUTPUT_FOFRAD_PARTICLES
      Pfof[NPfof].ID    = P[i].ID;
#endif      
      NPfof++;
    }
  }
  MPI_Allreduce(&NPfof, &NPtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //tot # of particles to group
/*  MPI_Allreduce(&localMinParticleMass, &globalMinParticleMass, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); //get min particle mass among all gas & stars
  globalMinParticleMass *= massconv; //to Msun */

  //if there are no particles to group, exit
  if(NPtot < FOF_MINGROUPLEN){
    if(ThisTask==0) printf("FOFRAD: not enough particles to group (%d)! (z=%0.2f)\n",NPtot,1./All.Time-1.);

    myfree(Pfof);

    //set Groups array to hold one default group
    if( Groups != NULL ) free(Groups);
    Groups = (fof_group_data *)malloc(sizeof(fof_group_data));
    Groups[0].Mstar = 0.0;
    Groups[0].Mgas  = 0.0;
    Groups[0].Mbar  = All.OrigGasMass * FOF_MINGROUPLEN * massconv;

    return 0;
  }
 
  // gather particle counts from each processor
  int local_bytes = NPfof * sizeof(fof_part_data);
  int all_bytes[NTask];
  
  // get bytes per processor
  MPI_Allgather(&local_bytes,1,MPI_INT,all_bytes,1,MPI_INT,MPI_COMM_WORLD);

  // setup offsets
  int tally = 0;
  int offsets[NTask];
  offsets[0] = 0;
  for(i=1; i<NTask; i++){
    tally += all_bytes[i-1];
    offsets[i] = tally;
  }

  // allocate and get FULL particle struct
  All_Pfof = (fof_part_data *) mymalloc("ALL_FOFRAD_Particles",NPtot*sizeof(fof_part_data));
  MPI_Allgatherv(Pfof,local_bytes,MPI_BYTE,All_Pfof,all_bytes,offsets,MPI_BYTE,MPI_COMM_WORLD);

  t1 = my_second(); 
  if(ThisTask==0) printf("FOFRAD: filter and exchange for %d particles took %g sec\n",NPtot,timediff(t0,t1));

  t0 = my_second();
  /* Sort in x-dir, break up into groups */
  qsort(All_Pfof, NPtot, sizeof(fof_part_data), fof_compare_xdir);
  All_Pfof[0].galID=1;
  for( i=1,Ngroup=1,Nsub=0,Mdark=0,xprev=All_Pfof[0].Pos[0]; i<NPtot; i++ ) {	// loop over particles
    if( All_Pfof[i].Pos[0]-xprev > LinkL) {	// if delta-x > LinkL, start new group
      Ngroup++;
      All_Pfof[i].galID=Ngroup;
      Nsub = 0;
      Mdark = 0;
    }
    else All_Pfof[i].galID=All_Pfof[i-1].galID;	// else assign to previous particle's group ID
    Mdark += All_Pfof[i].Mass*massconv;
    Nsub++;
    xprev = All_Pfof[i].Pos[0];
  }

  /* Within each x-dir group, sort in y, then in z, each time breaking up further into groups */
  for( idim=1; idim<NDIM; idim++ ) {	// loop over y,z directions
  for( i=1,jprev=0; i<=Ngroup; i++ ) {	// loop over groups identified in x-dir
    Nsub=0; while( All_Pfof[Nsub+jprev].galID == i ) Nsub++; // count group i particles
    if( Nsub >= FOF_MINGROUPLEN ) {	// if group is already too small, ignore
      if( idim==1 ) qsort(All_Pfof+jprev, Nsub, sizeof(fof_part_data), fof_compare_ydir); 
      if( idim==2 ) qsort(All_Pfof+jprev, Nsub, sizeof(fof_part_data), fof_compare_zdir); 
      for( j=jprev+1,xprev=All_Pfof[jprev].Pos[idim]; j<jprev+Nsub; j++ ) {  // split into new groups
        if( All_Pfof[j].Pos[idim]-xprev > LinkL ) {	// if delta-x > LinkL, start new group
          Ngroup++;
          All_Pfof[j].galID=Ngroup;
        }
        else All_Pfof[j].galID=All_Pfof[j-1].galID;	// else assign to previous particle's group ID
        xprev = All_Pfof[j].Pos[idim];
      }
    }
    jprev += Nsub;
  }
  qsort(All_Pfof, NPtot, sizeof(fof_part_data), fof_compare_galID);  // sort by group ID
  }
  t1 = my_second();
  if(ThisTask == 0) printf("FOFRAD: found %d groups, took %g sec\n",Ngroup,timediff(t0,t1));

  if( Groups != NULL ) free(Groups);
  Groups = (fof_group_data *)malloc((Ngroup+1)*sizeof(fof_group_data));

  //if no groups, exit
  if(Ngroup == 0){
    myfree(All_Pfof);
    myfree(Pfof);
    Groups[0].Mstar = 0.0;
    Groups[0].Mgas  = 0.0;
    Groups[0].Mbar  = All.OrigGasMass * FOF_MINGROUPLEN * massconv;
#ifdef SIGMASFR_WIND
    Groups[0].rw2 = 0.0;
    Groups[0].Sfr = 0.0;
#endif
    return 0;
  }

  //open group file
  FILE *FdGroupFile=NULL; 
  if(ThisTask == 0 && snapnum >= 0){
    char buf[200];
    sprintf(buf, "%sFOFRADGROUPS", All.OutputDir);
    mkdir(buf, 02755);
    sprintf(buf, "%sFOFRADGROUPS/groups_%03d.txt", All.OutputDir, snapnum);
    if(!(FdGroupFile = fopen(buf, "w"))){
      printf("error in opening file '%s'\n",buf);
      endrun(1);
    }
    fprintf(FdGroupFile,"#a\tgalID\tNpart\tNstar\tMbar\tMgas\tMstar\tcmx\tcmy\tcmz\n");
    fflush(FdGroupFile);

#ifdef OUTPUT_FOFRAD_PARTICLES
    FILE *FdParticleFile;
    sprintf(buf, "%sFOFRADGROUPS/particles_%03d.txt", All.OutputDir, snapnum);
    if(!(FdParticleFile = fopen(buf, "w"))){
      printf("error in opening file '%s'\n",buf);
      endrun(1);
    }
    fprintf(FdParticleFile,"#galID\tType\tPID\n");
    fflush(FdParticleFile);
#endif
  }

  t0 = my_second();
  /* output group counts */
  //fprintf(stdout,"# ID Npart Nstar  Mtot Mbary Mgas Mstar  xcm ycm zcm\n");
  for( i=1,jprev=0; i<=Ngroup; i++ ) {
    Mstar=Mgas=Mdark=Nsub=Nstar=0; 
    for( j=0; j<NDIM; j++ ) COM[j] = 0;
    while( All_Pfof[Nsub+jprev].galID == i ) {

      switch( All_Pfof[Nsub+jprev].Type ){
      case 0:
	Mgas += All_Pfof[Nsub+jprev].Mass;
	break;
      case 4:
	Mstar += All_Pfof[Nsub+jprev].Mass;
	Nstar++;
	break;
      default:
	printf("FOFRAD: Error! Detected particle type %d on task %d\n",All_Pfof[Nsub+jprev].Type,ThisTask);
	//endrun(6522);
	break;
      }

      for( j=0; j<NDIM; j++ ) COM[j] += All_Pfof[Nsub+jprev].Mass*All_Pfof[Nsub+jprev].Pos[j];
      Nsub++;
    }
    for( j=0; j<NDIM; j++ ) COM[j] /= (Mgas+Mstar+Mdark);

#ifdef SIGMASFR_WIND
    /* Compute radius from COM for each particle */
    for( j=jprev; j<jprev+Nsub; j++ ) {
	All_Pfof[j].r2 = fofrad_periodic(All_Pfof[j].Pos[0]-COM[0])*fofrad_periodic(All_Pfof[j].Pos[0]-COM[0]) 
		+ fofrad_periodic(All_Pfof[j].Pos[1]-COM[1])*fofrad_periodic(All_Pfof[j].Pos[1]-COM[1]) 
		+ fofrad_periodic(All_Pfof[j].Pos[2]-COM[2])*fofrad_periodic(All_Pfof[j].Pos[2]-COM[2]);
    }
    qsort(All_Pfof+jprev, Nsub, sizeof(fof_part_data), fof_compare_r2);  // sort by radius
    /* compute "wind radius", i.e. radius beyond which Sigma_Sfr<SIGMASFR_WIND */
    Sfr=Groups[i].rw2=0;
    for( j=jprev; j<jprev+Nsub; j++ ) {
      if( All_Pfof[j].Type == 0 ) {
	Sfr += All_Pfof[j].Sfr;
	if( Sfr/(M_PI*All_Pfof[j].r2*All.cf_atime*All.cf_atime+1.e-30) > SIGMASFR_WIND && Groups[i].rw2 < All_Pfof[j].r2 ) Groups[i].rw2 = All_Pfof[j].r2;
	//if( ThisTask==0 && Nsub>12 ) printf("SIGMASFR: %d %g %d %d %d %g %g %g %g %g\n",ThisTask,All.cf_atime,i,j-jprev,Nsub,All_Pfof[j].Sfr,Sfr,Sfr/(M_PI*All_Pfof[j].r2*All.cf_atime*All.cf_atime+1.e-30),All_Pfof[j].r2,Groups[i].rw2);
      }
    }
    for( j=jprev,sfrwind=0; j<jprev+Nsub; j++ ) if( All_Pfof[j].Type == 0 && All_Pfof[j].r2 < Groups[i].rw2 ) sfrwind+= All_Pfof[j].Sfr;
    Groups[i].Sfr   = Sfr;
    //if( ThisTask==0 && Nsub-Nstar>0 ) printf("SIGMASFR: z= %g i= %d Np= %d Ng= %d rw= %g M*= %g SFR= %g frac= %g\n",1./All.cf_atime-1,i,Nsub,Nsub-Nstar,sqrt(Groups[i].rw2),log10(Mstar*massconv+1.e-30),Sfr,sfrwind/(Sfr+1.e-30));
#endif 

    /* Assign group properties */
    Mgas *= massconv; Mstar *= massconv; Mdark *= massconv;	// to solar masses
    Groups[i].Mgas   = Mgas;
    Groups[i].Mstar  = Mstar;
    Groups[i].Mbar   = Mgas+Mstar;
    Groups[i].COM[0] = COM[0];
    Groups[i].COM[1] = COM[1];
    Groups[i].COM[2] = COM[2];

    if( Nsub >= FOF_MINGROUPLEN && ThisTask==0 ){
      if(snapnum >= 0){
	fprintf(FdGroupFile,"%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n",
		All.Time,i,Nsub,Nstar,Mgas+Mstar,Mgas,Mstar,COM[0],COM[1],COM[2]);
	fflush(FdGroupFile);
      }
      /*
      else
	printf("FOFRADGRP z=%g %d %d %d %g %g %g %g %g %g\n",
	       1./All.Time - 1.,i,Nsub,Nstar,Mgas+Mstar,Mgas,Mstar,COM[0],COM[1],COM[2]);    
      */
    }
    
    jprev += Nsub;
  }

  /* Set small groups to galID=0 */
  int Nvalid = 0;
  Mstar=Mgas=Mdark=0; 
  for( i=1,jprev=0; i<=Ngroup; i++ ) {
    Nsub=0; 
    while( All_Pfof[Nsub+jprev].galID == i ) {

      switch( All_Pfof[Nsub+jprev].Type ){
      case 0:
	Mgas += All_Pfof[Nsub+jprev].Mass;
	break;
      case 4:
	Mstar += All_Pfof[Nsub+jprev].Mass;
	break;
      default:
	printf("FOFRAD: Error! Detected particle type %d on task %d\n",All_Pfof[Nsub+jprev].Type,ThisTask);
	//endrun(6523);
	break;
      }

      Nsub++;
    }
    if( Nsub < FOF_MINGROUPLEN )
      for( j=jprev; j<jprev+Nsub; j++ ) All_Pfof[j].galID = 0;
    else{
#ifdef OUTPUT_FOFRAD_PARTICLES
      if(ThisTask == 0 && snapnum >= 0){
	for( j=jprev; j<jprev+Nsub; j++ )
	  if( All_Pfof[j].galID == i ){
#ifdef LONGIDS
	    fprintf(FdParticleFile,"%d %d %llu\n",i,All_Pfof[j].Type,All_Pfof[j].ID);
#else
	    fprintf(FdParticleFile,"%d %d %u\n",  i,All_Pfof[j].Type,All_Pfof[j].ID);
#endif
	  }
	fflush(FdParticleFile);
      }
#endif
      Nvalid += 1;
    }
    jprev += Nsub;
  }
  Mgas *= massconv; Mstar *= massconv; Mdark *= massconv;
  qsort(All_Pfof, NPtot, sizeof(fof_part_data), fof_compare_galID);  // sort by group ID

  Groups[0].Mstar = 0.0;
  Groups[0].Mgas  = 0.0;
  Groups[0].Mbar  = All.OrigGasMass * FOF_MINGROUPLEN * massconv;

  t1 = my_second();
  //if(ThisTask==0) printf("FOFRAD: data tally and group filter took %g sec\n",timediff(t0,t1));

  tend = my_second();
  if( ThisTask==0 ) {
    printf("FOFRAD: Found %d valid groups from %d particles, total Mgas= %g Mstar= %g tcpu= %g\n",
	   Nvalid,NPtot,Mgas,Mstar,timediff(tbeg,tend));
    
    if(snapnum >= 0){
      fclose(FdGroupFile);
#ifdef OUTPUT_FOFRAD_PARTICLES
      fclose(FdParticleFile);
#endif
    }
  }

  /* insert group mass info into particle list */
  for( i=0; i<NPtot; i++ )
    if(All_Pfof[i].Task == ThisTask)
      P[All_Pfof[i].localIndex].galID = All_Pfof[i].galID;
  
  //print out particle list
  /*
  if(ThisTask == 0){
    for( i=0; i<NPtot; i++ )
      if(All_Pfof[i].galID > 0)
	printf("FOFRADPART z=%g %u %d %d %g %g %g\n", 1./All.Time - 1.,P[All_Pfof[i].localIndex].ID,
	       All_Pfof[i].Type,All_Pfof[i].galID,All_Pfof[i].Pos[0],All_Pfof[i].Pos[1],All_Pfof[i].Pos[2]);
  }
  */

  myfree(All_Pfof);
  myfree(Pfof);

  return 0;
}

double celldistance(int i,int j)
{
  int idim;
  double dx[NDIM],r2;
  double fofrad_periodic();

  for( idim=0; idim<NDIM; idim++ ) {
    if( Groups[i].xlo[idim] > Groups[j].xhi[idim] ) dx[idim] = fofrad_periodic(Groups[i].xlo[idim] - Groups[j].xhi[idim]);
    else if( Groups[i].xhi[idim] < Groups[j].xlo[idim] ) dx[idim] = fofrad_periodic(Groups[j].xlo[idim] - Groups[i].xhi[idim]);
    else dx[idim] = 0.;
  }
  for( idim=0,r2=0.; idim<NDIM; idim++ ) r2 += dx[idim]*dx[idim];

  return r2;
}

double Omegaf(double Omega0, double Lambda0, double z)
{
    double answer,zplus12,zplus13;

    zplus12 = (1.0+z)*(1.0+z);
    zplus13 = zplus12*(1.0+z);
    answer = Omega0*zplus13/(Omega0*zplus13+(1.-Omega0-Lambda0)*zplus12+Lambda0);
    /*    fprintf(stderr,"Omegaf: %lg\n",answer); */
    return (answer);
}

double rhovir_over_rhobar(double Omega0, int Lambda_OPT, double z)
{
    double etaf,wf,Lambda0,answer;

    if (Omega0 == 1.0) {
        return(178.0);
    }
    if (Lambda_OPT) {
        Lambda0 = 1.0-Omega0;
        wf = 1./Omegaf(Omega0,Lambda0,z) - 1.;
        answer = 18.*(M_PI*M_PI)*(1.+0.4093*pow(wf,0.9052));
    } else {
        etaf = acosh(2.0/Omegaf(Omega0,0.0,z) - 1.0);
        answer =  4.0*(M_PI*M_PI)/(pow(sinh(etaf)-etaf,2));
        answer *= pow(cosh(etaf)-1,3);
    }
    if(ThisTask==0) printf("FOFRAD: overdensity at z=%g: %lg\n",1./All.Time-1,answer);
    return(answer);
}

float NFW_rho(float r,float conc)
{
        float rs;
        rs = 1./conc;
        return rs/(r*(1+r/rs)*(1+r/rs));
}

double NFW_massfactor()
{
	float conc=10;
	float delta=FOFRAD_HALO_RHOCRIT;
        float x,dx,omega0=0.3;
        float dm1,dm2,mvir,mass,rhovir,meanrho,virialdens;

        virialdens = rhovir_over_rhobar(omega0,1,1./All.Time-1);
        dx = 0.001;
        for(x=dx,mvir=0.;x<1.;x+=dx) {
            dm1 = 4*M_PI*x*x*NFW_rho(x,conc);
            dm2 = 4*M_PI*(x+dx)*(x+dx)*NFW_rho(x+dx,conc);
            mvir += 0.5*(dm1+dm2)*dx;
        }
        rhovir = mvir/(4./3*M_PI);
        fprintf(stderr,"%g %g\n",mvir,rhovir);
        for(x=dx,mass=0.;x<1.;x+=dx) {
            dm1 = 4*M_PI*x*x*NFW_rho(x,conc);
            dm2 = 4*M_PI*(x+dx)*(x+dx)*NFW_rho(x+dx,conc);
            mass += 0.5*(dm1+dm2)*dx;
            meanrho = mass/(4./3*M_PI*x*x*x);
            if( (meanrho/rhovir) < delta/virialdens ) break;
        }
        fprintf(stderr,"delta= %g  fmass= %g  fradius= %g\n",delta,mass/mvir,x);

    return mass/mvir;
}

#ifdef GALSF_SUBGRID_VBOOST
/* Boost wind velocity to account for potential between particle location and GALSF_SUBGRID_VBOOST*Rvir; for FIRE case, GALSF_SUBGRID_VBOOST=0.25. */
double wind_vel_boost(int i, double v_circ, double v_wind)
{
	double vw_new,v_c,s,g_c,c_halo,Phi_out,Phi_in,r200c;
        double HaloConcentrationNorm = 9.;  /* concentration c0 of a halo of unit mass */
        double HaloConcentrationSlope = -0.15;  /* slope n of mass concentration relation, namely c = c0 * M_200,crit^n */

	if( P[i].galID == 0 ) return v_wind; 	// don't change vel of ungrouped particles

	double Mhalo = Groups[P[i].galID].Mbar*All.Omega0/All.OmegaBaryon;
	if( Mhalo < FOF_MINGROUPLEN*All.OrigGasMass*All.Omega0/All.OmegaBaryon ) Mhalo = FOF_MINGROUPLEN*All.OrigGasMass*All.Omega0/All.OmegaBaryon;
	Mhalo *= All.HubbleParam*SOLAR_MASS/All.UnitMass_in_g;  // convert to system units

	double rhocrit = 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
        rhocrit *= All.Omega0/All.cf_a3inv + (1-All.Omega0-All.OmegaLambda)/All.cf_a2inv + All.OmegaLambda; /* physical critical density at redshift z */
	r200c  = pow(Mhalo / (4 * M_PI / 3.0 * 200 * rhocrit), 1.0 / 3.0); 
	s = (P[i].Pos[0]-Groups[P[i].galID].COM[0])*(P[i].Pos[0]-Groups[P[i].galID].COM[0])
	    + (P[i].Pos[1]-Groups[P[i].galID].COM[1])*(P[i].Pos[1]-Groups[P[i].galID].COM[1])
	    + (P[i].Pos[2]-Groups[P[i].galID].COM[2])*(P[i].Pos[2]-Groups[P[i].galID].COM[2]);
	//if( ThisTask==0 ) printf("GALSF_SUBGRID_VBOOST: %d %g %g %g %g %g %g %g\n",P[i].galID,Mhalo,rhocrit,r200c,sqrt(s),P[i].Pos[0]-Groups[P[i].galID].COM[0],P[i].Pos[0],Groups[P[i].galID].COM[0]);
	s = s/(r200c*r200c);	// distance in units of r200
	if( s > GALSF_SUBGRID_VBOOST*GALSF_SUBGRID_VBOOST ) return v_wind;	// if already outside boost radius, don't change velocity

	/* OK, so we are now going to calculate the new velocity */
	s = sqrt(s);
        //v_c = sqrt(All.G * Mhalo / r200c);    /* physical circular velocity at r_200,crit */
        v_c = v_circ;
	c_halo = HaloConcentrationNorm * pow(Mhalo, HaloConcentrationSlope);
	g_c = 1./(log(1+c_halo)-c_halo/(1+c_halo));	// NFW Phi calc from Lokas&Mamon 2001
	Phi_out = -g_c*log(1+c_halo*GALSF_SUBGRID_VBOOST)/GALSF_SUBGRID_VBOOST;
	Phi_in = -g_c*log(1+c_halo*s)/s;
	vw_new = (Phi_out-Phi_in)*v_c*v_c+0.5*v_wind*v_wind;  // new_KE = delta-Phi + old_KE
	vw_new = sqrt(2*vw_new);
	//if( ThisTask==0 ) printf("VBOOST: z= %g M*= %g Mh= %g vold= %g vnew= %g boost= %g s= %g\n",1./All.cf_atime-1,Groups[P[i].galID].Mbar,Mhalo,v_wind,vw_new,vw_new/v_wind,s);

	return vw_new;
}

#endif

#ifdef SIGMASFR_WIND
int fof_compare_r2(const void *a, const void *b)
{
  if(((fof_part_data *) a)->r2 < ((fof_part_data *) b)->r2) return -1;
  if(((fof_part_data *) a)->r2 > ((fof_part_data *) b)->r2) return +1;
  return 0;
}
#endif

int fof_compare_galID(const void *a, const void *b)
{
  if(((fof_part_data *) a)->galID < ((fof_part_data *) b)->galID) return -1;
  if(((fof_part_data *) a)->galID > ((fof_part_data *) b)->galID) return +1;
  return 0;
}

int fof_compare_xdir(const void *a, const void *b)
{
  if(((fof_part_data *) a)->Pos[0] < ((fof_part_data *) b)->Pos[0]) return -1;
  if(((fof_part_data *) a)->Pos[0] > ((fof_part_data *) b)->Pos[0]) return +1;
  return 0;
}

int fof_compare_ydir(const void *a, const void *b)
{
  if(((fof_part_data *) a)->Pos[1] < ((fof_part_data *) b)->Pos[1]) return -1;
  if(((fof_part_data *) a)->Pos[1] > ((fof_part_data *) b)->Pos[1]) return +1;
  return 0;
}

int fof_compare_zdir(const void *a, const void *b)
{
  if(((fof_part_data *) a)->Pos[2] < ((fof_part_data *) b)->Pos[2]) return -1;
  if(((fof_part_data *) a)->Pos[2] > ((fof_part_data *) b)->Pos[2]) return +1;
  return 0;
}

double fofrad_periodic(double x)
{
  if(x >= 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;
  return x;
}


#endif // FOFRAD
