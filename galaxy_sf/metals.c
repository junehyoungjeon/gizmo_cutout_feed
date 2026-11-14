#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif

/* Routines for doing metals in the instantaneous enrichment approximation */

/*
 * Romeel Dave, 4/8/15
 */


#if defined(METALS) && defined(GALSF_INSTANTANEOUS_METALS)

/* in case you're wondering, here are some conventions that may be useful for solar abundances
    All.SolarAbundances[0]=0.02;        // all metals (by mass); present photospheric abundances from Asplund et al. 2009 (Z=0.0134, proto-solar=0.0142) in notes;
                                        //   also Anders+Grevesse 1989 (older, but hugely-cited compilation; their Z=0.0201, proto-solar=0.0213)
    // note that the 'all metals' above is the only one where solar enters with any direct role; for everything else, 
    //  'solar' is totally arbitrary (the enrichment routines, etc, don't know what solar is 'supposed' to be): so these
    //  are here purely for convenience and to initialize non-zero metallicities
    //
    All.SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
    All.SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)
    All.SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)
    All.SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)
    All.SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
    All.SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
    All.SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)
    All.SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)
    All.SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
    All.SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)
*/

#define NMET_YIELDTABLE	4

#define DELAYED_FEEDBACK_NPART	16	/* distribute metals and energy from delayed feedback to at least this many surrounding non-wind gas particles */

static double yieldsII[NUM_METAL_SPECIES][NMET_YIELDTABLE],yieldsIa[NUM_METAL_SPECIES];
static double yieldmet[NMET_YIELDTABLE];
static double yields_factor[NMET_YIELDTABLE];	// =SUM(tracked metals)/total metals
static int firsttime=1;

double get_yield(double Zmet, int k, double age)
{
    int i,j;
    double agemax=0.03753; /* Gyr; age<agemax: Type II; else Type Ia yields */
    double yield,logZ;

    if( firsttime ) {
        yieldmet[0] = -4; yieldmet[1] = -3; yieldmet[2] = -2.4; yieldmet[3] = -1.7;
        if(NUM_METAL_SPECIES==11) {
        // All: He,C,N,O,Ne,Mg,Si,S,Ca,Fe.  yieldsIa holds mass of metals formed in 1 SNIa.
            /* SNIa */ /* from Iwamoto et al. 1999; 'W7' models */
            yieldsIa[0]=1.4;/* total metal mass; all H+He burned */
            yieldsIa[1]=0.0;/*He*/ yieldsIa[2]=0.049;/*C*/ yieldsIa[3]=1.2e-6;/*N*/ yieldsIa[4]=0.143;/*O*/
            yieldsIa[5]=0.0045;/*Ne*/ yieldsIa[6]=0.0086;/*Mg*/ yieldsIa[7]=0.156;/*Si*/
            yieldsIa[8]=0.087;/*S*/ yieldsIa[9]=0.012;/*Ca*/ yieldsIa[10]=0.743;/*Fe*/

            /* SNII (IMF-averaged... may not be the best approx on short timescales...) */
            // Woosley & Weaver 1995 'B' models, integrated (log-log interp) over Salpeter IMF, from 10-50 rescaled from 0.1-100 solar //
            //yields[0]=1.70;/*Z*/
            //yields[1]=4.03;/*He*/ yields[2]=0.117;/*C*/ yields[3]=0.0399;/*N*/ yields[4]=1.06;/*O*/
            //yields[5]=0.169;/*Ne*/ yields[6]=0.0596;/*Mg*/ yields[7]=0.0924;/*Si*/
            //yields[8]=0.0408;/*S*/ yields[9]=0.00492;/*Ca*/ yields[10]=0.0842;/*Fe*/
            /* note that the Mg yield here, and perhaps some other products, is systematically low compared to observations (by ~0.4 dex) */
            // Nomoto 2006 (arXiv:0605725) IMF-weighted yields for Salpeter 0.07-50 Mo 
            yieldsII[0][0]=0.0127; yieldsII[0][1]=0.0123; yieldsII[0][2]=0.0113; yieldsII[0][3]=0.0113;/*Z*/
            yieldsII[1][0]=0.023; yieldsII[1][1]=0.0226; yieldsII[1][2]=0.0222; yieldsII[1][3]=0.0205;/*He*/
            yieldsII[2][0]=7.51e-4; yieldsII[2][1]=4.14e-4; yieldsII[2][2]=5.17e-4; yieldsII[2][3]=6.74e-4;/*C*/
            yieldsII[3][0]=3.36e-6; yieldsII[3][1]=2.82e-5; yieldsII[3][2]=8.77e-5; yieldsII[3][3]=2.35e-4;/*N*/
            yieldsII[4][0]=8.03e-3; yieldsII[4][1]=8.05e-3; yieldsII[4][2]=6.58e-3; yieldsII[4][3]=6.14e-3;/*O*/
            yieldsII[5][0]=1.28e-3; yieldsII[5][1]=1.38e-3; yieldsII[5][2]=1.57e-3; yieldsII[5][3]=1.68e-3;/*Ne*/
            yieldsII[6][0]=5.31e-4; yieldsII[6][1]=2.55e-4; yieldsII[6][2]=5.04e-4; yieldsII[6][3]=6.08e-4;/*Mg*/
            yieldsII[7][0]=8.11e-4; yieldsII[7][1]=7.09e-4; yieldsII[7][2]=6.17e-4; yieldsII[7][3]=4.55e-4;/*Si*/
            yieldsII[8][0]=3.12e-4; yieldsII[8][1]=2.86e-4; yieldsII[8][2]=2.50e-4; yieldsII[8][3]=1.96e-4;/*S*/
            yieldsII[9][0]=3.58e-5; yieldsII[9][1]=3.34e-5; yieldsII[9][2]=3.02e-5; yieldsII[9][3]=2.39e-5;/*Ca*/
            yieldsII[10][0]=3.17e-4; yieldsII[10][1]=3.38e-4; yieldsII[10][2]=3.22e-4; yieldsII[10][3]=3.48e-4;/*Fe*/
	    for( i=0; i<NUM_METAL_SPECIES; i++ ) for( j=0; j<NMET_YIELDTABLE; j++ ) yieldsII[i][j] *= All.FactorSN/0.1;	// Assumes Salpeter; multiply by e.g. 1.8 to get yield wrt Kroupa/Chabrier IMF
	    for( i=2; i<NUM_METAL_SPECIES; i++ ) for( j=0; j<NMET_YIELDTABLE; j++ ) yields_factor[j] += yieldsII[i][j];
	    for( j=0; j<NMET_YIELDTABLE; j++ ) yields_factor[j] /= yieldsII[0][j];
	    if( ThisTask==0 ) {
		printf("YIELDS: Ratio of sum of tracked metals to total, for 4 ages:");
	    	for( j=0; j<NMET_YIELDTABLE; j++ ) printf("%g ",yields_factor[j]); 
		printf("\n");
	    }
        }
        else if(NUM_METAL_SPECIES==1) {
            yieldsIa[0]=1; 
            yieldsII[0][0]=0.0127; yieldsII[0][1]=0.0123; yieldsII[0][2]=0.0113; yieldsII[0][3]=0.0113;/*Z*/
	    for( j=0; j<NMET_YIELDTABLE; j++ ) yieldsII[0][j] *= All.FactorSN/0.1;
        }
	else {
	    if( ThisTask==0 ) printf("NUM_METAL_SPECIES= %d must be 1 or 11 in metals.c\n",NUM_METAL_SPECIES);
	    endrun(2384);
	}
	firsttime=0;
    }	// firsttime

    logZ = log10(Zmet+1.e-6);
    yield = 0.0;
    if( age <= agemax ) {
	if( logZ <= yieldmet[0] ) yield = yieldsII[k][0];
	else if( logZ >= yieldmet[NMET_YIELDTABLE-1] ) yield = yieldsII[k][NMET_YIELDTABLE-1];
	else {
	    i = 0;
	    while(  i < NMET_YIELDTABLE && logZ > yieldmet[i] ) i++;
	    yield = (yieldsII[k][i]*(logZ-yieldmet[i-1])+yieldsII[k][i-1]*(yieldmet[i]-logZ))/(yieldmet[i]-yieldmet[i-1]);
	//fprintf(stdout,"%d %d %g %g %g %g %g %g\n",k,i,yieldsII[k][i],(logZ-yieldmet[i-1]),yieldsII[k][i-1],(yieldmet[i]-logZ),(yieldmet[i]-yieldmet[i-1]),yield);
	}
    }
    else yield = yieldsIa[k];

    return yield; 
}

#if defined(GALSF_TYPEIA) || defined(GALSF_AGBFEEDBACK)

#ifdef GALSF_TYPEIA
#define TYPEIARATE_A 5.3e-14 /* param A with Mstar in Mo; Sullivan etal 06 */
#define TYPEIARATE_B 3.9e-4 /* param B with Mdot in Mo/yr; S06 */
#define TYPEIAENERGY 1.e51 /* ergs per Type Ia SN. */
#define TYPEIADELAY 0.7 /* Delay time before stars for Type Ia's, in Gyr; SB05 */
#endif

#ifdef GALSF_AGBFEEDBACK
float AGByield[2][24][150];
float AGBmdot[24][150];
#endif

#define PROMPT_TIME_RESOLUTION 1.e7 /* Frequency for doing prompt TypeIa feedback */
#define DELAYED_TIME_RESOLUTION 2.e7 /* Frequency for doing delayed SNIa+AGB feedback */

static int NstepDelayed=0;
static MyFloat TimeLastPrompt=0,TimeLastDelayed=0;	// time in Gyr since last call

void Delayed_Feedback(void)
{
  int i,j,k,n,nd=0;
  int nprompt=0,ndelayed=0,totprompt=0,totdelayed=0,nagb=0,totagb=0;
  int Nstep,imet,itime;
  int startnode, numngb, dummy;
  MyFloat dtime, tage, N_Ia=0, egySN=0, p, prob; 
  double egyPrompt=0, egyDelayed=0, totegyP, totegyD, egyAGB=0, Mmetsum=0, totMet=0, ufac, densfac, Mmet[NUM_METAL_SPECIES];
  double r2, h2, u, dx, dy, dz, hsml, hinv3, hinv4, wk=0, dwk=0, norm_j, sumwk; 
  double yield[NUM_METAL_SPECIES],AGBmet,Ztot,massj_fact;
  MyFloat *Ngbweight;
  MyFloat metal_loss[NUM_METAL_SPECIES], mass_loss;
  void InitAGBYields();

  if( NstepDelayed == 0 ) {	// initialization stuff
#ifdef GALSF_AGBFEEDBACK
    InitAGBYields();
#endif
    TimeLastDelayed=All.Time;
    TimeLastPrompt=All.Time;
    if( ThisTask==0 ) printf("All.OrigGasMass= %g (%g Mo) at a=%g\n",All.OrigGasMass,All.OrigGasMass*All.UnitMass_in_g/SOLAR_MASS/All.HubbleParam,All.Time);
    if( NUM_METAL_SPECIES < 10 ) endrun (7804);	// must be tracking full 11 metal species
    NstepDelayed = 1;
  }

  ufac = All.UnitPressure_in_cgs/All.UnitDensity_in_cgs;
  densfac = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;

  // get MIN/MAX NumPart among processors for logging purposes
  unsigned int MinNumPart = 1e6;
  unsigned int MaxNumPart = 0;
  MPI_Allreduce(&NumPart,&MinNumPart, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&NumPart,&MaxNumPart, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(ThisTask==0)
    printf("NUMPART a= %e  MIN= %d  MAX= %d\n",All.Time,MinNumPart,MaxNumPart);

/* Do prompt TypeIa feedback */
#ifdef GALSF_TYPEIA
  dtime = 1.e9 * evaluate_stellar_age_Gyr(TimeLastPrompt);	/* yr since prompt SNIa was done */
  if( dtime > PROMPT_TIME_RESOLUTION ) {
    for(i = 0, Mmetsum=0; i < NumPart; i++) {
      if(P[i].Type == 0 && SphP[i].Sfr > 0. ) {
        nprompt++;
        N_Ia = TYPEIARATE_B* SphP[i].Sfr * dtime / All.HubbleParam; /* # of type Ia's this timestep */
        egySN = N_Ia * TYPEIAENERGY / All.UnitEnergy_in_cgs;
        SphP[i].InternalEnergy += egySN / P[i].Mass;
        egyPrompt += egySN;
        for( k=0; k<NUM_METAL_SPECIES; k++) {
	  Mmet[k] = N_Ia * yieldsIa[k] * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g;	// metal mass to be added
	  P[i].Metallicity[k] = (P[i].Mass * P[i].Metallicity[k] + Mmet[k] ) / P[i].Mass;
	  if( k==0 ) Mmetsum += Mmet[k];
        }
        convert_u_to_temp(SphP[i].InternalEnergy*ufac, Particle_density_for_energy_i(i)*densfac,&SphP[i].Ne,i);  //update electron abundance
      }
    }
    MPI_Allreduce(&nprompt, &totprompt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&egyPrompt, &totegyP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Mmetsum, &totMet, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    TimeLastPrompt=All.Time;
    if( ThisTask == 0 && totprompt>0 ) fprintf(stdout,"DELAYED_FEEDBACK Prompt: a= %g dt/Myr= %g Nprompt=%d (%g) MZtot= %g\n",All.Time,dtime*1.e-6,totprompt,totegyP*All.UnitEnergy_in_cgs/(totprompt+1.e-20)/TYPEIAENERGY,totMet*All.UnitMass_in_g/(SOLAR_MASS * All.HubbleParam)/totprompt);
  }
#endif

/* Do delayed feedback from star particles, both TypeIa and AGB */
  dtime = 1.e9 * evaluate_stellar_age_Gyr(TimeLastDelayed);	/* yr since delayed was done */
  if( dtime > DELAYED_TIME_RESOLUTION ) {
    NstepDelayed++;
    Ngblist = (int *) mymalloc("delayed_feedback neighbor list",NumPart * sizeof(int));
    Ngbweight = (MyFloat *) mymalloc("delayed_feedback neighbor weights",NumPart * sizeof(MyFloat));

    for(i = 0, Mmetsum=0; i < NumPart; i++) {
      if(P[i].Type != 4 ) continue;
      tage = evaluate_stellar_age_Gyr(P[i].StellarAge);
/* Check whether delayed feedback should occur at tage */
      if( tage < 0.02 ) continue;	// no delayed feedback at <20 Myr
#ifndef GALSF_AGBFEEDBACK
      if( tage < TYPEIADELAY ) continue;	// if no AGB, no delayed SNIa before TYPEIADELAY Gyr
#endif
//      dtime = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a * (All.UnitTime_in_s / SEC_PER_YEAR); // dloga to yrs
      Nstep = (int) (tage/0.2)+1;	// set how often to do delayed feedback from a given star
      if( Nstep > 10 ) Nstep = 10;
      if( NstepDelayed % Nstep != 0 ) continue;	// only do this star every Nstep steps

#ifdef GALSF_TYPEIA
      if( tage > TYPEIADELAY ) {
        ndelayed++;
        N_Ia = (P[i].Mass/All.HubbleParam) * dtime * TYPEIARATE_A * (All.UnitMass_in_g / SOLAR_MASS) * Nstep; /* # of Type Ia's this timestep */
        egySN = N_Ia * TYPEIAENERGY / All.UnitEnergy_in_cgs;  /* Energy from each SNIa in code units*/
      }
#endif
#ifdef GALSF_AGBFEEDBACK
      nagb++;
      itime = (int)((log10(tage*1.e9)-7.18)/0.02);
      imet = (int) ((log10(P[i].Metallicity[0]+pow(10,-4.8))+4.8)/0.2); 
      if(imet<0 || itime < 0 || imet>22 || itime > 148){
	printf("P%d AGB WARNING: IMET OR ITIME OUT OF RANGE %8.6f id= %8d imet=%d itime=%d ",ThisTask,All.Time,P[i].ID,imet,itime);
	printf("%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %8.6e\n",P[i].Mass,P[i].Metallicity[0],P[i].Metallicity[2],P[i].Metallicity[4],tage*1.e9, dtime, AGBmdot[imet][itime]);
	if( imet > 22 ) imet = 22;
	else endrun(7803);
      }
/* Get mass loss over dtime from lookup tables */
      p = AGBmdot[imet][itime]*dtime;
      prob = 1 - exp(-p);
      mass_loss = prob * P[i].Mass * Nstep;
      if( mass_loss == 0 ) continue;	// no mass or metal loss; no point in doing this particle

/* Compute metal mass loss.  */
      for(k=2;k<NUM_METAL_SPECIES;k++) yield[k] = P[i].Metallicity[k]; 
/* If the AGByield table value is listed as 2 (i.e. >1), this signifies that the yield should be the original metallicity itself; only if it is <1 do we use the value from the AGByield tables */
      yield[1] = 0.36;	// helium mass fraction of AGB star ejecta
      if(AGByield[0][imet][itime]<1.0) yield[2] = AGByield[0][imet][itime];  	// carbon
      yield[3] = 0.0041;	// nitrogen
      if(AGByield[1][imet][itime]<1.0) yield[4] = AGByield[1][imet][itime];  	// oxygen
      yield[0]=0.0; for(k=2;k<NUM_METAL_SPECIES;k++) yield[0]+=yield[k];
      for(k=0;k<NUM_METAL_SPECIES;k++) metal_loss[k] = yield[k] * mass_loss;
	//printf("DELAYED %d %d %g %g %g %g %g %g %g\n",imet,itime,P[i].Metallicity[0],tage,mass_loss/P[i].Mass,mass_loss*All.UnitMass_in_g/(SOLAR_MASS * All.HubbleParam),yield[2],yield[2]*mass_loss*All.UnitMass_in_g/(SOLAR_MASS * All.HubbleParam),P[i].Metallicity[2]*mass_loss*All.UnitMass_in_g/(SOLAR_MASS * All.HubbleParam));
#endif

/* Find neighboring gas particles to add SN heat, metals */
      hsml = All.MinHsml;
      do {
	hsml *= 1.05;
	startnode = All.MaxPart;
	numngb = ngb_treefind_variable_threads(P[i].Pos,hsml,-1,&startnode,0,&dummy,&dummy,&dummy,Ngblist);
	for( n=0, nd=0, sumwk=0; n<numngb; n++ ) {
	  j = Ngblist[n];
	  Ngbweight[n] = 0;
#ifdef GALSF_SUBGRID_WINDS
	  if(SphP[j].DelayTime > 0) {
	    nd++;	// wind particle; ignore
	    continue;
	  }
#endif
	  if( numngb >= DELAYED_FEEDBACK_NPART ) {
	    /*
            dx = fabs(P[i].Pos[0] - P[j].Pos[0]);
            dy = fabs(P[i].Pos[1] - P[j].Pos[1]);
            dz = fabs(P[i].Pos[2] - P[j].Pos[2]);
	    */
	    dx = P[i].Pos[0] - P[j].Pos[0];
	    dy = P[i].Pos[1] - P[j].Pos[1];
	    dz = P[i].Pos[2] - P[j].Pos[2];
#ifdef PERIODIC   /*  distance to the closest image in the given box size  */
	    /*
            if(dx > 0.5 * All.BoxSize) dx = All.BoxSize-dx;
            if(dy > 0.5 * All.BoxSize) dy = All.BoxSize-dy;
            if(dz > 0.5 * All.BoxSize) dz = All.BoxSize-dz;
	    */
	    NEAREST_XYZ(dx,dy,dz,1);
#endif
            r2 = dx * dx + dy * dy + dz * dz;
	    h2 = hsml*hsml;
	    if( r2 < h2 ) {
	      u = sqrt(r2/h2);
	      hinv4 = 1/(h2*h2);
	      hinv3 = hinv4*hsml;
	      kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
	      sumwk += wk;
	      Ngbweight[n] = wk;
	    }
	    else {
	      nd++;	 /* particle is too far away; ignore */
	      Ngbweight[n] = 0;
	    }
	  }
	}
	numngb -= nd;
      } while( numngb < DELAYED_FEEDBACK_NPART && hsml < 0.5*All.BoxSize);

      //if(numngb<DELAYED_FEEDBACK_NPART) printf("DELAYED_FEEDBACK WARNING NOT ENOUGH NEIGHBORS: ID= %d numngb= %d nd= %d pos= %g %g %g\n",P[i].ID,numngb,nd,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
      if(numngb<DELAYED_FEEDBACK_NPART){
	printf("DELAYED_FEEDBACK SUB-OPTIMAL NEIGHBOR COUNT TASK %d: ID= %u numngb= %d NumPart= %d nd= %d pos= %g %g %g\n",
	       ThisTask,P[i].ID,numngb,NumPart,nd,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
	continue;  // SKIP!
      }

/* Loop over neighbors to distribute energy+metals in a kernel-weighted manner */
      for( n=0; n<numngb+nd; n++ ) {
	if( Ngbweight[n] > 0 ) {
	  j = Ngblist[n];
	  norm_j = Ngbweight[n]/sumwk;
#ifdef GALSF_AGBFEEDBACK	/* This must be done before Type Ia so that yields are correct */
#ifdef GALSF_AGBWINDHEATING   // this variable should be set to the AGB wind vel in km/s (Conroy+14)
	  for(k=0,egyAGB=0;k<3;k++) egyAGB += (P[i].Vel[k]-P[j].Vel[k])*(P[i].Vel[k]-P[j].Vel[k]);
	  egyAGB = 0.5*mass_loss*(egyAGB+GALSF_AGBWINDHEATING*GALSF_AGBWINDHEATING*1e10*All.Time/(All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s));	
	  SphP[j].InternalEnergy += norm_j*egyAGB/P[j].Mass;	// add kernel-weighted energy
#endif
          for( k=2,Ztot=0; k<NUM_METAL_SPECIES; k++) Ztot+=P[j].Metallicity[k];
	  massj_fact = 1./(P[j].Mass+norm_j*mass_loss);
          for( k=0; k<NUM_METAL_SPECIES; k++) {	// add AGB ejectae from i to particle j
	    AGBmet = (P[j].Metallicity[k]*P[j].Mass+norm_j*metal_loss[k])*massj_fact;
	    P[j].Metallicity[k] = AGBmet;
	  }
	  for(k=0;k<3;k++) P[j].Vel[k] = (P[j].Mass*P[j].Vel[k] + norm_j*mass_loss*P[i].Vel[k])*massj_fact;	// conserve momentum
	  P[j].Mass += norm_j*mass_loss;
	  P[i].Mass -= norm_j*mass_loss;	// conserve mass 
          for( k=2,Ztot=0; k<NUM_METAL_SPECIES; k++) Ztot+=P[j].Metallicity[k];
	  //if( P[j].Metallicity[0]<Ztot ) fprintf(stdout,"METAGB TROUBLE: %d %d %g %g %g %g\n",i,j,P[i].Metallicity[0],Ztot,P[j].Metallicity[0],mass_loss);
#endif
#ifdef GALSF_TYPEIA
          if( tage > TYPEIADELAY ) {
            SphP[j].InternalEnergy += norm_j * egySN / P[j].Mass;  // add kernel-weighted SNIa energy
	    egyDelayed += norm_j * egySN;
            for( k=2,Mmetsum=0; k<NUM_METAL_SPECIES; k++) {  // compute metal mass to be exchanged
	      Mmet[k] = norm_j * N_Ia * yieldsIa[k] * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g; // metal mass to be added to neighbor j
	      Mmetsum += Mmet[k];
	    }
	    massj_fact = 1./(P[j].Mass+Mmetsum);
            for( k=2,Mmetsum=0; k<NUM_METAL_SPECIES; k++) {  // now add metals to j
	      P[j].Metallicity[k] = (P[j].Mass * P[j].Metallicity[k] + Mmet[k] ) * massj_fact;
            }
	    P[j].Metallicity[0] = (P[j].Mass * P[j].Metallicity[0] + Mmetsum ) * massj_fact;
	    P[j].Mass += Mmetsum;
	    P[i].Mass -= Mmetsum;	// conserve mass 
	    for(k=0;k<3;k++) P[j].Vel[k] = (P[j].Mass*P[j].Vel[k] + norm_j*Mmetsum*P[i].Vel[k])*massj_fact;	// conserve momentum
	  }
          for( k=2,Ztot=0; k<NUM_METAL_SPECIES; k++) Ztot+=P[j].Metallicity[k];
#endif
          convert_u_to_temp(SphP[j].InternalEnergy*ufac, Particle_density_for_energy_i(j)*densfac,&SphP[j].Ne,j);  /* update electron abundance for new thermal energy */
	}
      }
      if( P[i].Mass < 0.2*All.OrigGasMass/GALSF_GENERATIONS ) endrun(7805);  // should never have this much mass loss
    }
    myfree(Ngbweight);
    myfree(Ngblist);
    MPI_Allreduce(&ndelayed, &totdelayed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nagb, &totagb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&egyDelayed, &totegyD, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Mmetsum, &totMet, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if( ThisTask == 0 && (totdelayed+totagb)>0 ) fprintf(stdout,"DELAYED_FEEDBACK Delayed a= %g dt/Myr= %g NIa= %d (%g) Nagb= %d MZIa= %g\n",All.Time,dtime*1.e-6,totdelayed,totegyD*All.UnitEnergy_in_cgs/(totdelayed+1.e-20)/TYPEIAENERGY,totagb,totMet/ (SOLAR_MASS * All.HubbleParam) * All.UnitMass_in_g/totagb);  fflush(stdout);

    TimeLastDelayed=All.Time;
  }	// done delayed feedback (dtime > DELAYED_TIME_RESOLUTION)

  return;
}

#ifdef GALSF_AGBFEEDBACK
void InitAGBYields(void)
{
  int i,k;
  char buf[100],line[300];
  float age;
  FILE *FdAGB;

  sprintf(buf, "AGBmdot.tab");
  if(ThisTask==0)printf("AGB mass loss table: %s\n",buf);
  if( (FdAGB = fopen(buf, "r")) == NULL){
    printf("No %s file, find this and put it in\n",buf);
    exit(-1000);
  }
  i = 0;
  while( fgets(line,300,FdAGB) != NULL ) {
    if( strstr(line,"#") != NULL ) continue;
    sscanf(line,"%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&age,&AGBmdot[0][i],&AGBmdot[1][i],&AGBmdot[2][i],&AGBmdot[3][i],&AGBmdot[4][i],&AGBmdot[5][i],&AGBmdot[6][i],&AGBmdot[7][i],&AGBmdot[8][i],&AGBmdot[9][i],&AGBmdot[10][i],&AGBmdot[11][i],&AGBmdot[12][i],&AGBmdot[13][i],&AGBmdot[14][i],&AGBmdot[15][i],&AGBmdot[16][i],&AGBmdot[17][i],&AGBmdot[18][i],&AGBmdot[19][i],&AGBmdot[20][i],&AGBmdot[21][i],&AGBmdot[22][i],&AGBmdot[23][i]);
    i++;
  }
  fclose(FdAGB);

  for(k=0;k<2;k++){
    if(k==0)sprintf(buf, "AGBcyield_v2.tab");
    if(k==1)sprintf(buf, "AGBoyield_v2.tab");
    if(ThisTask==0)printf("AGB yield table: %s\n",buf);
    if( (FdAGB = fopen(buf, "r")) == NULL){
      printf("No %s file, find this and put it in\n",buf);
      exit(-1000);
    }
    i=0;
    while( fgets(line,300,FdAGB) != NULL ) {
      if( strstr(line,"#") != NULL ) continue;
      sscanf(line,"%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&age,&AGByield[k][0][i],&AGByield[k][1][i],&AGByield[k][2][i],&AGByield[k][3][i],&AGByield[k][4][i],&AGByield[k][5][i],&AGByield[k][6][i],&AGByield[k][7][i],&AGByield[k][8][i],&AGByield[k][9][i],&AGByield[k][10][i],&AGByield[k][11][i],&AGByield[k][12][i],&AGByield[k][13][i],&AGByield[k][14][i],&AGByield[k][15][i],&AGByield[k][16][i],&AGByield[k][17][i],&AGByield[k][18][i],&AGByield[k][19][i],&AGByield[k][20][i],&AGByield[k][21][i],&AGByield[k][22][i],&AGByield[k][23][i]);
      i++;
    }
    fclose(FdAGB);
  }
}
#endif /* GALSF_AGBFEEDBACK */

#endif /* TYPEIA || GALSF_AGBFEEDBACK */

#endif /* METALS */

