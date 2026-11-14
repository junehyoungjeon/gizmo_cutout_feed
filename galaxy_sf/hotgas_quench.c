#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*!
 *  routines for doing quenching via keeping hot halo gas hot, as in Gabor+12
 */

#ifdef HOTGAS_QUENCH // master switch for compiling the routines below //

/* Quenching model params */
#define HOTGAS_QUENCH_MHALO	0.0	// halo mass above which to quench, in units of 10^12 Mo; if 0 then use Mitra+14 z-dependent parameterization
#define HOTGAS_QUENCH_DENS	0.0   	// in atoms/cm^3; only heat gas *below* this density; if set to 0, heats everything that is not self-shielded
#define HOTGAS_QUENCH_TFACT	1.2	// heat gas to this factor times the virial temperature
#define HOTGAS_QUENCH_TLOW	1.0  	// heat gas if it drops below this times Tvir
#define HOTGAS_QUENCH_SLOWHEAT	1000.0	// maximum factor to heat particle per timestep
#define HOTGAS_QUENCH_TLIM	0.	// only heat halo gas above this temperature

/* miscellaneous definitions */
static double yhelium = (1 - HYDROGEN_MASSFRAC) / (4 * HYDROGEN_MASSFRAC);
static double mhboltz = PROTONMASS / BOLTZMANN;         /* hydrogen mass over Boltzmann constant */
static int firsttime=1,firsthotgas=1;

void hotgas_quenching()
{
	int i,k,Nquench=0,Nquenchtot;
	double Mhalo,Mquench,VirTemp,uVir,u_old,T_old,HeatFactor,dEnergy=0,maxmass=0,dEnergytot;
	double massfac = All.UnitMass_in_g / (SOLAR_MASS*All.HubbleParam) * 1.e-12;	// converts mass to 1e12 Mo units
	double densfac = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS; // converts density to physical atoms/cm^3
	double ufac = (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs);  // to cgs

  	if( firsthotgas==1 ) {
	  if( ThisTask==0 ) printf("HOTGAS_QUENCH: MHALO=%g DENS=%g TFACT=%g TLOW=%g SLOWHEAT=%g TLIM=%g\n", HOTGAS_QUENCH_MHALO,HOTGAS_QUENCH_DENS,HOTGAS_QUENCH_TFACT,HOTGAS_QUENCH_TLOW,HOTGAS_QUENCH_SLOWHEAT,HOTGAS_QUENCH_TLIM);
	  firsthotgas = 0;
	}
	if( HOTGAS_QUENCH_MHALO == 0.0 ) Mquench = 0.96 + 0.48*(1./All.Time-1);	// Mitra+14
	else Mquench = HOTGAS_QUENCH_MHALO;

	for( i=0; i<NumPart; i++ ) {
	  if( P[i].Metallicity[0] != P[i].Metallicity[0] ) for(k=0;k<NUM_METAL_SPECIES;k++) P[i].Metallicity[k]=0.;
	  if( SphP[i].HostHaloMass == 0 ) continue; // only heat gas in halos
	  if( maxmass < SphP[i].HostHaloMass*massfac ) maxmass = SphP[i].HostHaloMass*massfac;
	  if( HOTGAS_QUENCH_DENS > 0.0 && Particle_density_for_energy_i(i)*densfac > HOTGAS_QUENCH_DENS ) continue; // only heat gas below this density

	  Mhalo = massfac * SphP[i].HostHaloMass; 	// halo mass in units of 1e12 Mo
	  if( Mhalo < Mquench ) continue;		// only heat halos above Mquench

#ifdef HOTGAS_PRECIP
	  if( precip(i) ) continue;	// only heat non-precipitating gas
#else
	  u_old = SphP[i].InternalEnergy;
	  T_old = GAMMA_MINUS1 * u_old * ufac * mhboltz * (1+4*yhelium) / (1+SphP[i].Ne+yhelium);
	  if( T_old < HOTGAS_QUENCH_TLIM ) continue;	// only heat gas below TLIM
	  if( HOTGAS_QUENCH_DENS == 0.0 && coldgasfrac(i,T_old,Particle_density_for_energy_i(i)*densfac,0) > 0.1 ) continue;  // if density threshold=0, only heat non-self-shielded gas
#endif

	  VirTemp = 9.52e+07 * pow(Mhalo*1.e-3, 0.66666); /* 8.2 KeV at 1e15, from Voit (2005) */
  	  if( T_old < HOTGAS_QUENCH_TLOW*VirTemp ) {
	    SphP[i].Ne = 1+2*yhelium;	// fully ionized now
	    uVir = VirTemp / (ufac*GAMMA_MINUS1*mhboltz*(1+4*yhelium)/(1+SphP[i].Ne+yhelium));
	    HeatFactor = HOTGAS_QUENCH_TFACT * uVir / SphP[i].InternalEnergy;
	    if( HeatFactor > HOTGAS_QUENCH_SLOWHEAT ) HeatFactor = HOTGAS_QUENCH_SLOWHEAT;
	    SphP[i].InternalEnergy *= HeatFactor;	// new u
	    SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
	    SphP[i].Pressure = get_pressure(i);
	    dEnergy += (SphP[i].InternalEnergy-u_old)*P[i].Mass;  // sum added energy
	    Nquench++;
	    //if( ThisTask==0 ) printf("HOTGAS_QUENCH: z= %g ID= %d Mh= %g n= %g T= %g %g Tv= %g\n",1./All.Time-1,P[i].ID,Mhalo,Particle_density_for_energy_i(i)*densfac,log10(T_old),LogTemp(SphP[i].InternalEnergy*ufac,SphP[i].Ne),log10(VirTemp));
	  }
	}
	MPI_Allreduce(&Nquench, &Nquenchtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&dEnergy, &dEnergytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if( ThisTask==0 && Nquenchtot>0 ) printf("HOTGAS_QUENCH: z= %g Mq= %g Nq= %d Eq= %g\n",1./All.Time-1,Mquench,Nquench,dEnergy*All.UnitEnergy_in_cgs);
	return;
}

#ifdef HOTGAS_PRECIP
int precip(int i)
{
	double tcool;
	tcool = CallGrackle(DMAX(All.MinEgySpec, SphP[i].InternalEnergy), Particle_density_for_energy_i(i)*All.cf_a3inv, 0.0, SphP[i].Ne, i, 1);
	tff = P[i].FreeFallTime;

	if( tcool < tff*HOTGAS_PRECIP ) return 1; // true if precipitating, e.g. tcool/tff < 10
	else return 0;	
}
#endif

#endif // HOTGAS_QUENCHING

