#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


#if defined(POPII_P2L) && defined(GALSF) && defined(METALS)

//ergs per event
#define ENERGY_TYPE2 1.0e51
#define ENERGY_PISN  1.0e52
#define ENERGY_BH    0.0
#define ENERGY_OTHER 0.0

// Bubble properties
#ifndef HII_DENSITY_ITERATION
#define HII_DENSITY_ITERATION 0
#endif
#define HII_DENSITY_POPII 1.0 //0.33
#define I_FRONT_RADIUS 0.146 //0.33     //kpc PHYSICAL
#define HII_redshift 10.0
#define TEMP_BUBBLE 2.0e4  //Kelvin
#define P2_ETOT 1.1e52
#define P2_METAL_RADIUS 0.8     //kpc PHYSICAL
#define P2_METAL_MASS 16.0     //Msun

// yield per type of SN event
#define YIELD_TYPE2 0.10
#define YIELD_PISN  0.50
#define YIELD_BH    0.0
#define YIELD_OTHER 0.0

#define LOWER 1.0     //lower limit to IMF
#define UPPER 150.0   //upper limit to IMF
#define DM 0.1        //dM for IMF table
//#define N_VALS 1490   //(UPPER-LOWER)/DM
#define N_VALS (int)((UPPER-LOWER)/DM)




void P2L(int i, int mode, double popII_mass, int LW_flag, double dt_years)
{
  double mass_of_star;
  int n, j, numngb_metal, numngb_Ifront, startnode, dummy;
  double energy_radius, metal_radius, energy_increase, total_metal_mass;
  double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
  double massj_fact, norm_j;
  double ufac = All.UnitPressure_in_cgs/All.UnitDensity_in_cgs;
  double mean_density, I_front_radius, temp_density;
  double n2code = PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
  double totmass, sumwk = 0.0;
  MyFloat *Ngbweight;
  double dx, dy, dz, ngb_r2;
#ifdef HII_EVOL
  double densfac = pow(1.0/(All.Time*(1.0+HII_redshift)), 3);
#else
  double densfac = 1.0;
#endif
  
  //printf("############################### POPII P2L called at PID: %d, heating_time = %e, P2L_flag = %0.01f ################################## \n",i,SphP[i].PE_heating_time,SphP[i].P2L_flag);
  //fflush(stdout);

  // Only identify particle within radius of initial star burst and set PE_heating_time.  NO FEEDBACK
  if(mode==0){
  
#ifndef NSFDBK
  energy_radius = I_FRONT_RADIUS*All.HubbleParam * pow(popII_mass / POPII_MASS /densfac/densfac, 1.0/3.0) / All.Time;  //to COMOVING
  startnode = All.MaxPart;
  Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
  numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);

/* reset radius
  I_front_radius = I_FRONT_RADIUS;
  for(int k=0; k<HII_DENSITY_ITERATION; k++){
    temp_density = 0.0;
    mean_density = 0.0;
    for(n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      temp_density += 1.0/(SphP[j].Density * All.cf_a3inv);
    }
    mean_density = numngb_Ifront/temp_density;
    I_front_radius = I_front_radius * pow(mean_density/HII_dens_ref, -2.0/3.0);
    myfree(Ngblist_Ifront);

    energy_radius = I_front_radius / All.Time;  //to COMOVING
    startnode = All.MaxPart;
    Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
    numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);
  }
*/
  /*                    loop through neighbors and set PE_heating_time             */
#ifdef UCHII_region
	double r2min = pow(energy_radius, 2.0);
	double UCHIIdens = HII_DENSITY_POPII;
	double r2UCHII = r2min;
  	for( n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      if(P[j].Type==0)
      {
    	dx = P[i].Pos[0] - P[j].Pos[0];
    	dy = P[i].Pos[1] - P[j].Pos[1];
    	dz = P[i].Pos[2] - P[j].Pos[2];
#ifdef PERIODIC   /*  distance to the closest image in the given box size  */
    	NEAREST_XYZ(dx,dy,dz,1);
#endif                                         
    	ngb_r2 = dx*dx + dy*dy + dz*dz;
    	if(ngb_r2 < r2min)
    	{
    		r2min = ngb_r2;
    		UCHIIdens = SphP[j].Density * All.cf_a3inv/n2code;
    	}
    }
    }
    r2UCHII = pow(energy_radius * pow(UCHIIdens/HII_DENSITY_POPII, -2.0/3.0), 2.0);
#endif

  for( n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      //if (SphP[j].Density * All.cf_a3inv <= )
      //
      if(P[j].Type==0)
      {
#if defined(DCBH) && defined(HII_FORBIDDEN) 
      	P[j].HighDens = 1;
#endif
#ifdef UCHII_region
    	dx = P[i].Pos[0] - P[j].Pos[0];
    	dy = P[i].Pos[1] - P[j].Pos[1];
    	dz = P[i].Pos[2] - P[j].Pos[2];
#ifdef PERIODIC   /*  distance to the closest image in the given box size  */
    	NEAREST_XYZ(dx,dy,dz,1);
#endif                                         
    	ngb_r2 = dx*dx + dy*dy + dz*dz;
    	if (ngb_r2 < r2UCHII) //&& (PPP[j].Hsml<))
    	{
    		P[j].InHII = 1;
/*
#ifdef ON_THE_FLY_FDBK
    		SphP[j].P2_heating_time = DMAX(dt_years, SphP[j].P2_heating_time);
#else	
			SphP[j].P2_heating_time = DMAX(POPII_PE_HEATING - P[i].SF_time, SphP[j].P2_heating_time);
#endif
*/
    	}
    	else if(P[j].InHII = 0)
#endif
#ifdef ON_THE_FLY_FDBK
    		SphP[j].P2_heating_time = DMAX(dt_years, SphP[j].P2_heating_time);
#else	
			SphP[j].P2_heating_time = DMAX(POPII_PE_HEATING - P[i].SF_time, SphP[j].P2_heating_time);
#endif
      }
  }
  myfree(Ngblist_Ifront);
#endif // NSFDBK

#ifdef LYMAN_WERNER   
#ifdef LW_LOCAL
  if(LW_flag==1){
  int numngb_LW;
  double rLW;
  double zred = 1.0/All.Time - 1.0;
  double eta_LW = 2e4;
  double J21_bg = LW_BACKGROUND(eta_LW, zred);
  double LW_radius = LW_MAX_RAD2;
  double LW_fac = 1.0;
  if(J21_bg>0){
  LW_radius = LW_rad_est(eta_LW/5.0, popII_mass, LW_P2T, J21_bg, LW_FAC);  
  LW_radius = DMIN(LW_radius, LW_MAX_RAD2)*All.HubbleParam/ All.Time;  //to COMOVING
  }
  
  startnode = All.MaxPart;
  Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
  numngb_LW = ngb_treefind_variable_threads(P[i].Pos, LW_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);
  
  for(n=0; n<numngb_LW; n++){
  	j = Ngblist_Ifront[n];
#ifdef GALSF_SUBGRID_WINDS
	if(P[j].Type==0)
	{
		if(SphP[j].DelayTime > 0) continue;
	}
#endif
	{
    dx = P[i].Pos[0] - P[j].Pos[0];
    dy = P[i].Pos[1] - P[j].Pos[1];
    dz = P[i].Pos[2] - P[j].Pos[2];
#ifdef PERIODIC   /*  distance to the closest image in the given box size  */
    NEAREST_XYZ(dx,dy,dz,1);
#endif                                         
    ngb_r2 = dx*dx + dy*dy + dz*dz;
    rLW = DMAX(sqrt(ngb_r2), P[i].Hsml) * All.Time; // to physical

#ifdef LW_OTF
	if (P[j].LW_time<dt_years)
	{
		P[j].JLW_local *= P[j].LW_time/dt_years;
		P[j].JLW_local += J21_local(eta_LW/5.0, popII_mass, LW_P2T, rLW);
		P[j].LW_time = dt_years;	
	}
	else
		P[j].JLW_local += J21_local(eta_LW/5.0, popII_mass, LW_P2T, rLW) * dt_years/P[j].LW_time;
#else
	if(LW_P2T > P[j].LW_time)
    	P[j].JLW_local = P[j].JLW_local * P[j].LW_time / LW_P2T;
    else
    	LW_fac = LW_P2T / P[j].LW_time;
    P[j].JLW_local += J21_local(eta_LW/5.0, popII_mass, LW_P2T, rLW) * LW_fac;
    P[j].LW_time = DMAX(LW_P2T - P[j].SF_time, P[j].LW_time);
#endif
    if(P[j].JLW_local < 0) P[j].JLW_local = 0.0;
    }
  }
  myfree(Ngblist_Ifront);
  }
#endif 
#endif
  
//#endif // NSFDBK
}//end mode=0


// if mode=1 then do the thermal and metal feedback at star death
if(mode==1){
  /*                              P2L METAL                            */
  double total_energy = popII_mass / POPII_MASS * P2_ETOT;
  double log_radius = -17.02 + 0.383 * log10(total_energy);
  metal_radius = pow(10, log_radius)*All.HubbleParam /1000/ All.Time;
  
  Ngblist = (int *) mymalloc("popII metal",NumPart * sizeof(int));  
  Ngbweight = (MyFloat *) mymalloc("popIII weights",NumPart * sizeof(MyFloat));
  
  startnode = All.MaxPart;
  numngb_metal = ngb_treefind_variable_threads(P[i].Pos, metal_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist);
  total_metal_mass = P2_METAL_MASS * Msun_to_Internal * popII_mass/POPII_MASS;

#ifdef POPII_KERNEL_DIST
  double wk=0.0, dwk=0.0;
  //double dx,dy,dz,ngb_r2,
  double r2,u,rinv4,rinv3;
  double radius_squared = metal_radius * metal_radius;

  //loop through neighbors and calculate kernel contribution
  for( n=0; n<numngb_metal; n++ ){
    j = Ngblist[n];  //index of neighbor n

    dx = P[i].Pos[0] - P[j].Pos[0];
    dy = P[i].Pos[1] - P[j].Pos[1];
    dz = P[i].Pos[2] - P[j].Pos[2];
#ifdef PERIODIC   /*  distance to the closest image in the given box size  */
    NEAREST_XYZ(dx,dy,dz,1);
#endif
                                                
    ngb_r2 = dx*dx + dy*dy + dz*dz;
    r2 = radius_squared;
    
    u = sqrt(ngb_r2/r2);
    rinv4 = 1/(r2*r2);
    rinv3 = rinv4*radius;
    kernel_main(u, rinv3, rinv4, &wk, &dwk, -1);
    sumwk += wk;
    Ngbweight[n] = wk;
  }
#else
  for( n=0; n<numngb_metal; n++ )
  {
  	j = Ngblist[n];
#ifndef POLLUTION
  	if(P[j].Type==0)
#else
	if(P[j].Type==0 || P[j].Type==3 || P[j].Type==4)
#endif
  	{
    	Ngbweight[n] = P[j].Mass;
    	sumwk += P[j].Mass;
    }
  }
  //sumwk = (double)numngb_metal;
#endif

  for( n=0; n<numngb_metal; n++){
      j = Ngblist[n];
      // do not add metals to the original star formation particle!!
      //if(j!=i){
#ifndef POLLUTION
  	if(P[j].Type==0)
#else
	if(P[j].Type==0 || P[j].Type==3 || P[j].Type==4)
#endif
      {
      //norm_j = 1./(double)numngb_metal;
      norm_j = Ngbweight[n]/sumwk;
      massj_fact = 1./(P[j].Mass+total_metal_mass * norm_j);

      //add to PopII metal tracer
      P[j].Metallicity[1] = (P[j].Mass * P[j].Metallicity[1] + total_metal_mass*norm_j) * massj_fact;
      //add to Total metal tracer
      P[j].Metallicity[0] = (P[j].Mass * P[j].Metallicity[0] + total_metal_mass*norm_j) * massj_fact;
      }
    }  
  myfree(Ngbweight);
  myfree(Ngblist);

#ifndef NSFDBK_SN
#ifdef POPII_SN
  /*                              P2L THERMAL                            */
  I_front_radius = I_FRONT_RADIUS*All.HubbleParam * pow(popII_mass / POPII_MASS /densfac/densfac, 1.0/3.0);
  //energy_radius = I_front_radius / All.Time;
  energy_radius = DMAX(I_front_radius / All.Time, metal_radius);  //to COMOVING
  startnode = All.MaxPart;
  Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
  numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);

  for(int k=0; k<HII_DENSITY_ITERATION; k++){
    temp_density = 0.0;
    mean_density = 0.0;
    totmass = 0.0;
    for(n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      if(P[j].Type==0){
	    temp_density += P[j].Mass/SphP[j].Density;
	    totmass += P[j].Mass;
	  }
    }
    if(totmass<=0) break; 
    mean_density = totmass * All.cf_a3inv/temp_density;
    I_front_radius = I_front_radius * pow(mean_density/(HII_DENSITY_POPII*n2code*densfac), -2.0/3.0);
    myfree(Ngblist_Ifront);

	energy_radius = DMIN(I_front_radius/All.Time, metal_radius);
    //energy_radius = DMAX(I_front_radius/All.Time, metal_radius);  //to COMOVING
    startnode = All.MaxPart;
    Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
    numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);
    if(I_front_radius >= metal_radius) break;
  }

  /*                    loop through neighbors and give THERMAL feedback              */
  for(n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      // no feedback on original star formation particle
    if(j!=i && P[j].Type==0 ){
      //if((SphP[j].Density * All.cf_a3inv) <= (All.HII_thresh*n2code))
      //if((SphP[j].Density * All.cf_a3inv) <= (BH_THRESH*n2code))
      {
      //SphP[j].InternalEnergy += energy_increase / All.UnitEnergy_in_cgs / P[j].Mass;
      SphP[j].Primordial_Chem[0] = 0.0;
      SphP[j].Primordial_Chem[1] = 1.0;
      SphP[j].Primordial_Chem[6] = 0.0;
      SphP[j].Primordial_Chem[7] = 1.0;
      SphP[j].Primordial_Chem[9] = 0.0;
      SphP[j].Primordial_Chem[10] = 1.0;
      SphP[j].Primordial_Chem[5] = 1.0/0.93;
      SphP[j].Primordial_Chem[2] = 0.0;
      SphP[j].Primordial_Chem[3] = 0.0;
      SphP[j].Primordial_Chem[4] = 0.0;
      SphP[j].Primordial_Chem[11] = 0.0;
      energy_increase = (TEMP_BUBBLE*BOLTZMANN)/(GAMMA_MINUS1*PROTONMASS *ufac);
      if (SphP[j].InternalEnergy < energy_increase)
        SphP[j].InternalEnergy = energy_increase; //test to set all to 5e4K
      SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;
      SphP[j].Pressure = get_pressure(j);
#if defined(DCBH) && defined(HII_FORBIDDEN) 
      	P[j].HighDens = 1;
#endif
      }
    }
  }
  myfree(Ngblist_Ifront);
#endif
#endif
  
}//end mode=1

}//end void

#endif



