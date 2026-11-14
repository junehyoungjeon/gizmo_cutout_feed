#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#ifdef FOFRAD
#include "../structure/fofrad.h"
extern fof_group_data *Groups;
#define VWVF_SCATTER 0.1  /* +/- this much around VariableWindVelFactor */
#endif

/*!
 *  routines for star formation in cosmological simulations
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for star formation and feedback have been
 *   replaced, and the algorithm is mostly new to GIZMO.
 */

#ifdef GALSF_SUBGRID_WINDS //RECOUPLE
static int nrec_dens=0,nrec_time=0,nrec_vel=0;
#endif

#ifdef GALSF // master switch for compiling the routines below //


#define MAXSFTEMP 5.5   // gas must have log(T)<MAXSFTEMP to form stars
#define WindInitialVelocityBoost 1.0 // (optional) boost velocity coupled (fixed momentum)

#if defined(GALSF_SFR_IMF_VARIATION) || defined(GALSF_SFR_IMF_SAMPLING)
/* function to determine what the IMF of a new star particle will be, based 
    on the gas properties of the particle out of which it forms */
void assign_imf_properties_from_starforming_gas(MyIDType i)
{
    
#ifdef GALSF_SFR_IMF_VARIATION
    double h = Get_Particle_Size(i) * All.cf_atime;
    double cs = Particle_effective_soundspeed_i(i) * All.cf_afac3; // actual sound speed in the simulation: might be unphysically high for SF conditions!
    cs = (1.9e4 / All.UnitVelocity_in_cm_per_s); // set to a minimum cooling temperature, for the actual star-forming conditions. for now, just use a constant //
    double dv2_abs = 0; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
    // squared norm of the trace-free symmetric [shear] component of the velocity gradient tensor //
    dv2_abs = ((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1])*(SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) +
                        (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2])*(SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) +
                        (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])*(SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
               (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] +
                         SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
                         SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) -
                        (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] +
                         SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
                         SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2]))) * All.cf_a2inv*All.cf_a2inv;
    double M_sonic = cs*cs*cs*cs / (All.G * dv2_abs * h);
    M_sonic *= All.UnitMass_in_g / All.HubbleParam / (1.989e33); // sonic mass in solar units //
    P[i].IMF_Mturnover = DMAX(0.01,DMIN(M_sonic,100.));
    P[i].IMF_Mturnover = 2.0; // 'normal' IMF in our definitions
    
    
    /* now we need to record all the properties we care to save about the star-forming gas, for the sake of later use: */
    int j,k;
    double NH = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1) * All.cf_a2inv;
    double dv2abs_tot = 0; /* calculate complete velocity dispersion (including hubble-flow correction) in physical units */
    for(j=0;j<3;j++)
    {
        for(k=0;k<3;k++)
        {
            double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
            if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
            dv2abs_tot += vt*vt;
        }
    }
    double acc=0,vel=0;
    for(k=0;k<3;k++)
    {
        double acc_tmp = P[i].GravAccel[k];
#ifdef PMGRID
        acc_tmp += P[i].GravPM[k];
#endif
        acc_tmp *= All.cf_a2inv;
        acc += acc_tmp * acc_tmp;
        vel += SphP[i].VelPred[k]*SphP[i].VelPred[k];
    }
    double b_mag = 0;
#ifdef MAGNETIC
    double gizmo2gauss = 4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam;
    for(k=0;k<3;k++) {b_mag += Get_Particle_BField(i,k)*Get_Particle_BField(i,k) * gizmo2gauss;}
#endif
    double rad_flux_uv = 1;
#ifdef GALSF_FB_LOCAL_UV_HEATING
    rad_flux_uv = SphP[i].RadFluxUV;
#endif
    double cr_energy_density = 0;
#ifdef COSMIC_RAYS
    cr_energy_density = SphP[i].CosmicRayEnergyPred * SphP[i].Density * All.cf_a3inv / P[i].Mass;
#endif

    P[i].IMF_FormProps[0] = P[i].IMF_Mturnover; // IMF turnover mass as defined above
    P[i].IMF_FormProps[1] = SphP[i].Density * All.cf_a3inv; // density
    P[i].IMF_FormProps[2] = SphP[i].InternalEnergyPred; // thermal internal energy (use to calculate temperature)
    P[i].IMF_FormProps[3] = Particle_effective_soundspeed_i(i) * All.cf_afac3; // sound speed (not trivially related to temperature if CRs, etc included)
    P[i].IMF_FormProps[4] = sqrt(dv2_abs); // shear velocity gradient (norm of shear gradient tensor)
    P[i].IMF_FormProps[5] = h; // particle length/size (inter-particle spacing)
    P[i].IMF_FormProps[6] = NH; // local gas surface density (our usual estimator) in the cloud where the particle formed
    P[i].IMF_FormProps[7] = sqrt(dv2abs_tot) * h; // total rms/turbulent velocity dispersion
    P[i].IMF_FormProps[8] = sqrt(acc); // gravitational acceleration
    P[i].IMF_FormProps[9] = sqrt(vel); // total velocity (use with acceleration to estimate shear omega, etc)
    P[i].IMF_FormProps[10] = sqrt(b_mag) * All.cf_a2inv; // magnetic field strength |B|
    P[i].IMF_FormProps[11] = rad_flux_uv; // incident UV flux normalized to MW 'canonical' (Habing) field value
    P[i].IMF_FormProps[12] = cr_energy_density; // cosmic ray energy density (if CRs are enabled)
    
#endif 
    
    
#ifdef GALSF_SFR_IMF_SAMPLING
    gsl_rng *random_generator_for_massivestars;
    random_generator_for_massivestars = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(random_generator_for_massivestars, P[i].ID+121);
    double mu = 0.01 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam / (1.989e33); // 1 O-star per 100 Msun
    unsigned int k = gsl_ran_poisson(random_generator_for_massivestars, mu);
    P[i].IMF_NumMassiveStars = (double)k;
#endif
    
}
#endif


/* return the light-to-mass ratio, for the IMF of a given particle, relative to the Chabrier/Kroupa IMF which 
    is otherwise (for all purposes) our 'default' choice */
double calculate_relative_light_to_mass_ratio_from_imf(MyIDType i)
{
#ifdef GALSF_SFR_IMF_VARIATION
    /* more accurate version from David Guszejnov's IMF calculations (ok for Mturnover in range 0.01-100) */
    double log_mimf = log10(P[i].IMF_Mturnover);
    return (0.051+0.042*(log_mimf+2)+0.031*(log_mimf+2)*(log_mimf+2)) / 0.31;
    // return pow(P[i].IMF_Mturnover/1.0,0.35);
#endif
#ifdef GALSF_SFR_IMF_SAMPLING
    double mu = 0.01 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam / (1.989e33); // 1 O-star per 100 Msun
    return P[i].IMF_NumMassiveStars / mu;
#endif
    return 1; // Chabrier or Kroupa IMF //
    // return 0.5; // Salpeter IMF down to 0.1 solar //
}


/* return the stellar age in Gyr for a given labeled age, needed throughout for stellar feedback */
double evaluate_stellar_age_Gyr(double stellar_tform)
{
    double age,a0,a1,a2,x0,x1,x2;
    if(All.ComovingIntegrationOn)
    {
        a0 = stellar_tform;
        a2 = All.Time;
        if(fabs(1-(All.Omega0+All.OmegaLambda))<=0.01)
        {
            /* use exact solution for flat universe */
            x0 = (All.Omega0/(1-All.Omega0))/(a0*a0*a0);
            x2 = (All.Omega0/(1-All.Omega0))/(a2*a2*a2);
            age = (2./(3.*sqrt(1-All.Omega0)))*log(sqrt(x0*x2)/((sqrt(1+x2)-1)*(sqrt(1+x0)+1)));
            age *= 1./All.Hubble_H0_CodeUnits;
        } else {
            /* use simple trap rule integration */
            a1 = 0.5*(a0+a2);
            x0 = 1./(a0*hubble_function(a0));
            x1 = 1./(a1*hubble_function(a1));
            x2 = 1./(a2*hubble_function(a2));
            age = (a2-a0)*(x0+4.*x1+x2)/6.;
        }
    } else {
        /* time variable is simple time, when not in comoving coordinates */
        age=All.Time-stellar_tform;
    }
    age *= 0.001*All.UnitTime_in_Megayears/All.HubbleParam; // convert to absolute Gyr
    return age;
}


/* return the (solar-scaled) light-to-mass ratio of an SSP with a given age; used throughout */
double evaluate_l_over_m_ssp(double stellar_age_in_gyr)
{
    // original SB99 tracks
    /*
    if(stellar_age_in_gyr < 0.0029)
    {
        return 1136.59;
    } else {
        double log_age = log10(stellar_age_in_gyr)-(-2.2681);
        return 478.63*pow(10.,-1.3625*log_age+0.115765*log_age*log_age);
        // could replace with piecewise linear functions; if this call in forcetree gets expensive //
    }
    */
    // updated SB99 tracks: including rotation, new mass-loss tracks, etc.
    if(stellar_age_in_gyr < 0.0035)
    {
        return 1136.59;
    } else {
        double log_age = log10(stellar_age_in_gyr/0.0035);
        return 1500.*pow(10.,-1.8*log_age+0.3*log_age*log_age-0.025*log_age*log_age*log_age);
    }
}



/* check whether conditions for star formation are fulfilled for a given particle */
int determine_sf_flag(int i)
{
    /*
     * f=1  normal cooling
     * f=0  star formation
     */
    
    double LogT=LogTemp(SphP[i].InternalEnergy*All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,SphP[i].Ne,i);
#if defined(POPIII_DISABLE_POPII) || defined(POPII_P2L)
     return 1;
#endif
     
    int flag = 1; /* default is normal cooling */ 
    //
    // this line modified by JDJ to add PopII metal threshold criteria to SF  Z/Zsun > 1e-4 and T<1e3 K
#ifndef POPIII_DISABLE_POPII
    if(SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh && P[i].Metallicity[0]/All.SolarAbundances[0] > POPII_METAL_THRESH && LogT<POPII_TEMP_THRESH) {flag = 0;}
#else
    if(SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh) {flag = 0;}
#endif
    
    if(All.ComovingIntegrationOn) {if(SphP[i].Density < All.OverDensThresh) flag = 1;}
    if(P[i].Mass <= 0) {flag = 1;}
#ifdef GALSF_JEANS_MIN_T
    if( LogT > 4.5 + 0.3333333 * log10(SphP[i].Density * All.cf_a3inv / All.MinRhoThresh) ) flag = 1;	 // This should be set to logTmaxcool in JeansMinT()
#else
    if(LogT >= MAXSFTEMP) flag = 1;
#endif

#ifdef GALSF_SUBGRID_WINDS
    if(SphP[i].DelayTime > 0)
        flag = 1; /* only normal cooling for particles in the wind */
#endif
    
    return flag;
}


/* simple routine to determine density thresholds and other common units for SF routines */
void set_units_sfr(void)
{
    All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
    All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
    
#ifdef GALSF_SFR_KMT
    All.UnitSurfaceDensity_in_Msun_per_pc2  = All.UnitDensity_in_cgs*All.UnitLength_in_cm*All.HubbleParam;  // g/cm^2
    All.UnitSurfaceDensity_in_Msun_per_pc2 /= SOLAR_MASS;           // Msun/cm^2
    All.UnitSurfaceDensity_in_Msun_per_pc2 *= pow(3.085678e18, 2);  // Msun/pc^2
#endif

#if defined(GALSF_EFFECTIVE_EQS) || defined(GALSF_SUBGRID_WINDS)
    double meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */
    All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
    All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
    meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
    All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
    All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif // GALSF_EFFECTIVE_EQS
}


/* Routine to actually determine the SFR assigned to an individual gas particle at each time */
double get_starformation_rate(int i)
{
    double rateOfSF,tsfr;
#ifdef GALSF_EFFECTIVE_EQS
    double factorEVP, egyhot, ne, tcool, x, y, cloudmass;
#endif
    // JDJ removed to make P2L work
    //if( determine_sf_flag(i) == 1 ) return 0;
    
#ifdef POPII_P2L
    tsfr = sqrt(POPII_SF_THRESH / (SphP[i].Density * All.cf_a3inv)) * All.MaxSfrTimescale;
#else
    tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * All.cf_a3inv)) * All.MaxSfrTimescale;
#endif
    
    if(tsfr<=0) return 0;
    rateOfSF = P[i].Mass / tsfr;
    //rateOfSF *= (1 - All.FactorSN);
    return rateOfSF;

#ifndef GALSF_EFFECTIVE_EQS
    /* 'normal' sfr from density law above */
    rateOfSF = P[i].Mass / tsfr;

#ifdef GALSF_INSTANTANEOUS_METALS
    rateOfSF *= (1 - All.FactorSN);
#endif
#else
    factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;
    //factorEVP = pow(SphP[i].Density * All.cf_a3inv / POPII_SF_THRESH, -0.8) * All.FactorEVP;
    egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
    ne = SphP[i].Ne;
    tcool = GetCoolingTime(egyhot, SphP[i].Density * All.cf_a3inv, &ne, i);
    y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
    x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
    cloudmass = x * P[i].Mass;
    rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;
    
    update_internalenergy_for_galsf_effective_eos(i,tcool,tsfr,x,rateOfSF); // updates entropies for the effective equation-of-state //
#endif // GALSF_EFFECTIVE_EQS
    
// #ifdef GALSF_SFR_KMT
//     double Sigma, Zmet, chi, s, fH2;
//     fH2    = 0.0;
//     Zmet   = P[i].Metallicity[0] / All.SolarAbundances[0];
// #ifdef GALSF_SFR_KMT_METALFLOOR
//     if( Zmet < GALSF_SFR_KMT_METALFLOOR ) Zmet = GALSF_SFR_KMT_METALFLOOR;
// #endif
//     // below is: Sigma = rho * HSML (must remove KMT #ifdef in hydra/gradients.c)
//     //Sigma  = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1);

//     // below is: Sigma = rho * rho / |gradRho|  (Sobolev)
//     Sigma  = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,0);

//     Sigma *= All.UnitSurfaceDensity_in_Msun_per_pc2 * All.cf_a2inv; // physical Msun/pc^2
//     chi    = 0.756 * (1. + 3.1 * pow(Zmet, 0.365));
//     s      = log(1. + 0.6*chi + 0.01*chi*chi) / (0.396 * Zmet * Sigma);
//     if( 0 < s && s < 2 )
//       fH2 = 1. - 0.75 * ( s / (1. + 0.25*s) );
//     rateOfSF *= fH2;

//     SphP[i].fH2   = fH2;
//     SphP[i].Sigma = Sigma / All.UnitSurfaceDensity_in_Msun_per_pc2 / All.cf_a2inv; // comoving system units 
// #endif

    
// #ifdef GALSF_SFR_MOLECULAR_CRITERION
//     /* Krumholz & Gnedin fitting function for f_H2 as a function of local properties */
//     double tau_fmol = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1) * All.cf_a2inv;
//     tau_fmol *= (0.1 + P[i].Metallicity[0]/All.SolarAbundances[0]);
//     if(tau_fmol>0) {
//         tau_fmol *= 434.78*All.UnitDensity_in_cgs*All.HubbleParam*All.UnitLength_in_cm;
//         y = 0.756*(1+3.1*pow(P[i].Metallicity[0]/All.SolarAbundances[0],0.365));
//         y = log(1+0.6*y+0.01*y*y)/(0.6*tau_fmol);
//         y = 1-0.75*y/(1+0.25*y);
//         if(y<0) y=0; if(y>1) y=1;
//         rateOfSF *= y;
//     } // if(tau_fmol>0)
// #endif // GALSF_SFR_MOLECULAR_CRITERION
    
    
// #ifdef GALSF_SFR_VIRIAL_SF_CRITERION
//     double dv2abs = 0; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
//     int j,k;
//     for(j=0;j<3;j++)
//         for(k=0;k<3;k++)
//         {
//             double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
//             if(All.ComovingIntegrationOn)
//                 if(j==k)
//                     vt += All.cf_hubble_a; /* add hubble-flow correction */
//             dv2abs += vt*vt;
//         }
//     /* add thermal support, although it is almost always irrelevant */
//     double cs_eff = Particle_effective_soundspeed_i(i);
//     double k_cs = cs_eff / (Get_Particle_Size(i)*All.cf_atime);
//     dv2abs += 2.*k_cs*k_cs; // account for thermal pressure with standard Jeans criterion (k^2*cs^2 vs 4pi*G*rho) //
    
//     //double alpha_vir = 0.2387 * dv2abs / (All.G * SphP[i].Density * All.cf_a3inv); // coefficient here was for old form, with only divv information
//     double alpha_vir = dv2abs / (8. * M_PI * All.G * SphP[i].Density * All.cf_a3inv); // 1/4 or 1/8 ? //

    
// #if !(EXPAND_PREPROCESSOR_(GALSF_SFR_VIRIAL_SF_CRITERION) == 1)
//      the above macro checks if GALSF_SFR_VIRIAL_SF_CRITERION has been assigned a numerical value 
// #if (GALSF_SFR_VIRIAL_SF_CRITERION > 0)
//     if(alpha_vir < 1.0)
//     {
//         /* check if Jeans mass is remotely close to solar; if not, dont allow it to form 'stars' */
//         double q = cs_eff * All.UnitVelocity_in_cm_per_s / (0.2e5);
//         double q2 = SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / (HYDROGEN_MASSFRAC*1.0e3*PROTONMASS);
//         double MJ_solar = 2.*q*q*q/sqrt(q2);
//         if(MJ_solar > 1000.) {alpha_vir = 100.;}
// #ifdef SINGLE_STAR_FORMATION
//         if(MJ_solar > 100.) {alpha_vir = 100.;}
// #endif
//     }
// #endif
// #if (GALSF_SFR_VIRIAL_SF_CRITERION > 1)
//     if(alpha_vir >= 1.0) {rateOfSF *= 0.0;}
// #endif
// #endif
    
//     if((alpha_vir<1.0)||(SphP[i].Density*All.cf_a3inv>100.*All.PhysDensThresh)) {rateOfSF *= 1.0;} else {rateOfSF *= 0.0015;}
//     // PFH: note the latter flag is an arbitrary choice currently set -by hand- to prevent runaway densities from this prescription! //
    
//     //  if( divv>=0 ) rateOfSF=0; // restrict to convergent flows (optional) //
//     //  rateOfSF *= 1.0/(1.0 + alpha_vir); // continuous cutoff w alpha_vir instead of sharp (optional) //
// #endif // GALSF_SFR_VIRIAL_SF_CRITERION
    
    return rateOfSF;
}



// #ifdef GALSF_EFFECTIVE_EQS
// /* compute the 'effective eos' cooling/heating, including thermal feedback sources, here */
// void update_internalenergy_for_galsf_effective_eos(int i, double tcool, double tsfr, double x, double rateOfSF)
// {
//     double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
//     double dtime = dt / All.cf_hubble_a; /*  the actual time-step */
//     double factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;
//     double trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
//     double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
//     double egyeff = egyhot * (1 - x) + All.EgySpecCold * x;
//     double egycurrent = SphP[i].InternalEnergy;

// #if defined(BH_THERMALFEEDBACK)
//     if((SphP[i].Injected_BH_Energy > 0) && (P[i].Mass>0))
//     {
//         egycurrent += SphP[i].Injected_BH_Energy / P[i].Mass;
//         if(egycurrent > egyeff)
//         {
//             tcool = GetCoolingTime(egycurrent, SphP[i].Density * All.cf_a3inv, &ne, i);
//             if(tcool < trelax && tcool > 0) trelax = tcool;
//         }
//         SphP[i].Injected_BH_Energy = 0;
//     }
// #endif // defined(BH_THERMALFEEDBACK)

//     /* now update the thermal variables */
//     SphP[i].InternalEnergy = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
//     SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
//     SphP[i].Pressure = get_pressure(i);
//     //SphP[i].dInternalEnergy = 0;
//     SphP[i].DtInternalEnergy = 0; /* HERE, it's ok, b/c effective EOS is designed to model new pressure even under compressions, 
//                                  (since we're zero'ing the second-half-step from the hydro step) */
// }
// #endif // GALSF_EFFECTIVE_EQS //




void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  unsigned int bits;
  double dt, dtime,dt_sec, mass_of_star, p, prob, rate_in_msunperyear, sfrrate, totsfrrate;
  double sum_sm, total_sm, sm=0, rate, sum_mass_stars, total_sum_mass_stars;
  double dt_popIII, sm_popIII, tot_sm_popIII, SFR_PopIII, SFR_PopIII_tot;
  double boxsize_mpc, p2_count, tot_p2_count;
  double old_sfrd, new_sfr,new_sfrd, cumulative_sfrd,inverse_years,sfr_time_average,sum_time,sum_popIII_mass,dt_years;
  int sfrd_ct,  popIII_flag, popIII_begin;
  int recouple();
#if defined(BH_POPIII_SEEDS) || defined(SINGLE_STAR_FORMATION) || defined(BH_BL)
  int num_bhformed=0, tot_bhformed=0;
  double GradRho;
#endif
#if defined(GALSF_FB_RPWIND_DO_IN_SFCALC) && defined(GALSF_FB_RPWIND_LOCAL)
  double total_n_wind,total_m_wind,total_mom_wind,total_prob_kick,avg_v_kick,momwt_avg_v_kick,avg_taufac;
  double totMPI_n_wind,totMPI_m_wind,totMPI_mom_wind,totMPI_prob_kick,totMPI_avg_v,totMPI_pwt_avg_v,totMPI_taufac;
  total_n_wind=total_m_wind=total_mom_wind=total_prob_kick=avg_v_kick=momwt_avg_v_kick=avg_taufac=0;
  totMPI_n_wind=totMPI_m_wind=totMPI_mom_wind=totMPI_prob_kick=totMPI_avg_v=totMPI_pwt_avg_v=totMPI_taufac=0;
#endif
    
#ifdef GALSF_SUBGRID_RECOUPLE
  int totnrec_dens,totnrec_time,totnrec_vel;
  nrec_dens=nrec_time=nrec_vel=0;
#endif
//#ifdef POPIII
  double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
//#endif
#ifdef BH_BL
  double BH1_count = 0;
  double BH2_count = 0;
  double tot_BH1_count=0;
  double tot_BH2_count=0;
  All.stellarBH=0.0;
#endif

#ifdef DENS_TERM
  double n_above_thres=0.0;
  double tot_n_above_thres=0.0;
#endif

double pSF, r, tFF;
if(All.TimeStep == 0){
  pSF = r = tFF = 0.0;
}
int SF_flag, coolflag;
double sfe;
  
#ifdef BH_DEBUG
  double max_p, max_tFF;
#endif
  
#ifdef GALSF_SUBGRID_WINDS
  double pvtau_return[4]; // not used
#endif
  
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinSfr[bin] = 0;
  new_sfr = 0.0;
  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;
  sm_popIII = tot_sm_popIII = 0;
  p2_count = tot_p2_count = 0;

  for(bits = 0; GALSF_GENERATIONS > (1 << bits); bits++);

//main loop
for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
      dtime = dt / All.cf_hubble_a; /*  the actual time-step */
      dt_sec = dtime * All.UnitTime_in_s / All.HubbleParam;  //physical seconds
      dt_years = dt_sec / SEC_PER_YEAR;
      //inverse_years = (All.cf_atime*All.cf_hubble_a)/(All.TimeStep*All.UnitTime_in_s/SEC_PER_YEAR);  // yr-1 per timestep

	coolflag = 1;
	flag =  1;  //default to cooling for particle (i.e. no PopII event)
	popIII_flag = 0; //default to no PopIII event  //default to no SF
	SF_flag = 0;

#ifdef PBH
    if((P[i].Type == 3) && (All.TimeStep == 0)) 
    {
       //P[i].BH_mass = All.mass_PBH;
       P[i].BH_mass = P[i].Mass / Msun_to_Internal;
       P[i].BH_type = 0;
       initBH(i);
    }
#endif

    if((P[i].Type == 0)&&(P[i].Mass>0)&&(dt>0))
    {
#ifdef GALSF_SUBGRID_WINDS
          if(SphP[i].DelayTime > 0) {
             SphP[i].DelayTime -= dtime;
#ifdef GALSF_SUBGRID_RECOUPLE
	     if( recouple(i) ) SphP[i].DelayTime = 0;
#endif
	   }
//#endif // GALSF_SUBGRID_WINDS
    
//#ifdef GALSF_SUBGRID_WINDS
    if(SphP[i].DelayTime > 0)
        coolflag = 1; /* only normal cooling for particles in the wind */
    else
#endif
	{
// #if defined(GALSF_EFFECTIVE_EQS)
//       if(flag==1)
// #endif
//         {
// #ifndef GALSF_TURNOFF_COOLING_WINDS
// #if defined(GALSF_SUBGRID_DS86WIND) || defined(GALSF_SUBGRID_HOTWIND)
// 	  if(SphP[i].DelayTime == 0)	/* if you have wind heating, don't cool wind particles */
// #endif
//             do_the_cooling_for_particle(i); /* actual cooling subroutine for particle i */
// #else
//             if(SphP[i].DelayTimeCoolingSNe<=0)
//             {
//                 do_the_cooling_for_particle(i); /* actual cooling subroutine for particle i */
//             } else {
//                 SphP[i].DelayTimeCoolingSNe -= dtime; /* 'counts down' until cooling is restored */
//             }
// #endif
//             SphP[i].Sfr = 0; /* will be reset below if flag==0 */
//         }

// do the cooling step first for all particle...then check for SF
// do_the_cooling_for_particle(i);
// SphP[i].Sfr = 0; /* will be reset below if flag==0 */

//PopIII SF criteria check: Metallicity < Metal_crit, Temp < 10^3, density > density_thresh, z < z_thresh
#ifdef DENS_TERM
    if(1./All.Time-1.0 < POPIII_z_THRESH){
      double dens_thres = All.density_threshold  * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
      if(SphP[i].Density * All.cf_a3inv >= dens_thres) n_above_thres += 1.0;
    }
#endif

#if defined(POPIII) && defined(GALSF) && defined(METALS)
      double logT = LogTemp(SphP[i].InternalEnergy*All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,SphP[i].Ne,i);
      double PopIII_dens_convert = POPIII_SF_THRESH  * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
      if(SphP[i].Density * All.cf_a3inv < PopIII_dens_convert) P[i].HighDens = 0;
      if(P[i].Metallicity[0]/All.SolarAbundances[0] <= POPIII_METAL_THRESH && 1./All.Time-1.0 < POPIII_z_THRESH)
      {
        if(SphP[i].Density * All.cf_a3inv >= PopIII_dens_convert && logT < POPIII_TEMP_THRESH)
        {
//#ifdef BH_BL
	        if(bits == 0)
            	number_of_stars_generated = 0;
	      	else
            	number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));

	      	mass_of_star = P[i].Mass / (GALSF_GENERATIONS - number_of_stars_generated);
            	if(number_of_stars_generated >= GALSF_GENERATIONS-1) mass_of_star=P[i].Mass;

#ifdef NSF_LIMIT_P3L
			if(number_of_stars_generated >= NSF_LIMIT_P3L)
			{
				pSF = 0.0;
				SphP[i].P3_heating_time = dt_years;//POPIII_PE_HEATING;
			}
			else
#endif
			{
				tFF = FreeFallTime(i);
				pSF = P[i].Mass / mass_of_star * (1 - exp(-dt_years/tFF * POPIII_SFE));// * P[i].Mass/Msun_to_Internal/POPIII_MASS;
			}
			//printf("%f: dt = %f Myr, tFF = %f Myr, p = %f\n", All.Time, dt_years/1e6, tFF/1e6, pSF);
			sum_sm += mass_of_star * pSF;
			r = get_random_number(P[i].ID + ThisTask + 1);
			if(r < pSF)
//#endif
			{
#if defined(GALSF_SUBGRID_WINDS) && defined(P3WIND)
//#ifndef NSFDBK
			assign_wind_kick_from_sf_routine(i, mass_of_star, dtime, pvtau_return);
            if (SphP[i].DelayTime > 0) coolflag = 1;
            else
//#endif
#endif
            {
            popIII_flag = 1;  //this keeps same particle from forming PopII in same timestep and used for type 
            P[i].P3L_flag = 1;
            flag = 0;
#ifdef ADCOOL
            if(number_of_stars_generated == (GALSF_GENERATIONS - 1)) coolflag = 0;
#endif
			P[i].parentID = P[i].ID;
            SF_flag = 1;
            sm_popIII += mass_of_star / Msun_to_Internal; 
            P[i].SF_dens = SphP[i].Density * All.cf_a3inv / PROTONMASS * (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
            P[i].H2_SF = SphP[i].Primordial_Chem[3];
#ifdef LYMAN_WERNER
#ifdef LW_LOCAL
            P[i].JLW_SF = P[i].JLW_local + LW_BACKGROUND(2e4, 1./All.Time-1.);
#else
			P[i].JLW_SF = LW_BACKGROUND(2e4, 1./All.Time-1.);
#endif
#endif
            P[i].StellarAge = All.Time; 
            P[i].Star_Type = 30000;
//#ifndef NSFDBK
#ifdef ADCOOL
			//P[i].SF_time += dt_years;
            P3L(i, mass_of_star/ Msun_to_Internal, 1, dt_years);
            //SphP[i].P3_heating_time = dt_years;
#endif
//#endif
			}
            }
        }
      }
#endif      

//  put P2L call here
#if defined(POPII_P2L) && defined(GALSF) && defined(METALS)
        double logT2 = LogTemp(SphP[i].InternalEnergy*All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,SphP[i].Ne,i);
        double PopII_dens_convert = POPII_SF_THRESH  * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
        
        //check SF criteria
        if(P[i].Metallicity[0]/All.SolarAbundances[0] >= POPII_METAL_THRESH 
            && SphP[i].Density * All.cf_a3inv >= PopII_dens_convert && logT2 < POPII_TEMP_THRESH && popIII_flag==0)
        {
//#ifdef BH_BL
	        if(bits == 0)
            	number_of_stars_generated = 0;
	      	else
            	number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));

	      	mass_of_star = P[i].Mass / (GALSF_GENERATIONS - number_of_stars_generated);
            	if(number_of_stars_generated >= GALSF_GENERATIONS-1) mass_of_star=P[i].Mass;
#ifdef NSF_LIMIT_P2L
			if(number_of_stars_generated >= NSF_LIMIT_P2L)// && P[i].P3L_flag<1)
			{	pSF = 0.0;
#ifndef NSFDBK
				SphP[i].P2_heating_time = dt_years;//POPII_PE_HEATING;
				//P[i].InHII = 1;
#endif
			}
			else
#endif
			{
			tFF = FreeFallTime(i);
			pSF = P[i].Mass / mass_of_star * (1 - exp(-dt_years/tFF * POPII_SFE));// * P[i].Mass*Msun_to_Internal/POPII_MASS;
			}
			sum_sm += mass_of_star * pSF;
			r = get_random_number(P[i].ID + ThisTask + 2);
			if(r < pSF)
//#endif
            {//Set P2 flag and PE_heating_time...call P2L mode=0 to identify heated particles
#if defined(GALSF_SUBGRID_WINDS) && defined(P2WIND)
#ifndef NSFDBK
            assign_wind_kick_from_sf_routine(i, mass_of_star, dtime, pvtau_return);
            if (SphP[i].DelayTime > 0) coolflag = 1;
            else
#endif
#endif
            {
            P[i].P2L_flag = 1;
            P[i].P3L_flag = 0;
            flag = 0;  
#ifdef ADCOOL
            if(number_of_stars_generated == (GALSF_GENERATIONS - 1)) coolflag = 0;
#endif
			P[i].parentID = P[i].ID;
            SF_flag = 1;
            p2_count += mass_of_star / Msun_to_Internal;
            P[i].SF_dens = SphP[i].Density * All.cf_a3inv / PROTONMASS * (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
            P[i].H2_SF = SphP[i].Primordial_Chem[3];
#ifdef LYMAN_WERNER
#ifdef LW_LOCAL
            P[i].JLW_SF = P[i].JLW_local + LW_BACKGROUND(2e4, 1./All.Time-1.);
#else
			P[i].JLW_SF = LW_BACKGROUND(2e4, 1./All.Time-1.);
#endif
#endif
            P[i].StellarAge = All.Time; //set age now so that SF data is correct.  This line commented out below
            P[i].Star_Type = 2;
//#ifndef NSFDBK
#ifdef ADCOOL
			//P[i].SF_time += dt_years;
            P2L(i,0, mass_of_star / Msun_to_Internal, 1, dt_years);
            //SphP[i].P3_heating_time = dt_years;
#endif
//#endif
			}
            //fflush(stdout);
            }
        }  
#endif

#ifdef BH_BL
#ifdef DCBH
      double BH_dens_convert = BH_THRESH  * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);

      // Check if particle is the most massive among its neighbors 
      int j,n,startnode,numngb,*Ngblist,dummy;
      int mostmassive = 0;
      MyFloat hsml;
    
      Ngblist = (int *) mymalloc("recouple neighbor list",NumPart * sizeof(int));
      hsml = PPP[i].Hsml;
      do {	// find neighbors,assign weights
        startnode = All.MaxPart;
        numngb = ngb_treefind_variable_threads(P[i].Pos,hsml,-1,&startnode,0,&dummy,&dummy,&dummy,Ngblist);
        for( n=0; n<numngb; n++ ) {
            j = Ngblist[n];
            if(P[j].Mass>P[i].Mass)
        {
            mostmassive = 1; // If not the most massive, set the flag
        } 
        }
      } while(startnode >= 0);
      myfree(Ngblist);
      // End of checking neighbors

      if(SF_flag==0 && SphP[i].Density * All.cf_a3inv >= BH_dens_convert && SphP[i].DelayTime <= 0)
      { 
      	double logT3 = LogTemp(SphP[i].InternalEnergy*All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,SphP[i].Ne,i);
	  	double T3 = pow(10, logT3);
      	if(T3 >= BH_TEMP_DOWN && T3 <= BH_TEMP_UP && P[i].Metallicity[0]/All.SolarAbundances[0] < BH_METAL_THRESH && P[i].HighDens==0 && mostmassive == 0)// && P[i].HighTemp==0)
      	{
#ifdef BH_H2_THRESH
      		if(SphP[i].Primordial_Chem[3] < BH_H2_THRESH)
#endif
      		{
      		printf("**************** Time to form a DCBH !!!  PID: %d\n",P[i].ID);
      		//P[i].Mass *= DCBH;
      		P[i].BH_mass = P[i].Mass / Msun_to_Internal*DCBH;
      		P[i].BH_type = 2;
      		P[i].H2_SF = SphP[i].Primordial_Chem[3];
#ifdef LYMAN_WERNER
#ifdef LW_LOCAL
            P[i].JLW_SF = P[i].JLW_local + LW_BACKGROUND(2e4, 1./All.Time-1.);
#else
			P[i].JLW_SF = LW_BACKGROUND(2e4, 1./All.Time-1.);
#endif
#endif
      		BH2_count += P[i].BH_mass;
      		num_bhformed += 1;
      		P[i].SF_dens = SphP[i].Density * All.cf_a3inv / PROTONMASS * (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
      		P[i].StellarAge = All.Time;
      		P[i].parentID = P[i].ID;
#ifdef BLACK_HOLES
			initBH(i);
#endif
      		flag = 0;
      		coolflag = 0;
      		}
      	}
      	/*
      	if(flag==1 && P[i].Metallicity[0]/All.SolarAbundances[0] >= POPII_METAL_THRESH && logT3 <= POPII_TEMP_THRESH)
      	{
      		sum_sm += mass_of_star;
            p2_count += mass_of_star/Msun_to_Internal;
            P[i].H2_SF = SphP[i].Primordial_Chem[3];
            P[i].JLW_SF = P[i].JLW_local + LW_BACKGROUND(2e4, 1./All.Time-1.);
            P2L(i,0, mass_of_star/Msun_to_Internal, 1); // mode=0 only determines particles in radius and sets PE heating time. 
            P[i].P2L_flag = 1;
            P[i].StellarAge = All.Time; //set age now so that SF data is correct.  This line commented out below
            P[i].Star_Type = 2;
            SF_flag = 1;
            flag = 0;
      	}
      	if(flag==1 && logT3 <= POPIII_TEMP_THRESH)
      	{
      		sum_sm += mass_of_star;
            sm_popIII += mass_of_star/Msun_to_Internal; 
            P[i].H2_SF = SphP[i].Primordial_Chem[3];
            P[i].JLW_SF = P[i].JLW_local + LW_BACKGROUND(2e4, 1./All.Time-1.);
            P3L(i, mass_of_star/Msun_to_Internal, 1);
            popIII_flag = 1;  //this keeps same particle from forming PopII in same timestep and used for type 
            P[i].P3L_flag = 1;
            P[i].StellarAge = All.Time; //set age now so that SF data is correct.  This line commented out below
            P[i].Star_Type += 30000;
            SF_flag = 1;
            popIII_flag = 1;
            flag = 0;
      	}
      	*/
#ifdef HII_FORBIDDEN
      	P[j].HighDens = 1;
#endif
      }
#endif
#endif
	}// if not a wind particle

#ifndef ADCOOL
	coolflag = flag;
#endif

if(coolflag == 1) {
    do_the_cooling_for_particle(i);
#ifdef ON_THE_FLY_FDBK
/*
#ifndef ADP_HII
if(SphP[i].DelayTime > 0)
#endif
{
    SphP[i].P2_heating_time = 0.0;
    SphP[i].P3_heating_time = 0.0;
#ifdef LW_OTF
    P[i].JLW_local = 0;
    P[i].LW_time = 0;
#endif
#ifdef UCHII_region
	P[i].InHII = 0;
#endif
}
*/
#endif  
/*#ifdef GALSF_SUBGRID_WINDS
    if(SphP[i].DelayTime > 0)
    {
    	SphP[i].P3_heating_time = 0.0;
    	SphP[i].P2_heating_time = 0.0;
    }
#endif*/
#ifdef JJ_profile
    CPU_Step[CPU_COOL] += measure_time();
#endif
    SphP[i].Sfr = 0; /* will be reset below if flag==0 */
}

/*
#ifdef BH_BL
if(BH_flag==1)
{
	P[i].Type = 3;
	TimeBinCountSph[P[i].TimeBin]--;
}
#endif
*/
        
if((flag == 0)&&(dt>0)&&(P[i].TimeBin))		/* active star formation (upon start-up, we need to protect against dt==0) */ 
    {
      if(SF_flag==1)
      {        
        /* 
	      if(bits == 0)
            number_of_stars_generated = 0;
	      else
            number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));

	      mass_of_star = P[i].Mass / (GALSF_GENERATIONS - number_of_stars_generated);
            if(number_of_stars_generated >= GALSF_GENERATIONS-1) mass_of_star=P[i].Mass;
		*/
		sfe = POPII_SFE;
		if(popIII_flag==1) sfe = POPIII_SFE;

          SphP[i].Sfr =  P[i].Mass * sfe/ tFF /Msun_to_Internal;
	      if(dtime>0) TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
/*
#ifdef POPII_P2L
              if(popIII_flag) 
              {
              	//P[NumPart + stars_spawned].Star_Type = 30000;  //add 30000 to indicate popIII star and conserve stellar population info
#ifdef ON_THE_FLY_FDBK
//#ifndef NSFDBK
#ifdef LW_OTF
				P3L(i, P[i].Mass/Msun_to_Internal, 1, dt_years);
#else
				P3L(i, P[i].Mass/Msun_to_Internal, 0, dt_years);
#endif
//#endif
#endif
              }
              else  
              {
              	//P[NumPart + stars_spawned].Star_Type = 2;
#ifdef ON_THE_FLY_FDBK
#ifdef LW_OTF
				P2L(i, 0, P[i].Mass/Msun_to_Internal, 1, dt_years);
#else
#ifndef NSFDBK
				P2L(i, 0, P[i].Mass/Msun_to_Internal, 0, dt_years);
#endif
#endif
#endif
              }
#endif
*/
if(number_of_stars_generated == (GALSF_GENERATIONS - 1))
	{
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;
		      sum_mass_stars += P[i].Mass;

		      P[i].Type = 4;
		      
		      TimeBinCountSph[P[i].TimeBin]--;
		      TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;

		      //P[i].StellarAge = All.Time;
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
                P[i].DensAroundStar = SphP[i].Density;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                P[i].Mass = SphP[i].MassTrue + SphP[i].dMass;
#endif                
    } /* closes final generation from original gas particle */
else
	{
		      /* here we spawn a new star particle */

		   if(NumPart + stars_spawned >= All.MaxPart)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPart);
			  fflush(stdout);
			  endrun(8888);
			}

		      P[NumPart + stars_spawned] = P[i];
		      P[NumPart + stars_spawned].Type = 4;
#ifdef POPII_P2L
              if(popIII_flag) 
              {
              	P[NumPart + stars_spawned].Star_Type = 30000;  //add 30000 to indicate popIII star and conserve stellar population info
              }
              else  
              {
              	P[NumPart + stars_spawned].Star_Type = 2;
              }
#endif
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
              P[NumPart + stars_spawned].DensAroundStar = SphP[i].Density;
#endif
		      NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
		      FirstActiveParticle = NumPart + stars_spawned;
		      NumForceUpdate++;

		      TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;
		      PrevInTimeBin[NumPart + stars_spawned] = i;
		      NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
		      if(NextInTimeBin[i] >= 0)
                  PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		      NextInTimeBin[i] = NumPart + stars_spawned;
		      if(LastInTimeBin[P[i].TimeBin] == i)
                  LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;

		      P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - bits));

		      P[NumPart + stars_spawned].Mass = mass_of_star;
		      P[i].Mass -= P[NumPart + stars_spawned].Mass;
              if(P[i].Mass<0) P[i].Mass=0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
              SphP[i].MassTrue -= P[NumPart + stars_spawned].Mass;
              if(SphP[i].MassTrue<0) SphP[i].MassTrue=0;
#endif
		      sum_mass_stars += P[NumPart + stars_spawned].Mass;
		      //P[NumPart + stars_spawned].StellarAge = All.Time;
		      force_add_star_to_tree(i, NumPart + stars_spawned);

		      stars_spawned++;
#ifdef POPIII
		      P[i].P3L_flag = 0;
#endif
#ifdef POPII_P2L
			  P[i].P2L_flag = 0;
#endif
	}//end else statement
	}//end of SF==1
#ifdef DCBH
	else
	{
		/* DCBH routines, developing */
		Stars_converted++;
		stars_converted++;
		P[i].Type = 3;
		TimeBinCountSph[P[i].TimeBin]--;
	}
#endif
	}//end of flag==0
    
#ifdef POPII_PE_HEATING
    if(SphP[i].P2_heating_time > 0.0) SphP[i].P2_heating_time -= dt_years; //counts down until PE heating turns off
#ifdef UCHII_region
	else
    P[i].InHII = 0;
#endif
#endif

#ifdef POPIII_PE_HEATING
    if(SphP[i].P3_heating_time > 0.0) SphP[i].P3_heating_time -= dt_years; //counts down until PE heating turns off
#ifdef UCHII_region
	else
    P[i].InHII = 0;
#endif
#endif

/*
#ifdef BH_BL
	double logT0 = LogTemp(SphP[i].InternalEnergy*All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,SphP[i].Ne,i);
	if(pow(10.0, logT0)>2*BH_TEMP_UP){
		P[i].HighTemp = 1;
	}
#endif
*/  
    } /* End of If Type = 0 */
    
    
    if((P[i].Type == 4)&&(P[i].Mass>0) && (SF_flag==0))
    {
#ifdef POPIII
        if(P[i].P3L_flag==1) {
        	if(P[i].SF_time >= POPIII_LIFE){
        	sm_popIII += Legacy_PopIII(i, P[i].Mass/Msun_to_Internal) - POPIII_MASS;
            P[i].P3L_flag = 0;
#ifdef BH_BL
			BH1_count += P[i].BH_mass;
			num_bhformed += All.stellarBH;
			All.stellarBH = 0;
			if(P[i].BH_mass > 0){
			P[i].Type = 3;}
#ifdef BH_SELF
			else{
			P[i].Mass = P[i].Mass * BH_SELF;}

#endif
			//P[i].StellarAge = All.Time;
#endif
#ifdef BLACK_HOLES
			initBH(i);
#ifdef EJECT_FRAC
			P[i].Mass *= (1 - EJECT_FRAC);
#endif
#ifdef BH_SELF
			if(BPP(i).BH_Mass>0)
			P[i].Mass = BPP(i).BH_Mass; //P[i].BH_mass * Msun_to_Internal;
#endif
#endif
#ifdef JJ_profile
            CPU_Step[CPU_P3L] += measure_time();
#endif
			}
#ifdef ON_THE_FLY_FDBK
			if(P[i].SF_time <= POPIII_PE_HEATING)
			//else
			{
//#ifndef NSFDBK
#ifdef LW_OTF
			P3L(i, P[i].Mass/Msun_to_Internal, 1, dt_years);
#else
			P3L(i, P[i].Mass/Msun_to_Internal, 0, dt_years);
#endif
//#endif
			}
#endif
			P[i].SF_time += dt_years;
        }
#endif

#ifdef POPII_P2L
        //once O+B stars die, then do metal + thermal enrichment (~10 Myr)
        if(P[i].P2L_flag==1) {
        	if(P[i].SF_time >= POPII_LIFE){
            P2L(i,1, P[i].Mass/Msun_to_Internal, 0, dt_years); //mode=1 does feedback.
            P[i].P2L_flag = 0; // turn off flag so it doesn't continue to enter this routine
#ifdef JJ_profile
            CPU_Step[CPU_P2L] += measure_time();
#endif
            }
#ifdef ON_THE_FLY_FDBK
			if(P[i].SF_time <= POPII_PE_HEATING)
			//else
			{
#ifdef LW_OTF
			P2L(i, 0, P[i].Mass/Msun_to_Internal, 1, dt_years);
#else
#ifndef NSFDBK
			P2L(i, 0, P[i].Mass/Msun_to_Internal, 0, dt_years);
#endif
#endif
			}
#endif
			P[i].SF_time += dt_years;
        }
#endif

    }
    
#ifdef LYMAN_WERNER
#ifdef LW_LOCAL
	if(P[i].LW_time > 0.0){
		P[i].LW_time -= dt_years;
	}
	else{
		P[i].JLW_local = 0.0;
	}
#endif
#endif
    
} /* end of main loop over active particles, huzzah! */

#ifdef DENS_TERM
All.reach_density_threshold=0.0;
if(1./All.Time-1.0 < POPIII_z_THRESH){
  MPI_Allreduce(&n_above_thres, &tot_n_above_thres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask==0){
    if(tot_n_above_thres>0)
      All.reach_density_threshold=1.0;
  }
      MPI_Bcast(&All.reach_density_threshold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
#endif

#if defined(GALSF_FB_RPWIND_DO_IN_SFCALC) && defined(GALSF_FB_RPWIND_LOCAL)
if(All.WindMomentumLoading)
{
    MPI_Reduce(&total_n_wind, &totMPI_n_wind, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_m_wind, &totMPI_m_wind, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_mom_wind, &totMPI_mom_wind, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_prob_kick, &totMPI_prob_kick, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_v_kick, &totMPI_avg_v, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&momwt_avg_v_kick, &totMPI_pwt_avg_v, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_taufac, &totMPI_taufac, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(totMPI_n_wind>0) {
            totMPI_avg_v /= totMPI_n_wind;
            totMPI_pwt_avg_v /= totMPI_prob_kick;
            totMPI_taufac /= totMPI_n_wind;
            totMPI_taufac = exp(totMPI_taufac);
            printf("Momentum Wind Feedback: Time=%g Nkicked=%g Mkicked=%g Momkicks=%g dPdtkick=%g V_avg=%g V_momwt_avg=%g Tau_avg=%g \n",
                   All.Time,totMPI_n_wind,totMPI_m_wind,totMPI_mom_wind,totMPI_prob_kick,totMPI_avg_v,
                   totMPI_pwt_avg_v,totMPI_taufac);
            fflush(stdout);
        }
        if(totMPI_n_wind>0) {
            fprintf(FdMomWinds, "%lg %g %g %g %g %g %g %g \n",
                    All.Time,totMPI_n_wind,totMPI_m_wind,totMPI_mom_wind,totMPI_prob_kick,totMPI_avg_v,
                    totMPI_pwt_avg_v,totMPI_taufac);
            fflush(FdMomWinds); 
        }
    }
}
#endif /* GALSF_FB_RPWIND_DO_IN_SFCALC */

    
#if defined(BH_POPIII_SEEDS) || defined(SINGLE_STAR_FORMATION) || defined(BH_BL) || defined(BLACK_HOLES)
  MPI_Allreduce(&num_bhformed, &tot_bhformed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_bhformed > 0)
  {
      if(ThisTask==0)
      {
      printf("BH/Sink formation: %d gas particles converted into BHs \n",tot_bhformed);
      fflush(stdout);
      }
      All.TotBHs += tot_bhformed;
  } // if(tot_bhformed > 0)
#endif

  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("SFR: spawned %d stars, converted %d gas particles into stars\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}
      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;
      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */
      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    } //(tot_spawned > 0 || tot_converted > 0)

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
        rate = total_sm / (All.TimeStep / (All.cf_atime*All.cf_hubble_a));
      else
        rate = 0;
      /* convert to solar masses per yr */
      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars);
      fflush(FdSfr);
    } // thistask==0
    CPU_Step[CPU_COOLINGSFR] += measure_time();

#if defined(POPIII) || defined(POPII_P2L)
    boxsize_mpc = All.BoxSize / All.HubbleParam / 1.e3;                                            // comoving Mpc
    sfr_time_average = SF_TIME;                                                                       // yr
    inverse_years = (All.cf_atime*All.cf_hubble_a)/(All.TimeStep*All.UnitTime_in_s/SEC_PER_YEAR);  // yr^-1 per timestep
#endif

/*
#ifdef BH_DEBUG
if(SF_flag==1){
	MPI_Allreduce(&pSF, &max_p, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&tFF, &max_tFF, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(ThisTask==0){
	fprintf(STSF, "%f: dt = %f Myr, tFF = %f Myr, p = \%f\n", All.Time, dt_years/1e6, max_tFF/1e6, max_p);
	}
}
#endif
*/

#ifdef BH_BL
	MPI_Allreduce(&BH1_count, &tot_BH1_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&BH2_count, &tot_BH2_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(ThisTask==0){
		if(All.TimeStep==0){
            fprintf(FdBH1fr, "%g %g %g %g\n", All.Time, All.BH1_TIME, All.BH1_FR, All.BH1_FRD);
            fflush(FdBH1fr);
            fprintf(FdBH2fr, "%g %g %g %g\n", All.Time, All.BH2_TIME, All.BH2_FR, All.BH2_FRD);
            fflush(FdBH2fr);
		}
		
		All.BH1_SUM += tot_BH1_count;
		All.BH1_TIME += (1./inverse_years);
		if(All.BH1_TIME > sfr_time_average){
			new_sfr = All.BH1_SUM / All.BH1_TIME;
			All.BH1_FR = new_sfr;
			All.BH1_FRD = new_sfr / (boxsize_mpc * boxsize_mpc * boxsize_mpc * All.ZoomInFactor); 
			//printf("\nwhat???%f\n", All.reach_density_threshold); 
			MPI_Bcast(&All.BH1_FRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            fprintf(FdBH1fr, "%g %g %g %g\n", All.Time, All.BH1_TIME, All.BH1_FR, All.BH1_FRD);
            fflush(FdBH1fr);
            All.BH1_TIME = All.BH1_TIME - sfr_time_average;
            All.BH1_SUM = 0.0;
		}
		
		All.BH2_SUM += tot_BH2_count;
		All.BH2_TIME += (1./inverse_years);
		if(All.BH2_TIME > sfr_time_average){
			new_sfr = All.BH2_SUM / All.BH2_TIME;
			All.BH2_FR = new_sfr;
			All.BH2_FRD = new_sfr / (boxsize_mpc * boxsize_mpc * boxsize_mpc * All.ZoomInFactor);  
			MPI_Bcast(&All.BH2_FRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            fprintf(FdBH2fr, "%g %g %g %g\n", All.Time, All.BH2_TIME, All.BH2_FR, All.BH2_FRD);
            fflush(FdBH2fr);
            All.BH2_TIME = All.BH2_TIME - sfr_time_average;
            All.BH2_SUM = 0.0;
		}
	}
	if(ThisTask>0){
		All.BH1_TIME += (1./inverse_years);
		if(All.BH1_TIME > sfr_time_average){
			MPI_Bcast(&All.BH1_FRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			All.BH1_TIME = All.BH1_TIME - sfr_time_average;
		}
		All.BH2_TIME += (1./inverse_years);
		if(All.BH2_TIME > sfr_time_average){
			MPI_Bcast(&All.BH2_FRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			All.BH2_TIME = All.BH2_TIME - sfr_time_average;
		}
	}
#endif

#ifdef POPIII
  /*                      determine Pop III star formation rate and sfrd over sfr_time_average                  */
    
    MPI_Allreduce(&sm_popIII, &tot_sm_popIII, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);             // popIII mass over timestep
    if(ThisTask==0){
        if(All.TimeStep == 0){
            fprintf(FdPopIIIsfr, "%g %g %g %g\n", All.Time, All.POPIII_TIME, All.POPIII_SFR, All.POPIII_SFRD);
            fflush(FdPopIIIsfr);
        }

        All.POPIII_SUM += DMAX(tot_sm_popIII, 0.0); // total mass over time average
        All.POPIII_TIME += (1./inverse_years);

        if(All.POPIII_TIME > sfr_time_average){
            new_sfr = All.POPIII_SUM / All.POPIII_TIME;
            All.POPIII_SFR = new_sfr;
            All.POPIII_SFRD = new_sfr / (boxsize_mpc * boxsize_mpc * boxsize_mpc * All.ZoomInFactor);                // Msun/yr/mpc^3
            MPI_Bcast(&All.POPIII_SFRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*  here I am not including little h in the SFR since it cancels between mass and time (not included in either);pyth
        little h is however included in the sfrd so I am outputting comoving Mpc^3 */
        
            fprintf(FdPopIIIsfr, "%g %g %g %g\n", All.Time, All.POPIII_TIME, All.POPIII_SFR, All.POPIII_SFRD);
            fflush(FdPopIIIsfr);
            All.POPIII_TIME -= sfr_time_average;
            All.POPIII_SUM = 0.0;
        } 
    }//end ThisTask==0
    
	if(ThisTask>0){
		All.POPIII_TIME += (1./inverse_years);
		if(All.POPIII_TIME > sfr_time_average){
			MPI_Bcast(&All.POPIII_SFRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			All.POPIII_TIME -= sfr_time_average;
		}
	}

#endif

#ifdef POPII_P2L
  /*                      determine Pop II star formation rate and sfrd over sfr_time_average                  */
    MPI_Allreduce(&p2_count, &tot_p2_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
               // popII mass over timestep
    if(ThisTask==0){
        if(All.TimeStep == 0){
            fprintf(FdPopIIsfr, "%g %g %g %g\n", All.Time, All.POPII_TIME, All.POPII_SFR, All.POPII_SFRD);
            fflush(FdPopIIsfr);
        }

        All.POPII_SUM += (tot_p2_count); //*POPII_MASS);                                // total mass over time average
        All.POPII_TIME += (1./inverse_years);

        if(All.POPII_TIME > sfr_time_average){
            new_sfr = All.POPII_SUM / All.POPII_TIME;
            All.POPII_SFR = new_sfr;
            All.POPII_SFRD = new_sfr / (boxsize_mpc * boxsize_mpc * boxsize_mpc * All.ZoomInFactor);                // Msun/yr/mpc^3
			MPI_Bcast(&All.POPII_SFRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*  here I am not including little h in the SFR since it cancels between mass and time (not included in either);
        little h is however included in the sfrd so I am outputting comoving Mpc^3 */
        
            fprintf(FdPopIIsfr, "%g %g %g %g\n", All.Time, All.POPII_TIME, All.POPII_SFR, All.POPII_SFRD);
            fflush(FdPopIIsfr);
            All.POPII_TIME -= sfr_time_average;
            All.POPII_SUM = 0.0;
        }

    }//end ThisTask==0
  
	if(ThisTask>0){
		All.POPII_TIME += (1./inverse_years);
		if(All.POPII_TIME > sfr_time_average){
			MPI_Bcast(&All.POPII_SFRD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			All.POPII_TIME -= sfr_time_average;
		}
	}

#endif




#ifdef GALSF_SUBGRID_RECOUPLE
  MPI_Allreduce(&nrec_time, &totnrec_time, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nrec_dens, &totnrec_dens, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nrec_vel, &totnrec_vel, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if( ThisTask==0 && totnrec_time+totnrec_dens+totnrec_vel>0 ) printf("RECOUPLE_STATS: z= %g ntime= %d %g  ndens= %d %g  nvel= %d %g\n",1./All.cf_atime-1,totnrec_time,1.*totnrec_time/(totnrec_time+totnrec_dens+totnrec_vel),totnrec_dens,1.*totnrec_dens/(totnrec_time+totnrec_dens+totnrec_vel),totnrec_vel,1.*totnrec_vel/(totnrec_time+totnrec_dens+totnrec_vel));
#endif
#ifdef JJ_profile
  CPU_Step[CPU_SFR] += measure_time();
#endif

} /* end of main sfr_cooling routine!!! */



#ifdef GALSF_SUBGRID_WINDS
//#ifdef GALSF_SUBGRID_RECOUPLE
int recouple(int i)	// Decide whether to recouple hydro for wind particle
{
    int flag=0;
#ifdef GALSF_SUBGRID_RECOUPLE
    int j,k,n,nd,startnode,numngb,*Ngblist,dummy;
    MyFloat hsml,sumwk,r2,h2,u,hinv3,hinv4,wk,dwk,dv[3],dx[3],dv2,du,norm_j,csound2,u_amb,*Ngbweight;
#endif

    if(SphP[i].DelayTime < 0) {
	nrec_time++;
	return 1;
    }
    if(SphP[i].Density * All.cf_a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh) {
	nrec_dens++;
	return 1;
    }

/* Do recoupling based on velocity relative to local sound speed */
#ifdef GALSF_SUBGRID_RECOUPLE
    if( SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh || P[i].ID > All.TotNumPart ) return 0; // don't recouple in ISM, and don't bother with non-local particles

/*    double densconv = All.HubbleParam*All.HubbleParam * All.cf_a3inv*HYDROGEN_MASSFRAC*All.UnitMass_in_g /
        (All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm);        // converts to H density
    double rhocrit = 1.88e-29 * All.HubbleParam * All.HubbleParam * All.cf_a3inv * All.OmegaBaryon * HYDROGEN_MASSFRAC; // physical critical density 
    if( SphP[i].Density*densconv < 200*rhocrit ) {
	nrec_dens++;
	return 1;   // recouple at typical halo density (should almost never get here)
    }*/
    Ngblist = (int *) mymalloc("recouple neighbor list",NumPart * sizeof(int));
    Ngbweight = (MyFloat *) mymalloc("recouple neighbor weights",NumPart * sizeof(MyFloat));
    hsml = PPP[i].Hsml;
    do {	// find neighbors,assign weights
        startnode = All.MaxPart;
        numngb = ngb_treefind_variable_threads(P[i].Pos,hsml,-1,&startnode,0,&dummy,&dummy,&dummy,Ngblist);
        for( n=0, nd=0, sumwk=0; n<numngb; n++ ) {
            j = Ngblist[n];
            Ngbweight[n] = 0;
	    if ( i == j || SphP[j].DelayTime > 0) {
                nd++;       // ignore wind particles for computing ambient quantities
                continue;
            }
	    for( k=0; k<3; k++ ) dx[k] = P[i].Pos[k] - P[j].Pos[k];
#ifdef PERIODIC
	    NEAREST_XYZ(dx[0],dx[1],dx[2],1);
#endif
            r2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
            h2 = hsml*hsml;
            if( r2 < h2 ) {
                u = sqrt(r2/h2);
                hinv4 = 1/(h2*h2);
                hinv3 = hinv4*hsml;
                kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
                sumwk += wk;
                Ngbweight[n] = wk;
            }
            else nd++;      /* particle is too far away; ignore */
         }
	 numngb -= nd;
	 hsml *= 1.05;	// in case hsml not large enough, increase till we get enough neighbors
    } while( numngb < All.DesNumNgb );

    for( k=0; k<3; k++ ) dv[k] = 0;
    for( n=0,du=0,u_amb=0; n<numngb+nd; n++ ) {
        j = Ngblist[n];
	if( i == j ) continue;
        norm_j = Ngbweight[n]/(sumwk+1.e-20);
	for( k=0; k<3; k++ ) dv[k] += norm_j*(P[i].Vel[k] - P[j].Vel[k]);
	du += norm_j * (SphP[i].InternalEnergy-SphP[j].InternalEnergy);
	u_amb += norm_j * SphP[j].InternalEnergy;
    }
    dv2 = (dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]) * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s * All.cf_a2inv;
    csound2 = GAMMA*GAMMA_MINUS1 * u_amb * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs);
/* once outside the ISM, recouple when v is small compared to sound speed */
    if( dv2 < csound2*GALSF_SUBGRID_RECOUPLE*GALSF_SUBGRID_RECOUPLE ) {
#if defined(GALSF_SUBGRID_HOTWIND) || defined(GALSF_SUBGRID_DS86WIND)
	if( SphP[i].InternalEnergy < u_amb ) 	// when winds are hot, only heat more at recoupling
#endif
	SphP[i].InternalEnergy = u_amb;	// set particle temperature to that of surrounding gas
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);
	flag = 1; 
    }
    myfree(Ngbweight);
    myfree(Ngblist);
	if( flag == 1 ) nrec_vel++;
	if( flag == 1 ) printf("RECOUPLE: P%d %d Mh=%.3g nH=%.3g Ta=%.3g dv=%g cs=%g du/u= %.3g\n",ThisTask,P[i].ID,log10(SphP[i].HostHaloMass*All.UnitMass_in_g / (SOLAR_MASS*All.HubbleParam) +1.e-20),log10(SphP[i].Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * HYDROGEN_MASSFRAC / PROTONMASS * All.cf_a3inv),LogTemp(u_amb * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs),SphP[i].Ne,i),sqrt(dv2)*1.e-5,sqrt(csound2)*1.e-5,du/SphP[i].InternalEnergy);
#endif

    return flag;	// if no criteria are true, don't recouple
}
//#endif
#endif


#if defined(GALSF_FB_RPWIND_DO_IN_SFCALC) || defined(GALSF_SUBGRID_WINDS) || defined(GALSF_SUBGRID_VARIABLEVELOCITY) || defined(GALSF_SUBGRID_VARIABLEVELOCITY_DM_DISPSERSION)
void assign_wind_kick_from_sf_routine(int i, double sm, double dtime, double pvtau_return[4])
{
    int j;
    double v,p,prob;
    double norm, dir[3];
#ifdef GALSF_FB_RPWIND_LOCAL
    double tau_IR=0;
    double m_gas_kernel,h,m_gas_kernel_i,vq;
    double h_kernel_i;
    double m_st_kernel,m_st_kernel_i,l_st_kernel,rho_to_launch;
#ifdef GALSF_FB_RPWIND_FROMSTARS
    MyDouble *pos;
    double star_age,log_age,lm_ssp;
    double unittime_in_gyr = 0.001*All.UnitTime_in_Megayears/All.HubbleParam;
    int startnode,numngb_inbox,n,dummy;
#endif
#endif
#ifdef GALSF_FB_RPWIND_FROMCLUMPS
    int ngb_run_cntr,kmin,k;
    double dmin1w,dmax1w,dmax2w;
#endif
#if defined(GALSF_WINDS_ISOTROPIC) || defined(GALSF_WINDS_RADIAL)
    double theta=0, phi;
#endif

#ifdef GALSF_SUBGRID_WINDS
    /* this is the simple, old standard wind model, with constant velocity & loading with SFR */
    p = All.WindEfficiency * sm / P[i].Mass;
    v = sqrt(2 * All.WindEnergyFraction*All.FactorSN*All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency);
    prob = 1 - exp(-p);
#endif
    

#ifdef GALSF_SUBGRID_DMDISPERSION
    /* wind model where launching scales with halo/galaxy bulk properties (as in Vogelsberger's simulations) */
    if(SphP[i].DM_VelDisp > 0 && sm > 0)
    {
        double wind_energy, wind_momentum, wind_mass;
        v = All.VariableWindVelFactor * SphP[i].DM_VelDisp;  /* physical wind velocity */
        //      if(v < 50.0) v = 50.0;
        wind_momentum = sm * All.VariableWindSpecMomentum;
        wind_energy = sm * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN);
        wind_mass = (wind_energy + sqrt(wind_energy * wind_energy + v * v * wind_momentum * wind_momentum)) / (v * v);
        /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
        p = wind_mass / P[i].Mass;
    }
    else
    {
        v = 0;
        p = 0;
    }
    prob = 1 - exp(-p);
#endif
    
    
    
#ifdef GALSF_SUBGRID_VARIABLEVELOCITY
    /* wind model where launching scales with halo/galaxy bulk properties (as in Romeel's simulations) */
    int canLaunch = 0;
    double v_circ=0, wind_mass;
#ifdef FOFRAD
    if(sm > 0 )
      canLaunch = 1;
#else
    if(SphP[i].HostHaloMass > 0 && sm > 0)
      canLaunch = 1;
#endif	// FOFRAD

    if(canLaunch)
    {
#ifdef FOFRAD
        // Baryonic Tully-Fisher relation from McGaugh 2012, plus z-scaling from MMW98 (PHYSICAL km/s)
        v_circ = (pow(Groups[P[i].galID].Mbar/102.329, 0.26178) * 
		  pow(hubble_function(All.Time)/All.Hubble_H0_CodeUnits, 0.3333));

	/* PHYSICAL wind velocity +/- VWVF_SCATTER */
	v = ((All.VariableWindVelFactor - VWVF_SCATTER) + (2.0 * VWVF_SCATTER * get_random_number(P[i].ID))) * v_circ;
#ifdef GALSF_SUBGRID_FIREVEL
#define GALSF_SUBGRID_FIREVEL_SLOPE 0.12	// Muratov+15: 0.12
	v = GALSF_SUBGRID_FIREVEL*pow(v_circ*All.cf_atime/200,GALSF_SUBGRID_FIREVEL_SLOPE)*pow(200,0.12-GALSF_SUBGRID_FIREVEL_SLOPE)*(1.-VWVF_SCATTER+2.*VWVF_SCATTER*get_random_number(P[i].ID)) * v_circ;	// eq.9 in Muratov+15, with tophat scatter
#endif
#ifdef GALSF_SUBGRID_VBOOST
	v = wind_vel_boost(i,v_circ,v);
#endif
#else	// FOFRAD
        double HaloConcentrationNorm = 9.;  /* concentration c0 of a halo of unit mass */
        double HaloConcentrationSlope = -0.15;  /* slope n of mass concentration relation, namely c = c0 * M_200,crit^n */

        double r200c, v_esc, c_halo, wind_energy, wind_momentum, wind_mass;
        double rhocrit = 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
        rhocrit *= All.Omega0/All.cf_a3inv + (1-All.Omega0-All.OmegaLambda)/All.cf_a2inv + All.OmegaLambda; /* physical critical density at redshift z */

        r200c  = pow(SphP[i].HostHaloMass / (4 * M_PI / 3.0 * 200 * rhocrit), 1.0 / 3.0);	/* physical r_200,crit value, assuming FoF mass = M_200,crit */
        v_circ = sqrt(All.G * SphP[i].HostHaloMass / r200c);	/* physical circular velocity at r_200,crit */
        c_halo = HaloConcentrationNorm * pow(SphP[i].HostHaloMass, HaloConcentrationSlope);
        v_esc  = v_circ * sqrt(2 * c_halo / (log(1 + c_halo) - c_halo / (1 + c_halo)));	/* physical escape velocity of halo */
        v = All.VariableWindVelFactor * v_esc;	/* physical wind velocity */

        wind_momentum = sm * All.VariableWindSpecMomentum;
        wind_energy = sm * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN);

        wind_mass = (wind_energy + sqrt(wind_energy * wind_energy + v * v * wind_momentum * wind_momentum)) / (v * v);
#endif 	// FOFRAD
 

#ifdef GALSF_SUBGRID_MDWIND
        /* Similar to Oppenheimer's momentum-driven winds; v_circ should be in PHYSICAL km/s */
        wind_mass = sm * GALSF_SUBGRID_MDWIND / v_circ;

#if defined(GALSF_SUBGRID_MDSLOPE)
#if defined(GALSF_SUBGRID_MDBREAK) 
        /* general form for having a different slope for eta for v_circ<GALSF_SUBGRID_MDBREAK.
 *            The slope below is the *negative* of the power of v_circ below the break; 
 *                       for instance, the energy-driven winds would set GALSF_SUBGRID_MDSLOPE=2 */
        if( v_circ < GALSF_SUBGRID_MDBREAK ) wind_mass *= pow(GALSF_SUBGRID_MDBREAK / v_circ, GALSF_SUBGRID_MDSLOPE-1.);
#else
        /* general form for having a different slope everywhere */
        wind_mass = sm * pow(GALSF_SUBGRID_MDWIND / v_circ,GALSF_SUBGRID_MDSLOPE);
#endif
#endif

#if defined(GALSF_SUBGRID_ETAZPOWER) && defined(GALSF_SUBGRID_ETAZMAX)
	//Muratov 15 (arXiv: 1501.03155) eta z-scalings
	//z_max     = 4.0;
	//eta_power = 1.25;
	//
	//Our adjustment
	//z_max     = 2.0;
	//eta_power = 1.75;

	double z_eta;
	z_eta = 1.0 / All.Time - 1.0;

	if( z_eta > GALSF_SUBGRID_ETAZMAX )
	  z_eta = GALSF_SUBGRID_ETAZMAX;

	wind_mass *= pow(1.0 + z_eta, GALSF_SUBGRID_ETAZPOWER);
#endif
#endif //GALSF_SUBGRID_MDWIND

#ifdef GALSF_SUBGRID_MSWIND
        double Mstar=Groups[P[i].galID].Mstar;
        double massconv = All.UnitMass_in_g/(1.989e33*All.HubbleParam); // converts to solar masses

        if( Mstar < All.OrigGasMass*0.5*FOF_MINGROUPLEN*massconv ) Mstar = All.OrigGasMass*0.5*FOF_MINGROUPLEN*massconv;  // set minimum Mstar for use in mass loading formula
        wind_mass = sm * GALSF_SUBGRID_MSNORM * pow(Mstar*1.e-10,GALSF_SUBGRID_MSWIND); // form based on Muratov+15 eq. 8:  MSWIND=-0.351 and MSNORM=3.55
#ifdef GALSF_SUBGRID_ETAHISLOPE
	if( Mstar > 1.e10 ) wind_mass *= pow(Mstar*1.e-10,GALSF_SUBGRID_ETAHISLOPE);
#endif
	//if( Groups[P[i].galID].Mstar < All.OrigGasMass*FOF_MINGROUPLEN ) wind_mass *= Groups[P[i].galID].Mstar/(All.OrigGasMass*FOF_MINGROUPLEN);  // attenuate winds before galaxy is fully resolved
#endif

#ifdef SIGMASFR_WIND
	int gID=P[i].galID;
	double r2;
	double fofrad_periodic();
	r2 = fofrad_periodic(P[i].Pos[0]-Groups[gID].COM[0])*fofrad_periodic(P[i].Pos[0]-Groups[gID].COM[0])
		+ fofrad_periodic(P[i].Pos[1]-Groups[gID].COM[1])*fofrad_periodic(P[i].Pos[1]-Groups[gID].COM[1])
		+ fofrad_periodic(P[i].Pos[2]-Groups[gID].COM[2])*fofrad_periodic(P[i].Pos[2]-Groups[gID].COM[2]);
	if( gID == 0 ) wind_mass = 0;	// particles not in a galaxy can't join wind
	else if( r2 > Groups[gID].rw2 ) wind_mass *= Groups[gID].rw2/r2;
	//if( ThisTask==0 && gID>0 ) printf("SIGMASFR_WIND: %g %d %g %g %g %g %g\n",All.cf_atime,gID,Groups[gID].Mstar,Groups[gID].Mstar/All.OrigGasMass,Groups[gID].Sfr/(M_PI*Groups[gID].rw2*All.cf_atime*All.cf_atime+1.e-30),sqrt(r2),sqrt(Groups[gID].rw2));
#endif


#ifdef OUTPUT_SUBGRID_WINDINFO
	FILE *FdWindFile;
	char buf[200];
	sprintf(buf, "%sSubgridWinds", All.OutputDir);
	mkdir(buf, 02755);
	sprintf(buf, "%sSubgridWinds/subgridWinds.%d.bin", All.OutputDir, ThisTask);
	if(!(FdWindFile = fopen(buf, "a"))){
	  printf("error in opening file '%s'\n",buf);
	  endrun(1);
	}
#ifdef FOFRAD
	double buffer;
	buffer = All.Time;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = Groups[P[i].galID].Mbar;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = Groups[P[i].galID].Mstar;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = v_circ;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = v;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = wind_mass / sm;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = (double)SphP[i].Density;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = LogTemp(SphP[i].InternalEnergy * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs), SphP[i].Ne,i);
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = (double)SphP[i].Ne;
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = (double)P[i].Metallicity[0];
	fwrite(&buffer, sizeof(double), 1, FdWindFile);
	buffer = (double)(SphP[i].N_windLaunches + 1);
	fwrite(&buffer, sizeof(double), 1, FdWindFile);

	/*
	fprintf(FdWindFile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",
		All.Time,Groups[P[i].galID].Mbar,Groups[P[i].galID].Mstar,
		v_circ,v*All.cf_atime,wind_mass/sm,
		SphP[i].Density,LogTemp(SphP[i].InternalEnergy * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs), SphP[i].Ne),
		SphP[i].Ne,P[i].Metallicity[0],SphP[i].N_windLaunches+1);
	*/
	/*
	if(ThisTask==0)
	  printf("FOFRADWIND a= %g Mbar= %g Mstar= %g Eta= %g vcirc= %g\n",
		 All.Time,Groups[P[i].galID].Mbar,Groups[P[i].galID].Mstar,wind_mass/sm,v_circ);
	*/
#else  // FOFRAD
	fprintf(FdWindFile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",
		All.Time,v*All.cf_atime,wind_mass/sm,SphP[i].Density,
		LogTemp(SphP[i].InternalEnergy * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs), SphP[i].Ne,i),SphP[i].Ne,
		P[i].Metallicity[0],SphP[i].N_windLaunches+1);
#endif // FOFRAD
	fflush(FdWindFile);
	fclose(FdWindFile);
#endif // output_subgrid_windinfo

        /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
        p = wind_mass / P[i].Mass;
    }
    else
    {
        v = 0;
        p = 0;
    }
    prob = 1 - exp(-p);
#endif
    
    
    
#ifdef GALSF_FB_RPWIND_LOCAL
    /* revised winds model (blowing winds out from local dense clumps) */
    double unitmass_in_e10solar = All.UnitMass_in_g/SOLAR_MASS/1.0e10/All.HubbleParam;
    double unitlength_in_kpc = All.UnitLength_in_cm/3.086e21/All.HubbleParam;
    double unitrho_in_e10solar_kpc3 = unitmass_in_e10solar/(unitlength_in_kpc*unitlength_in_kpc*unitlength_in_kpc);
    double unitvel_in_km_s = All.UnitVelocity_in_cm_per_s/1.0e5;
    if(All.WindMomentumLoading<=0)
    {
        p=v=0;
    } else {
        m_st_kernel=l_st_kernel=m_st_kernel_i=h_kernel_i=m_gas_kernel_i=m_gas_kernel=0;
        
#if !defined(GALSF_FB_RPWIND_FROMSTARS)
        /* we only weakly revise the model here, to scale velocities with density estimate */
        m_gas_kernel=0; h=PPP[i].Hsml;
        rho_to_launch = SphP[i].Density*All.cf_a3inv;
        v = WindInitialVelocityBoost * (45.0/unitvel_in_km_s) * pow(rho_to_launch*unitrho_in_e10solar_kpc3, 0.25);
        if(v<0.01/unitvel_in_km_s) v=0.01/unitvel_in_km_s;
        p = 116 * sm / (P[i].Mass*(v/unitvel_in_km_s));
#endif
        
#ifdef GALSF_FB_RPWIND_FROMSTARS
        /* here we have to do a loop to search for nearby stars */
        Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
        startnode = All.MaxPart;
        dummy=0;h=0;numngb_inbox=0;
        h=3.0*PPP[i].Hsml;
        pos=P[i].Pos;
        m_st_kernel=0; l_st_kernel=0; //l_st_kernel_nonrad=0;
        do {
            numngb_inbox = ngb_treefind_newstars(&pos[0],h,-1,&startnode,0,&dummy,&dummy);
            /* searches for all new stars inside h */
            if(numngb_inbox>0)
            {
                for(n=0; n<numngb_inbox; n++)
                {
                    j = Ngblist[n];
                    star_age = evaluate_stellar_age_Gyr(P[j].StellarAge);
                    m_st_kernel += P[j].Mass;
                    lm_ssp = evaluate_l_over_m_ssp(star_age) * calculate_relative_light_to_mass_ratio_from_imf(j);
                    l_st_kernel += lm_ssp*P[j].Mass;
                } /* for(n=0; n<numngb_inbox; n++) */
            } /* if(numngb_inbox>0) */
        } while(startnode >= 0);
        
        if(l_st_kernel>0) {
            m_gas_kernel_i = NORM_COEFF*SphP[i].Density*PPP[i].Hsml*PPP[i].Hsml*PPP[i].Hsml;
            rho_to_launch = SphP[i].Density*a3inv;
            m_st_kernel_i=m_st_kernel*(PPP[i].Hsml*PPP[i].Hsml*PPP[i].Hsml/(h*h*h));
            h_kernel_i=PPP[i].Hsml;
            
            /* calculate gas mass in same kernel properly to appropriately share the luminosity */
            m_gas_kernel=0; pos=P[i].Pos;
            startnode = All.MaxPart; dummy=0;
            do {
                numngb_inbox = ngb_treefind_variable(pos,h,-1,&startnode,0,&dummy,&dummy);
                if(numngb_inbox>0) for(n=0; n<numngb_inbox; n++) m_gas_kernel+=P[Ngblist[n]].Mass;
            } while(startnode >= 0);
            //printf("wind h %g numngb %d m_gas_kernel %g \n",h,numngb_inbox,m_gas_kernel);
            if(m_gas_kernel>0) {
                m_gas_kernel_i=m_gas_kernel;
                m_st_kernel_i=m_st_kernel;
                h_kernel_i=h;
            }
            rho_to_launch *= (1 + m_st_kernel_i/m_gas_kernel_i);
        }
        
#ifdef GALSF_FB_RPWIND_FROMCLUMPS
        /* ok, here we do a chained search of neighbors to find the local density peak,
         and use that to define the clump center (and escape velocity wrt that center) */
        
        vq=1; /* use as flag for whether to continue search iteration for a new density peak */
        k=i;  /* start at current particle */
        kmin=i; dmin1w=SphP[k].Density*1000.0;
        ngb_run_cntr=0; h=0;
        do {
            startnode = All.MaxPart; dummy=0;
            h=3.0*PPP[k].Hsml; pos=P[k].Pos;
            dmax1w=SphP[k].Density; vq=0;
            do {
                numngb_inbox=ngb_treefind_variable(pos,h,-1,&startnode,0,&dummy,&dummy);
                if(numngb_inbox>0) {
                    for(n=0; n<numngb_inbox; n++) {
                        j = Ngblist[n];
                        dmax2w=SphP[j].Density;
                        if((ngb_run_cntr==0)||(kmin==i)) {
                            if((dmax2w<dmin1w)&&(j != i)) {
                                dmin1w=dmax2w;
                                kmin=j;
                            }
                        } /* looks for nearby low-density particles to ID density gradient */
                        if(dmax2w>dmax1w)
                        {
                            dmax1w=dmax2w;
                            k=j;
                            vq=1; /* more dense particle found, continue chain */
                        }
                    } /* for(n=0; n<numngb_inbox; n++) */
                } /* if(numngb_inbox>0) */
            }while(startnode>=0);
            ngb_run_cntr++;
            h=0; for(j=0;j<3;j++) h+=(P[i].Pos[j]-P[k].Pos[j])*(P[i].Pos[j]-P[k].Pos[j]);
            h=sqrt(h);
        }while((vq==1)&&(ngb_run_cntr<20)&&(h<=20.0*PPP[i].Hsml));
        
        if(k==i)
        {
            m_gas_kernel = 8.*NORM_COEFF*SphP[i].Density*PPP[i].Hsml*PPP[i].Hsml*PPP[i].Hsml;
            h = 2.0*PPP[i].Hsml;
        }
        else
        {
            h=0;
            for(j=0;j<3;j++) h+=(P[i].Pos[j]-P[k].Pos[j])*(P[i].Pos[j]-P[k].Pos[j]);
            h=sqrt(h);
            if(h<2.0*PPP[i].Hsml) h=2.0*PPP[i].Hsml;
            if(h<2.0*PPP[k].Hsml) h=2.0*PPP[k].Hsml;
            
            m_gas_kernel=0;
            startnode = All.MaxPart; dummy=0;
            do {
                pos=P[k].Pos;
                numngb_inbox = ngb_treefind_variable(pos,h,-1,&startnode,0,&dummy,&dummy);
                if(numngb_inbox>0) for(n=0; n<numngb_inbox; n++) m_gas_kernel+=P[Ngblist[n]].Mass;
            } while(startnode >= 0);
        }
        
        float stcom_pos[3],cl_st_lum=0;
        m_st_kernel=0; //cl_st_lum_nonrad=0;
        startnode = All.MaxPart; dummy=0;
        do {
            pos=P[k].Pos;
            numngb_inbox = ngb_treefind_newstars(pos,h,-1,&startnode,0,&dummy,&dummy);
            stcom_pos[0]=0.0;stcom_pos[1]=0.0;stcom_pos[2]=0.0;
            if(numngb_inbox>0) { for(n=0; n<numngb_inbox; n++) {
                star_age = evaluate_stellar_age_Gyr(P[Ngblist[n]].StellarAge);
                vq = 10.1+500*pow(m_st_kernel_i*unitmass_in_e10solar,0.62)*
                pow(rho_to_launch*unitrho_in_e10solar_kpc3,-0.5);
                if(star_age <= vq) {
                    m_st_kernel+=P[Ngblist[n]].Mass;
                    lm_ssp = evaluate_l_over_m_ssp(star_age) * calculate_relative_light_to_mass_ratio_from_imf(Ngblist[n]);
                    cl_st_lum+=lm_ssp*P[Ngblist[n]].Mass;
                    
                    stcom_pos[0]+=P[Ngblist[n]].Pos[0]*lm_ssp*P[Ngblist[n]].Mass;
                    stcom_pos[1]+=P[Ngblist[n]].Pos[1]*lm_ssp*P[Ngblist[n]].Mass;
                    stcom_pos[2]+=P[Ngblist[n]].Pos[2]*lm_ssp*P[Ngblist[n]].Mass;
                }
            }
                stcom_pos[0]/=cl_st_lum;
                stcom_pos[1]/=cl_st_lum;
                stcom_pos[2]/=cl_st_lum;
                if(cl_st_lum>0&&m_gas_kernel>0) {
                    if(cl_st_lum/m_gas_kernel > l_st_kernel/m_gas_kernel_i) {
                        l_st_kernel=cl_st_lum;
                        m_gas_kernel_i=m_gas_kernel;
                        //l_st_kernel_nonrad=cl_st_lum_nonrad;
                    }
                }
            }
        } while(startnode >= 0);
#endif // GALSF_FB_RPWIND_FROMCLUMPS
        
        /* alright, now do calculations on the results to determine speed & loading of kicks */
        vq= WindInitialVelocityBoost*sqrt(All.G*(m_gas_kernel_i+m_st_kernel_i)/(h_kernel_i*All.cf_atime));
        v = WindInitialVelocityBoost*sqrt(All.G*(m_gas_kernel+m_st_kernel)/(h*All.cf_atime));
        if(vq>v) v=vq;
        /* compare the velocity from the central star cluster, using the observed cluster size-mass relation */
        if((m_st_kernel==0)&&(m_st_kernel_i>0)) m_st_kernel=m_st_kernel_i;
        vq=WindInitialVelocityBoost*(65.748/unitvel_in_km_s)*pow(m_st_kernel*unitmass_in_e10solar/0.0001,0.25);
        vq=vq*1.82;
        /* this corresponds to =G M_star/R_e for a 10^6 Msun cluster, scaling upwards from there;
         note that All.WindEnergyFraction will boost appropriately; should be at least sqrt(2), if want
         full velocities; in fact for Hernquist profile, the mass-weighted V_esc=1.82 times this */
        if (vq>v) v=vq;
        if(v<0.01/unitvel_in_km_s) v=0.01/unitvel_in_km_s;
        if(v>1000./unitvel_in_km_s) v=1000./unitvel_in_km_s;
        
        /* and now their mass-loading */
        /* winds driven by young stars, rather than proportional to SFR */
        p = 20.52*(dtime/unittime_in_gyr)*(l_st_kernel/m_gas_kernel_i)/(v/unitvel_in_km_s);
        
        myfree(Ngblist); // free neighbor lists used for wind kick searches
#endif // GALSF_FB_RPWIND_FROMSTARS
        
        
        // now calculate the 'boost factor' for the local properties //
        tau_IR = 1.91*SphP[i].Density*All.cf_a3inv*PPP[i].Hsml*All.cf_atime;
        if((m_gas_kernel>0)&&(h>0)) {
            vq=0.48*m_gas_kernel/(h*h*All.cf_atime*All.cf_atime);
            if(vq>tau_IR) tau_IR=vq;
        }
        tau_IR *= KAPPA_IR * All.UnitDensity_in_cgs*All.HubbleParam*All.UnitLength_in_cm;
        tau_IR *= (0.1+P[i].Metallicity[0]/All.SolarAbundances[0]);
        p *= (1.0 + All.WindMomentumLoading*tau_IR); //p += p_nonrad; // now nonrad is explicitly tracked in SNe/gasreturn
    } // closes(All.WindMomentumLoading>0)
    prob=p; if(prob>1) v*=p;
#endif // GALSF_FB_RPWIND_LOCAL
    
    
    
#if defined(GALSF_FB_RPWIND_DO_IN_SFCALC) && defined(GALSF_FB_RPWIND_LOCAL) // set values to return to previous step (for tabulation) //
    pvtau_return[0]=0.0;
#endif
#ifdef GALSF_FB_RPWIND_CONTINUOUS
    prob = 2;
#endif
    if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
    {
#ifdef GALSF_FB_RPWIND_CONTINUOUS
        v *= p;
#endif
        
#if defined(GALSF_FB_RPWIND_DO_IN_SFCALC) && defined(GALSF_FB_RPWIND_LOCAL) // set values to return to previous step (for tabulation) //
        pvtau_return[0]=1.0;
        pvtau_return[1]=p;
        pvtau_return[2]=v;
        pvtau_return[3]=tau_IR;
#endif
        
        // determine the wind acceleration orientation //
#ifdef GALSF_FB_RPWIND_FROMSTARS
        for(j=0;j<3;j++) dir[j]=-P[i].GradRho[j]; // default is along opacity gradient //
#endif
        
#if defined(GALSF_WINDS_POLAR) || !defined(GALSF_FB_RPWIND_FROMSTARS) || defined(GALSF_WINDS_RADIAL) // polar wind (defined by acel.cross.vel)
        dir[0] = P[i].GravAccel[1] * P[i].Vel[2] - P[i].GravAccel[2] * P[i].Vel[1];
        dir[1] = P[i].GravAccel[2] * P[i].Vel[0] - P[i].GravAccel[0] * P[i].Vel[2];
        dir[2] = P[i].GravAccel[0] * P[i].Vel[1] - P[i].GravAccel[1] * P[i].Vel[0];
#ifdef GALSF_WINDS_RADIAL  // winds always ejected outwards from galaxy; requires FOFRAD
	if( P[i].galID > 0 ) {	// check that particle is actually in a galaxy
	    theta = dir[0]*(P[i].Pos[0]-Groups[P[i].galID].COM[0])	// dot product -> cos(theta)
			+ dir[1]*(P[i].Pos[1]-Groups[P[i].galID].COM[1])
	    		+ dir[2]*(P[i].Pos[2]-Groups[P[i].galID].COM[2]);
            if( theta < 0 ) for(j=0;j<3;j++) dir[j] *= -1;
	}
        else if(get_random_number(P[i].ID + 5) < 0.5) {for(j=0;j<3;j++) dir[j]=-dir[j];}
	//printf("RADIAL: P%d %d %d %g %g %g %g %g %g %g\n",ThisTask,i,P[i].galID,dir[0],dir[1],dir[2],P[i].Pos[0]-Groups[P[i].galID].COM[0],P[i].Pos[1]-Groups[P[i].galID].COM[1],P[i].Pos[2]-Groups[P[i].galID].COM[2],theta);
#else
        if(get_random_number(P[i].ID + 5) < 0.5) {for(j=0;j<3;j++) dir[j]=-dir[j];}
#endif	// RADIAL
#endif
        
#ifdef GALSF_FB_RPWIND_FROMCLUMPS // in this case wind is directed from the local clump center //
        if(i != k)
        {
            for(j=0;j<3;j++) dir[j]=P[i].Pos[j]-P[k].Pos[j];
        } else {
            for(j=0;j<3;j++) dir[j]=P[kmin].Pos[j]-P[i].Pos[j];
        }
#endif
        
#ifdef GALSF_WINDS_ISOTROPIC // here the direction is random //
        theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
        phi = 2 * M_PI * get_random_number(P[i].ID + 4);
        dir[0] = sin(theta) * cos(phi); dir[1] = sin(theta) * sin(phi); dir[2] = cos(theta);
        if(get_random_number(P[i].ID + 5) < 0.5) {for(j=0;j<3;j++) dir[j]=-dir[j];}
#endif
        
        // now actually do the kick for the wind //
        for(j=0,norm=0;j<3;j++) norm+=dir[j]*dir[j];
        if(norm>0) {norm=sqrt(norm);} else {dir[0]=dir[1]=0; dir[2]=norm=1;}
        for(j = 0; j < 3; j++)
        {
            P[i].Vel[j] += v * All.cf_atime * dir[j]/norm;
            SphP[i].VelPred[j] += v * All.cf_atime * dir[j]/norm;
        }
#ifdef GALSF_SUBGRID_WINDS
            SphP[i].DelayTime = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
	    SphP[i].N_windLaunches += 1;

#ifdef FOFRAD
#ifdef GALSF_SUBGRID_DS86WIND
/* Heats wind in galaxies with v_circ<v_crit according to Dekel & Silk 86 formulation .
 * The idea is that in this regime SNe provide the energy to unbind, while rad pressure provides the
 * kinetic energy.  The remaining SN energy goes into heat, accounting for radiative cooling
 * during the SN snowplot phase. */
        double f_g=1, du=0, Tnew=0;
        double densfac = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv; /* converts density to physical g/cm^3 */
        double ufacinv = (All.UnitDensity_in_cgs/All.UnitPressure_in_cgs);  /* from cgs to system */

	f_g = (Groups[P[i].galID].Mbar-Groups[P[i].galID].Mstar)/Groups[P[i].galID].Mbar;
	if( f_g < 0.1 ) f_g = 0.1;
	du = 0.5*(GALSF_SUBGRID_DS86WIND*GALSF_SUBGRID_DS86WIND*pow(f_g,-0.26)-v_circ*v_circ)*All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s*All.cf_atime*All.cf_atime;
	if( du > 0 ) {
	    SphP[i].InternalEnergy += du*ufacinv;
            SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
            SphP[i].Pressure = get_pressure(i);
            Tnew = convert_u_to_temp(SphP[i].InternalEnergy/ufacinv, Particle_density_for_energy_i(i)*densfac,&SphP[i].Ne,i);  /* update electron abundance */
        //    if( ThisTask==0 ) printf("DS86WIND: z=%g Mb= %.2f fg= %.3f vc= %.1f vcrit= %.1f eta= %.2f T= %.3f %.3f\n",1./All.cf_atime-1,log10(Groups[P[i].galID].Mbar),f_g,v_circ,GALSF_SUBGRID_DS86WIND*pow(f_g,-0.13),wind_mass/sm,LogTemp(SphP[i].InternalEnergy/ufacinv-du,SphP[i].Ne),log10(Tnew));
	}
#endif	// DS86WIND
#endif	// FOFRAD

#ifdef GALSF_SUBGRID_HOTWIND
/* Heats wind with remaining SNe energy.  GALSF_SUBGRID_HOTWIND specifies fraction of supernova energy used to heat gas that is expelled in a wind (after kinetic energy has been accounted for).  0.0102778 is # of SNe (>8Mo) formed from 1 Mo of stars for Kroupa 2001 IMF. */
        double u_wind, u_SN, Tnew;
        double densfac = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv; /* converts density to physical g/cm^3 */
        double ufacinv = (All.UnitDensity_in_cgs/All.UnitPressure_in_cgs);  /* from cgs to system */

        u_wind = 0.5*v*v*All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s*All.cf_atime*All.cf_atime;
        u_SN = 1.e51 * (0.0102778/SOLAR_MASS) * sm/wind_mass;	// convert SN energy per unit gas mass into erg/g, distributed among eta=wind_mass/sm particles
        if( u_SN > u_wind ) {
            SphP[i].InternalEnergy += GALSF_SUBGRID_HOTWIND*(u_SN-u_wind)*ufacinv*(0.5+get_random_number(P[i].ID + 51));
            SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
            SphP[i].Pressure = get_pressure(i);
            Tnew = convert_u_to_temp(SphP[i].InternalEnergy/ufacinv, Particle_density_for_energy_i(i)*densfac,&SphP[i].Ne,i);  /* update electron abundance */
            //if( ThisTask==0 ) printf("HOTWIND: a=%g p= %g vw= %g eta= %g u= %g u_wind,SN= %g %g %g u/uold= %g Tnew= %g %g\n",All.Time,prob,v,wind_mass/sm,SphP[i].InternalEnergy,u_wind,u_SN,u_SN/u_wind,SphP[i].InternalEnergy/(SphP[i].InternalEnergy-(u_SN+u_wind)*ufacinv),Tnew,LogTemp(SphP[i].InternalEnergy/ufacinv,SphP[i].Ne));
        }
#endif	// HOTWIND
#endif	// GALSF_SUBGRID_WINDS
    } /* if(get_random_number(P[i].ID + 2) < prob) */
}
#endif // defined(GALSF_FB_RPWIND_DO_IN_SFCALC) || defined(GALSF_SUBGRID_WINDS) || defined(GALSF_SUBGRID_VARIABLEVELOCITY) //




#ifdef GALSF_EFFECTIVE_EQS
/* Routine to initialize quantities needed for the Spingel & Hernquist effective equation of state */
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, thresholdStarburst;
#ifdef COOL_METAL_LINES_BY_SPECIES
  int k; double Z[NUM_METAL_SPECIES]; for(k=0;k<NUM_METAL_SPECIES;k++) Z[k]=0;
#endif

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;
      egyhot = All.EgySpecSN / A0;
      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
      dens = 1.0e6 * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
      tcool = GetCoolingTime(egyhot, dens, &ne, -1);
      coolrate = egyhot / tcool / dens;
      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh = x / pow(1 - x, 2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      if(ThisTask == 0)
	{
	  printf("\nA0= %g  \n", A0);
	  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;
      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
      tcool = GetCoolingTime(egyhot, dens, &ne, -1);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;
	  peff = GAMMA_MINUS1 * dens * egyeff;
	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;
	  neff = -log(peff) * fac;
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
      tcool = GetCoolingTime(egyhot, dens, &ne, -1);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;
	  peff = GAMMA_MINUS1 * dens * egyeff;
	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;
#ifdef MODIFIEDBONDI
      All.BlackHoleRefDensity = thresholdStarburst;
      All.BlackHoleRefSoundspeed = sqrt(GAMMA * GAMMA_MINUS1 * egyeff);
#endif

      if(ThisTask == 0)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

#ifdef WINDS
#ifndef GALSF_FB_RPWIND_LOCAL
      if(All.WindEfficiency > 0)
	if(ThisTask == 0)
	  printf("Windspeed: %g\n",
		 sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency));
#endif
#endif
    }
}
#endif // GALSF_EFFECTIVE_EQS //

#ifdef FOF_SCALEDEPENDENT
double get_next_FoF_time(double haloDynamicalTime){
  double growth_timescale, redshift, da;
  redshift = 1.0 / All.Time - 1.0;
  growth_timescale = haloDynamicalTime * 2.5 / pow(1.0 + redshift,1.5);

  // large # to convert H(z) from km/s/Mpc --> 1/Gyr
  da = All.HubbleParam * 100. * All.Time * growth_timescale * 0.001022914943273906 *
    sqrt(All.OmegaLambda + All.Omega0 * pow(1. + redshift,3));

  return (All.Time+da);
}
#endif

//#ifdef BH_BL
double FreeFallTime(int i)
{
	double density = SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam;
	return sqrt(3*PI/32/density/GRAVITY)/SEC_PER_YEAR;
}
//#endif

#ifdef BLACK_HOLES
void initBH(int i)
{
	double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
	BPP(i).BH_Mass = P[i].BH_mass * Msun_to_Internal;
	BPP(i).BH_accreted = 0;
#ifdef BH_ALPHADISK_ACCRETION                                     // DAA: this should be outside of BH_SEED_STAR_MASS_FRACTION ...
    BPP(i).BH_Mass_AlphaDisk = 0;
#endif
    BPP(i).BH_Mdot = 0;
#ifdef BH_COUNTPROGS
    BPP(i).BH_CountProgs = 1;
#endif

#ifdef BH_BUBBLES
    P[i].BH_Mass_bubbles = P[i].BH_mass * Msun_to_Internal;
    P[i].BH_Mass_ini = P[i].BH_mass * Msun_to_Internal;
#ifdef UNIFIED_FEEDBACK
    P[i].BH_Mass_radio = P[i].BH_mass * Msun_to_Internal;
#endif
#endif
}
#endif

#endif // GALSF


