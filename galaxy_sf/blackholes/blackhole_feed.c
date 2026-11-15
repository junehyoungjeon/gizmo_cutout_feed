/*! \file blackhole_feed.c
 *  \brief This is where particles are marked for gas accretion.
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   on 1/9/15 by Paul Torrey (ptorrey@mit.edu) for clarity by parsing the existing code into
 *   smaller files and routines. Some communication and black hole structures were modified
 *   to reduce memory usage. Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
#include "blackhole_local.h"

#if (!defined(BH_WIND_BETA)) && defined(BH_BAL_KICK)
#define BH_WIND_BETA 1.0
#endif

void blackhole_feed_loop(void)
{
    int i, j, k, ndone_flag, ndone;
    int ngrp, recvTask, place, nexport, nimport, dummy;
    MPI_Status status;
    
    /* allocate buffers to arrange communication */
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct blackholedata_in) +
                                                             sizeof(struct blackholedata_out) +
                                                             sizemax(sizeof(struct blackholedata_in),
                                                                     sizeof(struct blackholedata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    /* Let's determine which particles may be swallowed by whom, and the weights for feedback */
    i = FirstActiveParticle;
    do
    {
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        
        /* do local particles and prepare export list */
        for(nexport = 0; i >= 0; i = NextActiveParticle[i])                      // DAA: can this be replaced by a loop over N_active_loc_BHs
            if(P[i].Type == 3)
                if(blackhole_feed_evaluate(i, 0, &nexport, Send_count) < 0)
                    break;
        
        MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        BlackholeDataGet = (struct blackholedata_in *) mymalloc("BlackholeDataGet", nimport * sizeof(struct blackholedata_in));
        BlackholeDataIn = (struct blackholedata_in *) mymalloc("BlackholeDataIn", nexport * sizeof(struct blackholedata_in));
        
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
            
            for(k = 0; k < 3; k++)
            {
                BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
                BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                BlackholeDataIn[j].Jgas_in_Kernel[k] = P[place].GradRho[k];
#endif
            }
            
#ifdef GW_TIME
#ifdef OTF_BULGE
			BlackholeDataIn[j].star_rms_vel = BlackholeTempInfo[P[place].IndexMapToTempStruc].star_rms_vel;
			BlackholeDataIn[j].rho_star = BlackholeTempInfo[P[place].IndexMapToTempStruc].rho_star;
#endif
#endif
            
#if defined(BH_GRAVCAPTURE_GAS)
            BlackholeDataIn[j].mass_to_swallow_edd = BlackholeTempInfo[P[place].IndexMapToTempStruc].mass_to_swallow_edd;
#endif
            BlackholeDataIn[j].Hsml = PPP[place].Hsml;
            BlackholeDataIn[j].Mass = P[place].Mass;
            BlackholeDataIn[j].BH_Mass = BPP(place).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
            BlackholeDataIn[j].BH_Mass_AlphaDisk = BPP(place).BH_Mass_AlphaDisk;
#endif
#if defined(BH_PHOTONMOMENTUM) 	|| defined(BH_BAL_WINDS)
            BlackholeDataIn[j].BH_disk_hr = P[place].BH_disk_hr;
#endif
            BlackholeDataIn[j].Density = BPP(place).DensAroundStar;
            BlackholeDataIn[j].Mdot = BPP(place).BH_Mdot;
#ifndef WAKEUP
            BlackholeDataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
            BlackholeDataIn[j].Dt = P[place].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
            BlackholeDataIn[j].ID = P[place].ID;
            #ifdef BH_FB_CONE
            {
                double Jv[3], nv[3]; int kk;
                for(kk=0;kk<3;kk++) Jv[kk] = BlackholeTempInfo[P[place].IndexMapToTempStruc].Jgas_in_Kernel[kk];
                bh_fb_unit_axis(Jv, nv);
                for(kk=0;kk<3;kk++) BlackholeDataIn[j].BH_fb_axis[kk] = nv[kk];
            }
#endif
            memcpy(BlackholeDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        
        /* exchange particle data */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_E,
                                 &BlackholeDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                                 recvTask, TAG_BH_E, MPI_COMM_WORLD, &status);
                }
            }
        }
        myfree(BlackholeDataIn);
        BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult",nimport * sizeof(struct blackholedata_out));
        BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));
        
        /* now do the particles that were sent to us */
        for(j = 0; j < nimport; j++)
            blackhole_feed_evaluate(j, 1, &dummy, &dummy);
        
        
        if(i < 0)
            ndone_flag = 1;
        else
            ndone_flag = 0;
        
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_F,
                                 &BlackholeDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct blackholedata_out),
                                 MPI_BYTE, recvTask, TAG_BH_F, MPI_COMM_WORLD, &status);
                }
            }
        } // for(ngrp = 1; ngrp < (1 << PTask); ngrp++) //
        
        /* add the result to the particles */
        for(j = 0; j < nexport; j++)
        {
            place = DataIndexTable[j].Index;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
            BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_angle_weighted_kernel_sum += BlackholeDataOut[j].BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_FB_CONE
            {int kk; for(kk=0;kk<3;kk++) BlackholeTempInfo[P[place].IndexMapToTempStruc].BH_fb_sum[kk] += BlackholeDataOut[j].BH_fb_sum[kk];}
#endif
        }
        
        myfree(BlackholeDataOut);
        myfree(BlackholeDataResult);
        myfree(BlackholeDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
}







/* do loop over neighbors to get quantities for accretion */
int blackhole_feed_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
    int startnode, numngb, j, k, n, listindex = 0;
    MyIDType id;
    MyFloat *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass;
    double h_i2, r2, r, u, hinv, hinv3, wk, dwk, vrel, vesc, dpos[3];
//#if defined() || defined()
    double meddington;
    double eta;
    double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
//#endif
    
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
    double medd_max_accretable, mass_to_swallow_edd, eddington_factor;
#endif
    
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    double norm, theta, BH_disk_hr, *Jgas_in_Kernel;
    double BH_angle_weighted_kernel_sum=0;
#endif

#if defined(BH_BAL_KICK) && !defined(BH_GRAVCAPTURE_GAS)
    double f_accreted=0; 
#endif
    
#ifdef BH_THERMALFEEDBACK
    double energy;
#endif
#ifdef BH_FB_CONE
    double fb_axis[3]={0,0,0}, fb_sum[3]={0,0,0}, fb_c, fb_f; int fb_axis_ok=0;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat bh_mass_alphadisk;
#endif
#if defined(BH_SWALLOWGAS)
    double w=0,p=0,mass_markedswallow=0,bh_mass_withdisk=0;
#endif
    
#ifdef GW_TIME
	MyFloat sigma_inf=0, rho_inf=0;
	double ecc=0, eps=0, hang2=0, mu=0, mred=0, mtot=0;
	double eps_cgs=0;
	double dvel[3], ang[3];
	//double mc_vel[3]; mc_pos[3];
#endif
    /* these are the BH properties */
    if(mode == 0)
    {
        pos = P[target].Pos;
        rho = P[target].DensAroundStar;       // DAA: DensAroundStar is not defined in BHP->BPP...
        mdot = BPP(target).BH_Mdot;
#ifndef WAKEUP
        dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
        dt = P[target].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
        h_i = PPP[target].Hsml;
        mass = P[target].Mass;
        bh_mass = BPP(target).BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BPP(target).BH_Mass_AlphaDisk;
#endif
        velocity = P[target].Vel;
        id = P[target].ID;
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = P[target].GradRho;
        BH_disk_hr = P[target].BH_disk_hr;
#endif
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
        mass_to_swallow_edd = BlackholeTempInfo[P[target].IndexMapToTempStruc].mass_to_swallow_edd;
#endif
#ifdef GW_TIME
#ifdef OTF_BULGE
		sigma_inf = BlackholeTempInfo[P[target].IndexMapToTempStruc].star_rms_vel;
		rho_inf = BlackholeTempInfo[P[target].IndexMapToTempStruc].rho_star;
#endif
#endif
    }
    else
    {
        pos = BlackholeDataGet[target].Pos;
        rho = BlackholeDataGet[target].Density;
        mdot = BlackholeDataGet[target].Mdot;
        dt = BlackholeDataGet[target].Dt;
        h_i = BlackholeDataGet[target].Hsml;
        mass = BlackholeDataGet[target].Mass;
        bh_mass = BlackholeDataGet[target].BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
        bh_mass_alphadisk = BlackholeDataGet[target].BH_Mass_AlphaDisk;
#endif
        velocity = BlackholeDataGet[target].Vel;
        id = BlackholeDataGet[target].ID;
#if defined(BH_PHOTONMOMENTUM)  || defined(BH_BAL_WINDS)
        Jgas_in_Kernel = BlackholeDataGet[target].Jgas_in_Kernel;
        BH_disk_hr = BlackholeDataGet[target].BH_disk_hr;
#endif
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
        mass_to_swallow_edd = BlackholeDataGet[target].mass_to_swallow_edd;
#endif
#ifdef GW_TIME
#ifdef OTF_BULGE
		sigma_inf = BlackholeDataGet[target].star_rms_vel;
		rho_inf = BlackholeDataGet[target].rho_star;
#endif
#endif
    }
    #ifdef BH_FB_CONE
    if(mode == 0) {
        double Jv[3]; for(k=0;k<3;k++) Jv[k] = BlackholeTempInfo[P[target].IndexMapToTempStruc].Jgas_in_Kernel[k];
        bh_fb_unit_axis(Jv, fb_axis);
    } else {
        for(k=0;k<3;k++) fb_axis[k] = BlackholeDataGet[target].BH_fb_axis[k];
    }
    fb_axis_ok = (fb_axis[0]*fb_axis[0] + fb_axis[1]*fb_axis[1] + fb_axis[2]*fb_axis[2] > 0.5);
#endif
    if((mass<0)||(h_i<=0)) return -1;
    
    /* initialize variables before SPH loop is started */
    h_i2 = h_i * h_i;
    hinv = 1 / h_i;
    hinv3 = hinv * hinv * hinv;
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
    meddington = (4 * M_PI * GRAVITY * Clight * PROTONMASS / (All.BlackHoleRadiativeEfficiency * Clight * Clight * THOMPSON)) * (bh_mass/All.HubbleParam) * All.UnitTime_in_s;
    medd_max_accretable = All.BlackHoleEddingtonFactor * meddington * dt;
    eddington_factor = mass_to_swallow_edd / medd_max_accretable;   /* if <1 no problem, if >1, need to not set some swallowIDs */
#endif
    
#if defined(BH_SWALLOWGAS)
    bh_mass_withdisk = bh_mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += bh_mass_alphadisk;
#endif
#endif

#if defined(BH_BAL_KICK) && !defined(BH_GRAVCAPTURE_GAS)
  
#ifdef BH_FDBK_Z_THRESH
  f_accreted = 1.0;
  if (All.Time>1.0/(1.0+BH_FDBK_Z_THRESH))
#endif
  {
#ifdef BH_BAL_KICK_MOMENTUM_FLUX
    /* DAA: increase the effective mass-loading of BAL winds to reach the desired momentum flux given the outflow velocity "All.BAL_v_outflow" chosen
       --> appropriate for cosmological simulations where particles are effectively kicked from ~kpc scales
           (i.e. we need lower velocity and higher mass outflow rates compared to accretion disk scales) - */
    f_accreted = 1. / ( 1. + BH_BAL_KICK_MOMENTUM_FLUX * All.BlackHoleRadiativeEfficiency * (Clight / All.UnitVelocity_in_cm_per_s) / All.BAL_v_outflow );
#else
#ifdef ROMULUS_WIND
	double eta_w;
	double eps_w;
	meddington = (4*M_PI*GRAVITY*Clight*PROTONMASS/(All.BlackHoleRadiativeEfficiency*Clight*Clight*THOMPSON))*bh_mass*All.UnitTime_in_s/All.HubbleParam;	
	eta = mdot / meddington;
	eps_w = eta * All.BH_Awind * All.BlackHoleWindEfficiency/(1 + All.BH_Awind * eta);
    double v_kick = All.BAL_v_outflow*1e5;
#ifdef BH_MASS_WIND
	double BH_mass_wind = BH_MASS_WIND * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g;
	v_kick *= pow(bh_mass_withdisk / BH_mass_wind, BH_WIND_BETA);
#endif
	eta_w = 2*eps_w * pow(Clight / v_kick, 2);
    f_accreted = 1./(1. + eta_w);
#endif
#endif
  }
#endif //#ifdef BH_BAL_KICK
    
    /* Now start the actual SPH computation for this BH particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = BlackholeDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    BH_angle_weighted_kernel_sum = 0;
#endif
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);
            if(numngb < 0) return -1;
            
            for(n = 0; n < numngb; n++)
            {
                j = Ngblist[n];
                if(P[j].Mass > 0)
                {
                    for(k=0;k<3;k++) dpos[k] = pos[k] - P[j].Pos[k]; 
#ifdef PERIODIC
                    NEAREST_XYZ(dpos[0],dpos[1],dpos[2],1);
#endif
                    r2=0; for(k=0;k<3;k++) r2+=dpos[k]*dpos[k];
                    
                    if(r2 < h_i2)
                    {
                        
                        vrel = 0;
                        for(k=0;k<3;k++) vrel += (P[j].Vel[k] - velocity[k])*(P[j].Vel[k] - velocity[k]);
                        vrel = sqrt(vrel) / All.cf_atime;       /* do this once and use below */
                        vesc = sqrt(2.0*All.G*(mass+P[j].Mass)/(sqrt(r2)*All.cf_atime));// + pow(10.e5/All.UnitVelocity_in_cm_per_s,2));
                        r = sqrt(r2);
#ifdef BH_MERGER
                        /* check_for_bh_merger.  Easy.  No Edd limit, just a pos and vel criteria. */
                        if((id != P[j].ID) && (P[j].Mass > 0) && (P[j].Type == 3) && (r < 2*All.ForceSoftening[3]) && (bh_mass>0) && (P[j].BH_mass>0) )	/* we may have a black hole merger */
                        {
                                if(vrel > BH_CSND_FRAC_BH_MERGE * vesc)
                                {
#ifdef BH_OUTPUT_MOREINFO           // DAA: BH merger info will be saved in a separate output file 
                                    printf("ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n",
                                            ThisTask, All.Time, id, P[j].ID, vrel, vesc);
#else
                                    fprintf(FdBlackHolesDetails,
                                            "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n",
                                            ThisTask, All.Time, id, P[j].ID, vrel, vesc);
#endif
                                }
                                else
                                {
#ifdef GW_TIME
                                	mtot = mass + P[j].Mass;
                                	mred = mass*P[j].Mass/(mass+P[j].Mass);
                                	eps = 0.5*vrel * vrel * mred - All.G*(mass*P[j].Mass)/(r*All.cf_atime);
                                	eps_cgs = eps*pow(All.UnitVelocity_in_cm_per_s,2)*All.UnitMass_in_g;
                                	eps /= mred; 
                                	mu = All.G*mtot;
                                	for(k=0;k<3;k++) dvel[k] = velocity[k]-P[j].Vel[k];
                                	ang[0] = dvel[1]*dpos[2] - dvel[2]*dpos[1];
                                	ang[1] = dvel[2]*dpos[0] - dvel[0]*dpos[2];
                                	ang[2] = dvel[0]*dpos[1] - dvel[1]*dpos[0];
                                	for(k=0;k<3;k++) hang2 += ang[k]*ang[k];
                                	ecc = 1 + 2*eps*hang2/(mu*mu);
                                	if(ecc < 0 && ecc >= 1)
                                	{
                                		fprintf(FdBlackHolesDetails,"Error in eccentricity for P[j.]ID=%llu to be swallowed by id=%llu \n",
                                           (unsigned long long) P[j].ID, (unsigned long long) id);
                                        ecc = 0;
                                	}
                                	else 
                                	{
                                		ecc = sqrt(ecc);
#ifdef HD_BINARY
#ifndef M_FS
#define M_FS 1
#endif
									if(-eps*mred*(1-ecc)/(1+ecc) > 0.5*pow(HD_BINARY/All.UnitVelocity_in_cm_per_s,2) * M_FS *Msun_to_Internal)
#endif
#endif
								{
                                    printf("MARKING_BH_MERGER: P[j.]ID=%llu to be swallowed by id=%llu \n",
                                           (unsigned long long) P[j].ID, (unsigned long long) id);

                                    //if(P[j].SwallowID < id && P[j].ID < id) // makes it so only one swallows the other
                                    // DAA: makes it so that the most massive BH swallows the other - simplifies analysis
                                    if((P[j].SwallowID == 0) && (BPP(j).BH_Mass < bh_mass)) 
                                    {
                                        P[j].SwallowID = id;
#ifdef GW_TIME	
#ifdef OTF_BULGE
									if(rho_inf <= 0 || sigma_inf <= 0)
#endif
									{
										rho_inf = rho_inf_est(mtot);
										sigma_inf = sigma_inf_est(mtot);
									}
#ifdef OTF_BULGE
#ifdef GW_DENS_BOOST
									else
									{
										rho_inf *= pow(rho_inf / (GW_DENS_BOOST*PROTONMASS/HYDROGEN_MASSFRAC), GW_BETA);
										if(rho_inf > rho_inf_est(mtot) * GW_RHO_MAX)
											rho_inf = rho_inf_est(mtot) * GW_RHO_MAX;
									}
#endif
#endif
									if (isnan(ecc)) ecc = 0;
                                        P[j].ecc = ecc;
                                        P[j].eps = eps_cgs;
                                        P[j].GW_time = coalescence_time(bh_mass, BPP(j).BH_Mass, sigma_inf, rho_inf, ecc);
#endif
                                    }
                                }
#ifdef GW_TIME
                                }
#endif
                                }
                        } // if(P[j].Type == 3) //
#endif
                        
                        
                        /* This is a similar loop to what we already did in blackhole_environment, but here we stochastically
                         reduce GRAVCAPT events in order to (statistically) obey the eddington limit */
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
                        if(P[j].Type != 3)
                        {
                            double cs=0; if(P[j].Type==0) {cs=Particle_effective_soundspeed_i(j);}
                            if(vrel < vesc){ /* bound */
                                if( All.ForceSoftening[3]*(1.0-(vrel*vrel+cs*cs)/(vesc*vesc))/r > 1.0 )
                                { /* apocenter within 2.8*epsilon (softening length) */
                                    
#ifdef BH_GRAVCAPTURE_NONGAS
                                    /* simply swallow non-gas particle if BH_GRAVCAPTURE_NONGAS enabled */
                                    if((P[j].Type != 0) && (P[j].SwallowID < id)) P[j].SwallowID = id;
#endif
                                    
#ifdef BH_GRAVCAPTURE_GAS
                                    /* now deal with gas */
                                    if (P[j].Type == 0){
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                        /* if Eddington-limited and NO alpha-disk, do this stochastically */
                                        p = 1/eddington_factor;
#if defined(BH_BAL_WINDS) || defined(BH_BAL_KICK)
										if(f_accreted>0) p /= f_accreted;
                                        else p /= All.BAL_f_accretion; // we need to accrete more, then remove the mass in winds
#endif
                                        w = get_random_number(P[j].ID);
                                        if(w < p) {
                                            printf("MARKING_BH_FOOD: P[j.]ID=%llu to be swallowed by id=%llu \n",
                                                   (unsigned long long) P[j].ID, (unsigned long long) id);
                                            if(P[j].SwallowID < id) P[j].SwallowID = id;
                                        } else { /* w < p */
                                            printf("MARKING_BH_FOOD (will be rejected): P[j.]ID=%llu to be swallowed by id=%llu \n",
                                                   (unsigned long long) P[j].ID, (unsigned long long) id);
                                            //if(P[j].SwallowID < id)  P[j].SwallowID = id; // rejected
                                        }
#else //if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                        /* in other cases, just swallow the particle */
                                        if(P[j].SwallowID < id) P[j].SwallowID = id;            
#endif //else defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                    } //if (P[j].Type == 0)
#endif //ifdef BH_GRAVCAPTURE_GAS
                                    
                                } // if( All.ForceSoftening[5]*(1.0-vrel*vrel/(vesc*vesc))/sqrt(r2) > 1.0 )
                            } // if(vrel < vesc)
                        } //if(P[j].Type != 3)
#endif // if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
                        
                        
                        
                        
                        /* now is the more standard accretion only of gas, according to the mdot calculated before */
                        if(P[j].Type == 0)
                        {
                            /* here we have a gas particle */
                            u = r * hinv;
                            kernel_main(u,hinv3,hinv*hinv3,&wk,&dwk,-1);
                            
// DAA: this below is only meaningful if !defined(BH_GRAVCAPTURE_GAS)...
//#ifdef BH_SWALLOWGAS
#if defined(BH_SWALLOWGAS) && !defined(BH_GRAVCAPTURE_GAS)
                            /* compute accretion probability */
                            if((bh_mass_withdisk - (mass + mass_markedswallow))>0)
                                p = (bh_mass_withdisk - (mass + mass_markedswallow)) * wk / rho;
                            else
                                p = 0;
                            
/* DAA: for stochastic winds (BH_BAL_KICK) we remove a fraction of mass from gas particles prior to kicking
 * --> need to increase the probability here to balance black hole growth   
 */
#ifdef BH_BAL_KICK
                            if(f_accreted>0) 
                            {
                                p /= f_accreted;

                                /* DAA: compute outflow probability when "bh_mass_withdisk < mass"
                                    - we don't need to enforce mass conservation in this case 
                                    - relevant only in low-res sims where the BH seed mass is much lower than the gas particle mass 
                                    - TODO: replace mdot below by mdot_alphadisk if BH_ALPHADISK_ACCRETION? */
/*#ifdef BH_MASS_RES
	double BH_mass_res = BH_MASS_RES * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g;
	if(bh_mass_withdisk > BH_mass_res)
#endif */
                                {
                                if((bh_mass_withdisk - mass) < 0)
                                    p = ( (1-f_accreted)/f_accreted ) * mdot * dt * wk / rho;
                                }
                            }
#endif

                            
                            w = get_random_number(P[j].ID);
                            if(w < p)
                            {
                                printf("MARKING_BH_FOOD: j %d w %g p %g TO_BE_SWALLOWED \n",j,w,p);
                                if(P[j].SwallowID < id) 
                                {
                                   P[j].SwallowID = id;
#ifdef BH_BAL_KICK
                                   mass_markedswallow += P[j].Mass*f_accreted;
#else
                                   mass_markedswallow += P[j].Mass;
#endif
                                }
                            } // if(w < p)
#endif // BH_SWALLOWGAS

                            
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
                            /* calculate the angle-weighting for the photon momentum */
                            if((mdot>0)&&(dt>0)&&(r>0)&&(P[j].SwallowID==0))
                            {
                                /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                                norm=0; for(k=0;k<3;k++) norm+=(dpos[k]/r)*Jgas_in_Kernel[k];
                                norm=fabs(norm); theta=acos(norm);
                                BH_angle_weighted_kernel_sum += bh_angleweight_localcoupling(j,BH_disk_hr,theta);
                            }
#endif
                            
#ifdef BH_THERMALFEEDBACK
#ifdef BH_FB_CONE
                            if(P[j].Mass > 0 && r > 0)
                            {
                                fb_f = 1.0;
                                if(fb_axis_ok)
                                {
                                    //fb_c = 0; for(k=0;k<3;k++) fb_c += (dpos[k]/r) * fb_axis[k];
                                    //fb_f = bh_cone_weight(fabs(fb_c));
                                    fb_c = 0; for(k=0;k<3;k++) fb_c -= (dpos[k]/r) * fb_axis[k];   /* minus: dpos points gas -> BH */
                                    fb_f = bh_cone_weight(fb_c);                                    /* no fabs */
                                }
                                fb_sum[1] += P[j].Mass * wk;
                                if(fb_f > 0) {fb_sum[0] += P[j].Mass * wk * fb_f; fb_sum[2] += 1.0;}
                            }
#else

#ifdef BH_FDBK_Z_THRESH
  if (All.Time>1.0/(1.0+BH_FDBK_Z_THRESH))
#endif
  {
						if (All.BlackHoleFeedbackFactor>0)
						{
							meddington = (4*M_PI*GRAVITY*Clight*PROTONMASS/(All.BlackHoleRadiativeEfficiency*Clight*Clight*THOMPSON))*bh_mass*All.UnitTime_in_s/All.HubbleParam;
#ifdef UNIFIED_FEEDBACK
                            if(mdot > All.RadioThreshold * meddington)
#endif
                            {
                            	eta = mdot / meddington;
                                energy = All.BlackHoleFeedbackFactor*All.BlackHoleRadiativeEfficiency * mdot*dt * pow(Clight/All.UnitVelocity_in_cm_per_s,2) * (All.BH_Arad * eta)/(1. + All.BH_Arad * eta);
                                if(rho > 0)
                                    SphP[j].Injected_BH_Energy += (wk/rho) * energy * P[j].Mass;
                            }
                        }
  }
#endif
#endif
                            
                        } // if(P[j].Type == 0)
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    } // if(r2 < h_i2)
                } // if(P[j].Mass > 0)
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = BlackholeDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // mode==1
    } // while(startnode >= 0) (outer of the double-loop)
    
    
    
    /* Now collect the result at the right place */
    if(mode == 0)
    {
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_angle_weighted_kernel_sum += BH_angle_weighted_kernel_sum;  /* need to correct target index */
#endif
#ifdef BH_FB_CONE
        //{int kk; for(kk=0;kk<3;kk++) BlackholeDataResult[target].BH_fb_sum[kk] = fb_sum[kk];}
        {int kk; for(kk=0;kk<3;kk++) BlackholeTempInfo[P[target].IndexMapToTempStruc].BH_fb_sum[kk] += fb_sum[kk];}
#endif
    }
    else
    {
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
        BlackholeDataResult[target].BH_angle_weighted_kernel_sum = BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_FB_CONE
        {int kk; for(kk=0;kk<3;kk++) BlackholeDataResult[target].BH_fb_sum[kk] = fb_sum[kk];}
#endif
    }
    return 0;
} /* closes bh_evaluate routine */
#ifdef GW_TIME
double coalescence_time(MyFloat M1, MyFloat M2, MyFloat sigma_inf, MyFloat rho_inf, double ecc)
{
	double Internal_to_g  = All.UnitMass_in_g / All.HubbleParam;
	double m1m2, m, sigma, rho, Fe, agw, e2 = ecc*ecc;
	Fe = pow(1-e2, -3.5) * (1 + (73./24.)*e2 + (37./96.)*e2*e2);
	m1m2 = M1 * M2 * pow(Internal_to_g, 2);
	m = (M1+M2) * Internal_to_g;
	sigma = sigma_inf;// * All.UnitVelocity_in_cm_per_s;
	rho = rho_inf;// * All.UnitDensity_in_cgs; 
	// unit convertion is done in normalize_temp_info_struct
	agw = pow( 64*GRAVITY*GRAVITY*sigma*m1m2*m*Fe/(5*pow(Clight,5)*GW_TIME*rho), 0.2);
	return sigma/(GRAVITY*GW_TIME*rho*agw) / SEC_PER_YEAR;
}

double sigma_inf_est(double mtot)
{
	double M = mtot * All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS;
	double sig = pow(M/1e9/0.309, 1/4.38) * 2e7;
	return sig;
}

double rho_inf_est(double mtot)
{
	double M = mtot * All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS;
	double r0, reff, rinf, Ms;
	Ms = pow(M/1e9/0.49, 1/1.16) * 1e11;
	reff = DMAX(2.95*pow(Ms/1e6, 0.596), 34.8*pow(Ms/1e6, 0.399));
	r0 = reff/0.75 * (pow(2, 1/(3-DP_GAMMA)) - 1);
	rinf = r0/(pow(Ms/(2*M), 1/(3-DP_GAMMA)) - 1);
	double rho = (3-DP_GAMMA) * Ms / (4*M_PI) * r0/(pow(rinf, DP_GAMMA)*pow(rinf+r0, 4-DP_GAMMA));
	rho *= SOLAR_MASS / pow(CM_PER_MPC/1e6, 3);
	return rho;
}
#endif
