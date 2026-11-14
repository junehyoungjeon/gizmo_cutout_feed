#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../../allvars.h"
#include "../../proto.h"

/*! \file blackhole_util.c
 *  \brief util routines for memory (de)allocation and array setting for black holes
 */
/*
 * This file was largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for black hole accretion and feedback have been
 *   replaced, and the algorithm for their coupling is new to GIZMO.  This file was modified
 *   by Paul Torrey (ptorrey@mit.edu) for clairity.  It was rearranged and parsed into
 *   smaller files and routines. The main functional difference is that BlackholeTempInfo
 *   is now allocated only for N_active_loc_BHs, rather than NumPart (as was done before).  Some
 *   extra index gymnastics are required to follow this change through in the MPI comm routines.
 *   Cleanup, de-bugging, and consolidation of routines by Xiangcheng Ma
 *   (xchma@caltech.edu) followed on 05/15/15; re-integrated by PFH.
 */

/* function for allocating temp BH data struc needed for feedback routines*/
void blackhole_start(void)
{
    int i, Nbh;
    
    /* count the num BHs on this task */
    N_active_loc_BHs=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==3)
        {
            P[i].IndexMapToTempStruc = N_active_loc_BHs;         /* allows access via BlackholeTempInfo[P[i].IndexMapToTempStruc] */
            N_active_loc_BHs++;                     /* N_active_loc_BHs now set for BH routines */
        }
    }
    
    /* allocate the blackhole temp struct -- defined in blackhole.h */
    if(N_active_loc_BHs>0)
    {
        BlackholeTempInfo = (struct blackhole_temp_particle_data *) mymalloc("BlackholeTempInfo", N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data));
    } else {
        BlackholeTempInfo = (struct blackhole_temp_particle_data *) mymalloc("BlackholeTempInfo", 1 * sizeof(struct blackhole_temp_particle_data));
    }
    
    memset( &BlackholeTempInfo[0], 0, N_active_loc_BHs * sizeof(struct blackhole_temp_particle_data) );
    
    Nbh=0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type ==3)
        {
            BlackholeTempInfo[Nbh].index = i;               /* only meaningful field set here */
            Nbh++;
        }
    }
    
    /* all future loops can now take the following form:
     for(i=0; i<N_active_loc_BHs; i++)
     {
     i_old = BlackholeTempInfo[i].index;
     ...
     }
     */
    
}


/* function for freeing temp BH data struc needed for feedback routines*/
void blackhole_end(void)
{
#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    {
        int bin, myrank;
        double mdot, mdot_in_msun_per_year;
        double mass_real, total_mass_real, medd, total_mdoteddington;
        double mass_holes, total_mass_holes, total_mdot;
        /*
        double mmax = 0, macc = 0;
        struct {
        	double val;
        	int rank;
        } in, out;
        */
        
        /* sum up numbers to print for summary of the BH step (blackholes.txt) */
        mdot = mass_holes = mass_real = medd = 0;
        for(bin = 0; bin < TIMEBINS; bin++)
        {
            if(TimeBinCount[bin])
            {
                mass_holes += TimeBin_BH_mass[bin];
                mass_real += TimeBin_BH_dynamicalmass[bin];
                mdot += TimeBin_BH_Mdot[bin];
                medd += TimeBin_BH_Medd[bin];
                /*
                if (TimeBin_BH_maxmass[bin]>mmax)
                {
                	mmax = TimeBin_BH_maxmass[bin];
                	macc = TimeBin_BH_maxmacc[bin];
                }
                */
            }
        }
        /*
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        in.val = mmax;
        in.rank = myrank;
        */
        MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&medd, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        /*
        MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        if(ThisTask==0) MPI_Bcast(&out, 1, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);
        if(myrank == out.rank)
        {
        	MPI_Bcast(&macc, 1, MPI_DOUBLE, myrank, MPI_COMM_WORLD);
        }
        */
        
        if((ThisTask == 0) && (total_mass_real > 0))
        {
            /* convert to solar masses per yr */
            mdot_in_msun_per_year = total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
            total_mdoteddington *= 1.0 / bh_eddington_mdot((double) All.TotBHs);
            //total_mdoteddington = total_mdot / bh_eddington_mdot(total_mass_holes);
            //if(All.TotBHs > 0) total_mdoteddington /= All.TotBHs;
            fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
                    All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
                    total_mass_real, total_mdoteddington);
        }
//#ifdef IO_REDUCED_MODE  DAA-IO: this is redundant
//        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
//#endif
//        {fflush(FdBlackHoles);} 
        fflush(FdBlackHoles);

//#ifndef IO_REDUCED_MODE   DAA-IO: BH_OUTPUT_MOREINFO overrides IO_REDUCED_MODE
#if !defined(IO_REDUCED_MODE) || defined(BH_OUTPUT_MOREINFO)
        fflush(FdBlackHolesDetails);
#ifdef BH_OUTPUT_GASSWALLOW
        fflush(FdBhSwallowDetails);
#endif
#ifdef BH_OUTPUT_MOREINFO
        fflush(FdBhMergerDetails);
#ifdef BH_WIND_KICK
        fflush(FdBhWindDetails);
#endif
#endif
#endif
    }
    myfree(BlackholeTempInfo);
}




/* simple routine to add quantities to BlackholeTempInfo */
void out2particle_blackhole(struct blackhole_temp_particle_data *out, int target, int mode)
{
    int k;
    ASSIGN_ADD(BlackholeTempInfo[target].BH_InternalEnergy,out->BH_InternalEnergy,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Mstar_in_Kernel,out->Mstar_in_Kernel,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].Mgas_in_Kernel,out->Mgas_in_Kernel,mode);
#ifdef BH_DEBUG
    ASSIGN_ADD(BlackholeTempInfo[target].Ngas_in_Kernel,out->Ngas_in_Kernel,mode);
#endif
    ASSIGN_ADD(BlackholeTempInfo[target].Malt_in_Kernel,out->Malt_in_Kernel,mode);
    //ASSIGN_ADD(BlackholeTempInfo[target].Sfr_in_Kernel,out->Sfr_in_Kernel,mode);
    for(k=0;k<3;k++)
    {
        ASSIGN_ADD(BlackholeTempInfo[target].Jgas_in_Kernel[k],out->Jgas_in_Kernel[k],mode);
        ASSIGN_ADD(BlackholeTempInfo[target].Jstar_in_Kernel[k],out->Jstar_in_Kernel[k],mode);
        ASSIGN_ADD(BlackholeTempInfo[target].Jalt_in_Kernel[k],out->Jalt_in_Kernel[k],mode);
    }
    
#ifdef GW_TIME
#ifdef OTF_BULGE
	ASSIGN_ADD(BlackholeTempInfo[target].star_rms_vel,out->star_rms_vel,mode);
#ifdef BH_DEBUG
	ASSIGN_ADD(BlackholeTempInfo[target].Nstar_in_Kernel,out->Nstar_in_Kernel,mode);
#endif
#endif
#endif
    
#ifdef BH_DYNFRICTION
    ASSIGN_ADD(BlackholeTempInfo[target].DF_rms_vel,out->DF_rms_vel,mode);
    for(k=0;k<3;k++)
    {
    	ASSIGN_ADD(BlackholeTempInfo[target].DF_mean_vel[k],out->DF_mean_vel[k],mode);
#ifdef ROMULUS_DF
		ASSIGN_ADD(BlackholeTempInfo[target].vCOM[k],out->vCOM[k],mode);
#endif
    }
    if(mode==0) {BlackholeTempInfo[target].DF_mmax_particles = out->DF_mmax_particles;}
        else {if(out->DF_mmax_particles > BlackholeTempInfo[target].DF_mmax_particles) {BlackholeTempInfo[target].DF_mmax_particles = out->DF_mmax_particles;}}
#endif

#ifdef MOD_BONDI
	ASSIGN_ADD(BlackholeTempInfo[target].vtheta_norm,out->vtheta_norm,mode);
	ASSIGN_ADD(BlackholeTempInfo[target].vbulk,out->vbulk,mode);
	for(k=0;k<3;k++)
		ASSIGN_ADD(BlackholeTempInfo[target].Jmb_in_Kernel[k],out->Jmb_in_Kernel[k],mode);
#endif

#if defined(BH_BONDI) || defined(BH_DRAG) || (BH_GRAVACCRETION >= 5)
    for(k=0;k<3;k++) {ASSIGN_ADD(BlackholeTempInfo[target].BH_SurroundingGasVel[k],out->BH_SurroundingGasVel[k],mode);}
#endif

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(BH_WIND_CONTINUOUS)
    for(k=0;k<3;k++) {ASSIGN_ADD(BlackholeTempInfo[target].GradRho_in_Kernel[k],out->GradRho_in_Kernel[k],mode);}
#endif
#if (BH_GRAVACCRETION == 8)
    ASSIGN_ADD(BlackholeTempInfo[target].hubber_mdot_bondi_limiter,out->hubber_mdot_bondi_limiter,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].hubber_mdot_vr_estimator,out->hubber_mdot_vr_estimator,mode);
    ASSIGN_ADD(BlackholeTempInfo[target].hubber_mdot_disk_estimator,out->hubber_mdot_disk_estimator,mode);
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    ASSIGN_ADD(BlackholeTempInfo[target].mass_to_swallow_edd, out->mass_to_swallow_edd, mode);
#endif

}

int ngb_treefind_blackhole_range(MyDouble searchcenter[3], MyFloat hsml1, MyFloat hsml2, int target, int *startnode, int mode, int *nexport, int *nsend_local)
{
    int numngb, no, p, task, nexport_save;
    struct NODE *current;
    MyDouble dx, dy, dz, dist1, dist2, r2;
    
#ifdef PERIODIC
    MyDouble xtmp;
#endif
    nexport_save = *nexport;
    
    numngb = 0;
    no = *startnode;
    
    while(no >= 0)
    {
        if(no < All.MaxPart)	/* single particle */
        {
            p = no;
            no = Nextnode[no];
            
            /* make sure we get all the particle types we need */
#ifdef MOD_BONDI
            if(P[p].Type != 0)
                continue;
#else            
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVACCRETION) || defined(BH_GRAVCAPTURE_NONGAS) || defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_DYNFRICTION)
            if(P[p].Type < 0)
                continue;
#else
            if(P[p].Type != 0 && P[p].Type != 3)
                continue;
#endif
#endif
            dist1 = hsml1;
            dist2 = hsml2;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dx > dist2 || dx < dist1)
                continue;
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dy > dist2 || dy < dist1)
                continue;
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dz > dist1 || dz < dist1)
                continue;
            r2 = dx * dx + dy * dy + dz * dz;
            if(r2 > dist2 * dist2 || r2 < dist1 * dist1)
                continue;
            
            Ngblist[numngb++] = p;
        }
        else
        {
            if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
            {
                if(mode == 1)
                    endrun(12312);
                
                if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                {
                    Exportflag[task] = target;
                    Exportnodecount[task] = NODELISTLENGTH;
                }
                
                if(Exportnodecount[task] == NODELISTLENGTH)
                {
                    if(*nexport >= All.BunchSize)
                    {
                        *nexport = nexport_save;
                        for(task = 0; task < NTask; task++)
                            nsend_local[task] = 0;
                        for(no = 0; no < nexport_save; no++)
                            nsend_local[DataIndexTable[no].Task]++;
                        return -1;
                    }
                    Exportnodecount[task] = 0;
                    Exportindex[task] = *nexport;
                    DataIndexTable[*nexport].Task = task;
                    DataIndexTable[*nexport].Index = target;
                    DataIndexTable[*nexport].IndexGet = *nexport;
                    *nexport = *nexport + 1;
                    nsend_local[task]++;
                }
                
                DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                
                if(Exportnodecount[task] < NODELISTLENGTH)
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                
                no = Nextnode[no - MaxNodes];
                continue;
            }
            
            current = &Nodes[no];
            
            if(mode == 1)
            {
                if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                {
                    *startnode = -1;
                    return numngb;
                }
            }
            
            no = current->u.d.sibling;	/* in case the node can be discarded */
            
            dist2 = hsml2 + 0.5 * current->len;;
            dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dx > dist2)
                continue;
            dy = NGB_PERIODIC_LONG_Y(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dy > dist2)
                continue;
            dz = NGB_PERIODIC_LONG_Z(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dz > dist2)
                continue;
            /* now test against the minimal sphere enclosing everything */
            dist2 += FACT1 * current->len;
            if(dx * dx + dy * dy + dz * dz > dist2 * dist2)
                continue;
            
            no = current->u.d.nextnode;	/* ok, we need to open the node */
        }
    }
    
    *startnode = -1;
    return numngb;
}

int ngb_treefind_blackhole(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local)
{
    int numngb, no, p, task, nexport_save;
    struct NODE *current;
    MyDouble dx, dy, dz, dist;
    
#ifdef PERIODIC
    MyDouble xtmp;
#endif
    nexport_save = *nexport;
    
    numngb = 0;
    no = *startnode;
    
    while(no >= 0)
    {
        if(no < All.MaxPart)	/* single particle */
        {
            p = no;
            no = Nextnode[no];
            
            /* make sure we get all the particle types we need */
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVACCRETION) || defined(BH_GRAVCAPTURE_NONGAS) || defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_DYNFRICTION)
            if(P[p].Type < 0)
                continue;
#else
            if(P[p].Type != 0 && P[p].Type != 3)
                continue;
#endif
            dist = hsml;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
            if(dz > dist)
                continue;
            if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;
            
            Ngblist[numngb++] = p;
        }
        else
        {
            if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
            {
                if(mode == 1)
                    endrun(12312);
                
                if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                {
                    Exportflag[task] = target;
                    Exportnodecount[task] = NODELISTLENGTH;
                }
                
                if(Exportnodecount[task] == NODELISTLENGTH)
                {
                    if(*nexport >= All.BunchSize)
                    {
                        *nexport = nexport_save;
                        for(task = 0; task < NTask; task++)
                            nsend_local[task] = 0;
                        for(no = 0; no < nexport_save; no++)
                            nsend_local[DataIndexTable[no].Task]++;
                        return -1;
                    }
                    Exportnodecount[task] = 0;
                    Exportindex[task] = *nexport;
                    DataIndexTable[*nexport].Task = task;
                    DataIndexTable[*nexport].Index = target;
                    DataIndexTable[*nexport].IndexGet = *nexport;
                    *nexport = *nexport + 1;
                    nsend_local[task]++;
                }
                
                DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                
                if(Exportnodecount[task] < NODELISTLENGTH)
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                
                no = Nextnode[no - MaxNodes];
                continue;
            }
            
            current = &Nodes[no];
            
            if(mode == 1)
            {
                if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                {
                    *startnode = -1;
                    return numngb;
                }
            }
            
            no = current->u.d.sibling;	/* in case the node can be discarded */
            
            dist = hsml + 0.5 * current->len;;
            dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
            if(dz > dist)
                continue;
            /* now test against the minimal sphere enclosing everything */
            dist += FACT1 * current->len;
            if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;
            
            no = current->u.d.nextnode;	/* ok, we need to open the node */
        }
    }
    
    *startnode = -1;
    return numngb;
}

