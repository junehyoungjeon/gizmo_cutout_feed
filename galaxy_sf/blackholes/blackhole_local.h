/*! \file blackhole_local.h
 *  \brief comm structures to be used throughout blackhole routines
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

#ifdef BH_FB_CONE
#ifndef BH_FB_CONE_NMIN
#define BH_FB_CONE_NMIN 4
#endif
/* unit jet axis from kernel angular momentum; (0,0,0) if undefined.
   BH_FB_CONE_FLIP selects the -J pole instead of +J */
static inline void bh_fb_unit_axis(double *J, double *n)
{
#ifdef BH_FB_AXIS_X
    /* fixed axis in the simulation frame; J is ignored */
    double ax = BH_FB_AXIS_X, ay = BH_FB_AXIS_Y, az = BH_FB_AXIS_Z;
    double a = sqrt(ax*ax + ay*ay + az*az);
    if(a > 0) {n[0]=ax/a; n[1]=ay/a; n[2]=az/a;} else {n[0]=n[1]=n[2]=0;}
#else
    double a = sqrt(J[0]*J[0] + J[1]*J[1] + J[2]*J[2]);
#ifdef BH_FB_CONE_FLIP
    if(a > 0) {n[0]=-J[0]/a; n[1]=-J[1]/a; n[2]=-J[2]/a;} else {n[0]=n[1]=n[2]=0;}
#else
    if(a > 0) {n[0]=J[0]/a; n[1]=J[1]/a; n[2]=J[2]/a;} else {n[0]=n[1]=n[2]=0;}
#endif
#endif
}
/* c = SIGNED cos(angle between r_hat (BH -> gas) and the axis): one cone around +axis only */
static inline double bh_cone_weight(double c)
{
    if(c > 1) c = 1;
    if(c < -1) c = -1;
#ifdef BH_FB_CONE_GAUSS
    double th = acos(c), s = BH_FB_CONE * (M_PI/180.0);   /* th runs 0..pi, so the far hemisphere is suppressed automatically */
    return exp(-0.5*th*th/(s*s));
#else
    return (c >= cos(BH_FB_CONE * (M_PI/180.0))) ? 1.0 : 0.0;
#endif
}
#endif

/* quantities that pass IN to the 'blackhole_evaluate' routines */
static struct blackholedata_in
{
    int orig_index;
    int trans_index;
    
    MyDouble Pos[3];
    MyFloat Density;
    MyFloat Mdot;
    MyFloat Dt;
    MyFloat Hsml;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    MyFloat AGS_Hsml;
#endif   
    MyFloat Mass;
    MyFloat BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat BH_Mass_AlphaDisk;
#endif
    MyFloat Vel[3];
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS) || defined(BH_BAL_KICK_COLLIMATED)
    MyFloat Jgas_in_Kernel[3];
#endif
#ifdef BH_FB_CONE
    MyFloat BH_fb_axis[3];
    MyFloat BH_fb_sum[3];
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat BH_disk_hr;
    MyFloat BH_angle_weighted_kernel_sum;
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    MyFloat mass_to_swallow_edd;
#endif
#if defined(BH_DYNFRICTION) && defined(ROMULUS_DF)
	MyFloat DF_mean_vel[3];
	MyFloat vCOM[3];
#endif
#if defined(GALSF) || defined(BLACK_HOLES)
	MyIDType parentID;
#endif
#ifdef GW_TIME
	MyFloat StellarAge;
	MyFloat SF_dens;
#ifdef OTF_BULGE
	MyFloat star_rms_vel;
	MyFloat rho_star;
#endif
#endif
}
*BlackholeDataIn, *BlackholeDataGet;


/* quantities that pass OUT of the 'blackhole_evaluate' routines */
static struct blackholedata_out
{
    
    int orig_index;
    int trans_index;
    
    MyLongDouble Mass;
    MyLongDouble BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
    MyLongDouble BH_Mass_AlphaDisk;
#endif
    MyLongDouble accreted_Mass;
    MyLongDouble accreted_BH_Mass;
    MyLongDouble accreted_momentum[3];
#ifdef BH_COUNTPROGS
    int BH_CountProgs;
#endif
#ifdef GALSF
	MyIDType parentID;
	MyFloat SF_dens;
    MyFloat Accreted_Age;
#endif
#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
    MyFloat BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_FB_CONE
    MyFloat BH_fb_sum[3];
#endif
#ifdef BH_BUBBLES
    MyFloat BH_Mass_bubbles;
    MyFloat BH_Mass_ini;
#ifdef UNIFIED_FEEDBACK
    MyFloat BH_Mass_radio;
#endif
#endif
}
*BlackholeDataResult, *BlackholeDataOut;


