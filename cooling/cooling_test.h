#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/*
 * This file contains the definitions for the cooling.c routines
 *
 * This file was originally part of the GADGET3 code developed by
 *   Volker Springel (volker.springel@h-its.org). The code has been modified by
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer, int target);
void   InitCool(void);
void   InitCoolMemory(void);
void   IonizeParams(void);
void   IonizeParamsFunction(void);
void   IonizeParamsTable(void);
double INLINE_FUNC LogTemp(double u, double ne);
void   MakeCoolingTable(void);
void   ReadIonizeParams(char *fname);
void   SetZeroIonization(void);
void   TestCool(void);

void   find_abundances_and_rates(double logT, double rho, double *ne_guess, int target);
double convert_u_to_temp(double u, double rho, double *ne_guess, int target);
double CoolingRate(double logT, double rho, double *nelec, int target);
//double CoolingRateFromU(double u, double rho, double *ne_guess, int target, double dt);
double CoolingRateFromU(double u, double rho, double ne, int target, double dt);
double DoCooling(double u_old, double rho, double dt, double *ne_guess, int target);
double GetCoolingTime(double u_old, double rho,  double *ne_guess, int target);
//double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess, int target);
double coldgasfrac(int i, double T, double rho, int flag);
//added by JDJ for primordial cooling based on Volker code 
double hlcool(int i, double redshift, double temp, double xnh);
double hmcool(int i, double Tcmb, double temp, double xnh);
void P_Rates(double tempr,double xnh,double J_21,double del_r);
void P_Chem(int target,double dt_tot,double xn, double rho, double tempr,double J_21,double del_r);
//void P_Chem(int target,double dt_tot,double xn, double rho, double tempr,double J_21,double del_r);
double get_fshield(int target, double temp, double xn, double rho);
//added by JDJ for metal cooling
double JJ_metal_cooling(int i, double temp, double xnh);
double PE2_heating(double xnh);
double PE3_heating(double xnh);
double UVB_heating(double xnh, int target);
double lam_2level(double gamma_H_21,double gamma_H_12,double gamma_e_21,double gamma_e_12,double A21,double delta_E21,double ne,double nh,double ntot);

int JeansMinT(double *u, double rho, double *ne_guess, int target);
#ifdef GRACKLE
void InitGrackle(void);
double CallGrackle(double u_old, double rho, double dt, double *ne_guess, int target, int mode);
#endif
