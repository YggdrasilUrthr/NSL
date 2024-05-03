/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"
// Std vector
#include <vector>

int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw;
double vtail, ptail, bin_size, nbins, sd, max_r;
double walker[m_props];
std::vector<double> gr;

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
std::vector<double> gr_blk_av;
double glob_av[m_props], glob_av2[m_props];
std::vector<double> gr_glob_av, gr_glob_av2;
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp, stima_gri;
double err_pot, err_press, err_kin, err_etot, err_temp;
std::vector<double> err_gr;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nblk, restart;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
