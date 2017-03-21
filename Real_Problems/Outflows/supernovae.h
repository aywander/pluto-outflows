//
// Created by Alexander Y. Wagner on 6/21/16.
//

#ifndef PLUTO_SUPERNOVAE_H
#define PLUTO_SUPERNOVAE_H

#include "pluto_usr.h"

void SetSupernovaePhysics();


typedef struct{

    double pow;                     // Supernovae energy
    double dur;                     // Supernovae duration
    double rad;                     // Starburst radius
    double hgt;                     // Supernovae/Starburst height
    double en1;                     // Single Supernova energy
    double ene;                     // Starburst total energy
    //int is_starburst;               // Single Supernovae or Starburst
    long num_tot;                   // Number of supernovae in this timestep
    long num_left;                  // Number of supernovae remaining
    long num_dt;                    // Number of supernovae in this timestep
    double x1[NSTARS_MAX];          // x1 position of supernovae
    double x2[NSTARS_MAX];          // x2 position of supernovae
    double x3[NSTARS_MAX];          // x3 position of supernovae

} Supernovae;

extern Supernovae sn;

void StrewSupernovae();

void InjectSupernovae(double *result, const double x1, const double x2, const double x3, double sn_inj);

#endif //PLUTO_SUPERNOVAE_H
