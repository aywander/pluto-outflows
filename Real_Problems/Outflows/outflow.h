//
// Created by Alexander Y. Wagner on 4/7/16.
//

#ifndef PLUTO_OUTFLOW_H
#define PLUTO_OUTFLOW_H

/* Structure for outflow state parameters. All are kept in code units. */
typedef struct {
    double pow;               // Power of outflow as given by input parameter
    double mdt;               // Mass related parameter of outflow as given by input parameter
    double spd;               // Speed related parameter of outflow as given by input parameter
    double prs;               // Pressure as calculated from pow, mdt, and speed
    double rho;               // Density as calculated from pow, mdt, and speed
    double eth;               // Enthalpy inj. rate as calculated from pow, mdt, and speed
    double kin;               // Initial kinetic energy inj. rate as calculated from pow, mdt, and speed
    int  is_on;               // Whether outflow is on or not
} OutflowState;

extern OutflowState os;

void OutflowPrimitives(double *out_primitives, double x1, double x2, double x3);

void OutflowVelocity(double *out_primitives, double speed, double x1, double x2, double x3);

void SetOutflowState(OutflowState *ofs);

void SetJetState(OutflowState *ofs);

void SetUfoState(OutflowState *ofs);

void OutflowOutput();

int FeedbackTrigger(double power, double speed);

#endif //PLUTO_OUTFLOW_H
