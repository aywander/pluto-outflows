#ifndef init_tools_h
#define init_tools_h

double Lorentz2Speed(const double lorentz);

double Speed2Lorentz(const double vel);


/* Regarding the fill arrays:
 *
 * For 64 bit linux system (typically)
 * type     nbytes
 * char     1
 * int      4
 * float    8
 * double   16
 * pointer  8
 * Tip: Use an array of the smallest type 
 * in the struct to create the fill. Then
 * you can just count in multiples of that
 * type.
 *
 * */

/* Structure containing normalization factors such that 
 * ... quantity*<..._norm> = ... quantity in cgs units */
typedef struct {
    double l_norm;
    double dens_norm;
    double v_norm;
    double t_norm;
    double power_norm;
    double eflux_norm;
    double eint_norm;
    double pres_norm;
    double area_norm;
    double temp_norm;
    double mdot_norm;
    double newton_norm;
    double pot_norm;
    double acc_norm;
    double n_norm;
    double m_norm;
} VarNorm;

extern VarNorm vn;


extern double ini_cgs[USER_DEF_PARAMETERS];
extern double ini_code[USER_DEF_PARAMETERS];


void SetBaseNormalization();

void SetIniNormalization();


void PrintInitData01();

void PrintGridStruct(Grid *grid, int show_for_rank, int k, int j, int i);


#endif
