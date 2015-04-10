

double Lorentz2Vel(const double lorentz);
double Vel2Lorentz(const double vel);

double Profile(const double x1, const double x2, const double x3);

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
  double fill[4];
} VarNorm;

typedef struct {
  double ang;               // Half opening angle.
  double rad;               // Radius of cone.
  double dir;               // Angle direction of cone.
  double dbh;               // Location of rotation axis (BH) relative to (0,0,0).
                            // Can be positive or negative.
  double cbh;               // Height of cone from cap up to rotation axis.
                            // Mathematically always positive.
  double omg;               // Precession angular velocity.
  double phi;               // Precession starting angle.
  double sph;               // Radius of spherical region if INTERNAL_BOUNDARY = YES
  double orig;              // Height of beginning of domain as measured 
                            // from the origin along flowaxis.
                            // = g_domBeg[FLOWAXIS(IDIR,JDIR,KDIR)];.
                            // Can be positive or negative.
  double area;              // Area through which flux is calculated.
  double cone_height;       // Height of cone.
                            // Always positive.
  double cone_apex;         // Location of apex relative to (0,0,0)
                            // for a cone whose apex is on the flow axis.
                            // Can be used for a cone after RotateGrid2Nozzle transform.
                            // Can be positive or negative.
  int isfan;                // Is the nozzle a fan (conical) or bullet-shaped (parallel)
  int fill[3];             // Useless, just to make the struct a power of 2.
} Nozzle;


extern VarNorm vn;
extern double ini_cgs[32];
extern double ini_code[32];

extern Nozzle nozzle;

void OutflowPrimitives(double* out_primitives, 
                       const double x1, const double x2, const double x3);

void OutflowVelocity(double * out_primitives, double speed,
    const double x1, const double x2, const double x3);

void JetPrimitives(double* jet_primitives,
                   const double x1, const double x2, const double x3);

void UfoPrimitives(double* ufo_primitives, 
                   const double x1, const double x2, const double x3);


void HotHaloPrimitives(double* halo, 
                       const double x1, const double x2, const double x3);

int CloudCubePixel(int* el, const double x1,
                   const double x2,
                   const double x3);

void ReadFractalData();

void GetFractalData(double* cloud, 
    const double x1, const double x2, const double x3);

void CloudApodize(double* cloud, 
                  const double x1, const double x2, const double x3);

int CloudExtract(double* cloud,
                 const double* halo,
                 const int* pixel,
                 const double x1, const double x2, const double x3);

int CloudPrimitives(double* cloud,
                    const double x1, const double x2, const double x3);

void CloudVelocity(double* cloud,
                   const double x1, const double x2, const double x3);

int WarmTcrit(double* const warm);

void SetBaseNormalization();

void SetIniNormalization();

void SetNozzleConeGeometry();

void PrintInitData01(const double* ufo_primitives,
                     const double* halo_primitives);

int InNozzleRegion(const double x1, const double x2, const double x3);

int InNozzleSphere(const double x1, const double x2, const double x3);

int SphereIntersectsDomain(Grid *grid);


int RotateGrid2Nozzle(double const cx1, double const cx2, double const cx3,
                      double* cx1p, double* cx2p, double* cx3p);

int RotateNozzle2Grid(double const cx1, double const cx2, double const cx3,
                      double* cx1p, double* cx2p, double* cx3p);

