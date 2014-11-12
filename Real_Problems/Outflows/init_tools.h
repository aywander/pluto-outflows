

double Lorentz2Vel(const double lorentz);
double Vel2Lorentz(const double vel);

double Profile(const double x1, const double x2, const double x3);
double Profile_sharp(const double x1, const double x2, const double x3);

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
} VarNorm;

extern VarNorm vn;

extern double ini_cgs[32];
extern double ini_code[32];


void OutflowPrimitives(double* out_primitives, 
    const double x1, const double x2, const double x3);

void OutflowVelocity(double * out_primitives, double speed,
    const double x1, const double x2, const double x3, const double th);

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

int WarmTcrit(double* const warm);

void SetBaseNormalization();

void SetIniNormalization();

void PrintInitData01(const double* ufo_primitives,
                            const double* halo_primitives);

int InNozzleBase(const double x1, const double x2, const double x3);

int InNozzleRegion(const double x1, const double x2, const double x3);

int RotateGrid2Nozzle(const double r, const double x, double* rp, double* xp);

int RotateNozzle2Grid(const double r, const double x, double* rp, double* xp);

int TranslateRotateGrid2Nozzle(const double r, const double x, double* rp, double* xp);

int TranslateRotateNozzle2Grid(const double r, const double x, double* rp, double* xp);


