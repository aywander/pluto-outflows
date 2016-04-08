//---for multi clouds---//

struct InGrid {
    double *x1;
    double *x2;
    double *x3;
    int id_nx1;
    int id_nx2;
    int id_nx3;
};

struct cld_domain {
    double x1c;
    double x2c;
    double x3c;
    double v1;
    double v2;
    double v3;
};

void readgridfile(struct InGrid *);
void Init_multiclds(double *, double, double, double, struct cld_domain);
void Read_Multicld(char *);
int MultiCloudPrimitives(double*, const double, const double, const double, struct cld_domain);

double ran1(long int *);
double gasdev(long int *);
double Simpson_ext(double,double *,int);
void gen_cldlist(int,int);

