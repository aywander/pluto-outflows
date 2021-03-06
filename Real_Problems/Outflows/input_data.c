/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Provide basic functionality for reading input data files.

  Collects a number of functions for opening, reading and interpolating
  initial conditions from user-supplied binary data file(s).
  The geometry and dimensions of the input grid can be different from 
  the actual grid employed by PLUTO, as long as the coordinate geometry
  transformation has been implemented.
  The input grid and data files should employ the same format and 
  conventions employed by PLUTO.
  Specifically:
  
  - Gridfile: coordinates should be written using the standard
              PLUTO 4 grid format.
  - Datafile: variables should be written in single or multiple
              binary data file using eighter single or double precision.
              The file extension must be ".flt" or ".dbl" for the former and
              the latter, respectively.
     
  InputDataOpen() initializes the inputData structure for a given
  input field (e.g. density) iby reading grid size and coordinates,
  geometry, precision, endianity, etc..

  InputDataInterpolate() is used to interpolate the given field from the
  input coordinates to the desired coordinate location using bi- or
  tri-linear interpolation to fill the data array.

  The input data is stored in a buffer by reading ::ID_NZ_MAX planes at
  a time to save computational memory.
  The value of ::ID_NZ_MAX can be changed from your personal \c definitions.h.
  This task is performed in the InputDataReadSlice() function.

  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Nov 17, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* AYW --
 * To use my geometry macros */
#include "pluto_usr.h"
/* -- AYW */

#define ID_NVAR_MAX 64
#ifndef ID_NZ_MAX
 #define ID_NZ_MAX   4  /**< Maximum size (in the 3rd dimension) of the
                             input buffer. */
#endif

typedef struct inputData_{
  char fname[64];
  size_t dsize;
  int nx1;
  int nx2; 
  int nx3;
  int geom;
  int swap_endian;
  int klast;        /**< The index of the last read plane */
  double *x1;
  double *x2;
  double *x3;
  double ***Vin;    /**< Input buffer array (containing at most ::ID_NZ_MAX
                         planes at a time */
  long int offset;
} inputData;

static inputData id_stack[ID_NVAR_MAX];


/* AYW --
 * To store input data extents */

double g_idBoxBeg[3];  /**< Lower limits of the input data domain. */
double g_idBoxEnd[3];  /**< Upper limits of the input data domain. */
double g_idnx1, g_idnx2, g_idnx3;
int g_idnvar, *g_idvarindx; /**< Number of variables */

/* -- AYW */


/* AYW -- Helper functions */

double shift_idx(double x, const double xmin, const double xmax,
                 const double id_x_beg, const double id_x_end);

void InputDataCoordTransformPoint(int indx, double *x1, double *x2, double *x3);

/* -- AYW */

/* ********************************************************************* */
int InputDataOpen(char *data_fname, char *grid_fname, char *endianity, int pos)
/*!
 * Initialize access to input data file by assigning values to 
 * grid-related information (geometry, number of points, etc...).
 * This function should be called just once for input-data initialization.
 *
 * \param [in] data_fname    the name of the input binary data file 
 * \param [in] grid_fname    the name of the corresponding grid file
 * \param [in] endianity     an input string ("little" or "big") giving
 *                           the byte-order of how the input data file
 *                           was originally written.
 *                           If an empty string is supplied, no change is
 *                           made.
 * \param [in] pos           an integer specifying the position of
 *                           the variable in the file.
 *
 *********************************************************************** */
{
  char *sub_str, sline[256];
  const char delimiters[] = " \t\r\f\n";
  int    i, ip, success;
  int    indx; /* The 1st available element in the id_stack structure */
  size_t str_length;
  double xl, xr;
  inputData *id;
  fpos_t file_pos;
  FILE *fp;

  print ("> Input data:\n\n");

/* --------------------------------------------------------------------- */
/*! 0. Find the 1st available (NULL pointer) element in the stack        */
/* --------------------------------------------------------------------- */

  indx = 0;
  while (id_stack[indx].Vin != NULL){
    indx++;
    if (indx == ID_NVAR_MAX){
      print ("! InputDataOpen(): max number of variables exceeded\n");
      print ("!                  indx = %d\n",indx);
      QUIT_PLUTO(1);
    }
  }

  print ("  Allocating memory for struct # %d\n",indx);
  id = id_stack + indx;

/* --------------------------------------------------------------------- */
/*! 1. Scan grid data file and try to determine the grid geometry 
       (\c id->geom). Search for tag "GEOMETRY:" and read the word that
       follows.                                                          */
/* --------------------------------------------------------------------- */

  fp = fopen(grid_fname,"r");
  if (fp == NULL){
    print ("! InputDataOpen(): grid file %s not found\n",grid_fname);
    QUIT_PLUTO(1);
  }
  success = 0;
  while(!success){
    fgets(sline,512,fp);
    sub_str = strtok(sline, delimiters);
    while (sub_str != NULL){
      if (!strcmp(sub_str,"GEOMETRY:")) {
        sub_str = strtok(NULL, delimiters);     
        success = 1;
        break;
      }
      sub_str = strtok(NULL, delimiters);     
    }  
  }
  
  if      (!strcmp(sub_str,"CARTESIAN"))   id->geom = CARTESIAN;
  else if (!strcmp(sub_str,"CYLINDRICAL")) id->geom = CYLINDRICAL;
  else if (!strcmp(sub_str,"POLAR"))       id->geom = POLAR;
  else if (!strcmp(sub_str,"SPHERICAL"))   id->geom = SPHERICAL;
  else{
    print ("! InputDataOpen(): unknown geometry\n");
    QUIT_PLUTO(1);
  }
    
  print ("  Input grid file:       %s\n", grid_fname);
  print ("  Input grid geometry:   %s\n", sub_str);

/* --------------------------------------------------------------------- */
/*! 2. Move file pointer to the first line of grid.out that does not
       begin with a "\c #".                                              */
/* --------------------------------------------------------------------- */

  success = 0;
  while(!success){
    fgetpos(fp, &file_pos);
    fgets(sline,512,fp);
    if (sline[0] != '#') success = 1;
  }

  fsetpos(fp, &file_pos);
  
/* --------------------------------------------------------------------- */
/*! 3. Read number of points, allocate grid arrays and store input
       grid coordinates into structure members \c id->nx1, \c id->x1,
       etc...                                                            */
/* --------------------------------------------------------------------- */
   
    /* AYW -- 2013-04-19 12:08 JST
    * Store read data extents in global arrays */

  fscanf (fp,"%d \n",&(id->nx1));
  id->x1 = ARRAY_1D(id->nx1, double);
  for (i = 0; i < id->nx1; i++){
    fscanf(fp,"%d  %lf %lf\n", &ip, &xl, &xr);
    id->x1[i] = 0.5*(xl + xr);
        if (i == 0)           g_idBoxBeg[IDIR] = xl;   // AYW added
        if (i == id->nx1 - 1) g_idBoxEnd[IDIR] = xr;   // AYW added
  }

  fscanf (fp,"%d \n",&(id->nx2));
  id->x2 = ARRAY_1D(id->nx2, double);
  for (i = 0; i < id->nx2; i++){
    fscanf(fp,"%d  %lf %lf\n", &ip, &xl, &xr);
    id->x2[i] = 0.5*(xl + xr);
      if (i == 0)           g_idBoxBeg[JDIR] = xl;   // AYW added
      if (i == id->nx2 - 1) g_idBoxEnd[JDIR] = xr;   // AYW added
  }

  fscanf (fp,"%d \n",&(id->nx3));
  id->x3 = ARRAY_1D(id->nx3, double);
  for (i = 0; i < id->nx3; i++){
    fscanf(fp,"%d  %lf %lf\n", &ip, &xl, &xr);
    id->x3[i] = 0.5*(xl + xr);
      if (i == 0)          g_idBoxBeg[KDIR] = xl;   // AYW added
      if (i == id->nx3 - 1) g_idBoxEnd[KDIR] = xr;   // AYW added
  }
  fclose(fp);

/* -- reset grid with 1 point -- */

  if (id->nx1 == 1) id->x1[0] = 0.0;
  if (id->nx2 == 1) id->x2[0] = 0.0;
  if (id->nx3 == 1) id->x3[0] = 0.0;

    /* AYW --
    * Store read data number of points in global vars */
    g_idnx1 = id->nx1;
    g_idnx2 = id->nx2;
    g_idnx3 = id->nx3;
    /* -- AYW */

  print ("  Input grid extension:  x1 = [%12.3e, %12.3e] (%d points)\n",
             id->x1[0], id->x1[id->nx1-1], id->nx1);
  print ("\t\t\t x2 = [%12.3e, %12.3e] (%d points)\n",
             id->x2[0], id->x2[id->nx2-1], id->nx2);
  print ("\t\t\t x3 = [%12.3e, %12.3e] (%d points)\n",
             id->x3[0], id->x3[id->nx3-1], id->nx3);
  
/* ---------------------------------------------------------- */
/*! 4. Check if binary data file exists and allocate memory
       buffer \c id->Vin                                      */
/* ---------------------------------------------------------- */

  char   ext[] = "   ";

  sprintf (id->fname,"%s",data_fname);
  fp = fopen(id->fname, "rb");
  if (fp == NULL){
    print ("! InputDataOpen(): file %s does not exist\n", data_fname);
    QUIT_PLUTO(1);
  }
  fclose(fp);

  id->Vin = ARRAY_3D(ID_NZ_MAX, id->nx2, id->nx1, double);

/* ---------------------------------------------------------- */
/*! 5. Set endianity (\c id->swap_endian)                     */
/* ---------------------------------------------------------- */
 
  if ( (!strcmp(endianity,"big")    &&  IsLittleEndian()) ||
       (!strcmp(endianity,"little") && !IsLittleEndian())) {
    id->swap_endian = YES;
  }else{
    id->swap_endian = NO;
  }
  
  print ("  Input data file:       %s (endianity: %s) \n", 
           data_fname, endianity);
  
/* ---------------------------------------------------------- */
/*! 6. Set data size (\c id->dsize) by looking at the file
       extension (\c dbl or \c flt).                          */
/* ---------------------------------------------------------- */

  str_length = strlen(data_fname);
  for (i = 0; i < 3; i++) ext[i] = data_fname[str_length-3+i];

  if (!strcmp(ext,"dbl")){
    print ("  Precision:             (double)\n");
    id->dsize = sizeof(double);
  } else if (!strcmp(ext,"flt")) {
    print ("  Precision:\t\t single\n");
    id->dsize = sizeof(float);
  } else {
    print ("! InputDataRead: unsupported data type '%s'\n",ext);
    QUIT_PLUTO(1);
  }

/* ---------------------------------------------------------- */
/*! 7. Compute offset (\c id->offset in bytes) and
       initialize \c id->klast (last read plane) to -1        */
/* ---------------------------------------------------------- */

  id->offset = (long)pos*id->dsize*id->nx1*id->nx2*id->nx3;
  id->klast  = -1;

  print ("  offset = %ld\n",id->offset);
  print ("\n");

  return indx;  /* The index of the id_stack[] array */
}


/* ********************************************************************* */
void InputDataClose(int indx)
/*!
 * Free memory and reset structure.
 *********************************************************************** */
{
  inputData *id = id_stack + indx;
  FreeArray1D((void *)id->x1);
  FreeArray1D((void *)id->x2);
  FreeArray1D((void *)id->x3);
  FreeArray3D((void *)id->Vin);
  id->Vin = NULL;
}
/* ********************************************************************* */
double InputDataInterpolate(int indx, double x1, double x2, double x3)
/*!
 * Perform bi- or tri-linear interpolation on external
 * dataset to compute vs[] at the given point {x1,x2,x3}.
 *
 * \param [in] indx   input data array element (handle)
 * \param [in] x1     coordinate point at which at interpolates are desired
 * \param [in] x2     coordinate point at which at interpolates are desired
 * \param [in] x3     coordinate point at which at interpolates are desired
 *
 * \return The interpolated value.
 *
 * The function performs the following tasks. 
 *********************************************************************** */
{
  int il = 0, jl = 0, kl = 0;
  int ih, jh, kh;
  int im, jm, km;
  int i, j, k, nv;
  float  uflt;
  double udbl;
  double xx, yy, zz, v;
  double **Vlo, **Vhi;
  inputData *id = id_stack + indx;
  static FILE *fp;

/* --------------------------------------------------------------------- */
/*! - Convert PLUTO coordinates to input grid geometry if necessary.     */
/* --------------------------------------------------------------------- */

     /* AYW -- My version */
    InputDataCoordTransformPoint(indx, &x1, &x2, &x3);
    /* -- AYW */

//#if GEOMETRY == CARTESIAN
//  if (id->geom == GEOMETRY) {
//
//   /* same coordinate system: nothing to do */
//
//  }else if (id->geom == CYLINDRICAL) {
//    double R, z, phi;
//    R   = sqrt(x1*x1 + x2*x2);
//    phi = atan2(x2,x1);
//    if (phi < 0.0) phi += 2.0*CONST_PI;
//    z   = x3;
//
//    x1 = R; x2 = z; x3 = phi;
//  }else if (id->geom == POLAR) {
//    double R, phi, z;
//    R   = sqrt(x1*x1 + x2*x2);
//    phi = atan2(x2,x1);
//    if (phi < 0.0) phi += 2.0*CONST_PI;
//    z   = x3;
//
//    x1 = R; x2 = phi; x3 = z;
//  }else if (id->geom == SPHERICAL){
//    double r, theta, phi;
//    r     = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
//    r     = sqrt(r);
//    theta = acos(x3/r);
//    phi   = atan2(x2,x1);
//    if (phi   < 0.0) phi   += 2.0*CONST_PI;
//    if (theta < 0.0) theta += 2.0*CONST_PI;
//
//    x1 = r; x2 = theta; x3 = phi;
//  }else{
//    print ("! InputDataInterpolate(): invalid or unsupported coordinate transformation.\n");
//    QUIT_PLUTO(1);
//  }
//#elif GEOMETRY == CYLINDRICAL
//  if (id->geom == GEOMETRY) {
//
//    /* same coordinate system: nothing to do */
//
//  }else if (id->geom == SPHERICAL) {
//    double r, theta, phi;
//    r     = D_EXPAND(x1*x1, + x2*x2, + 0.0);
//    r     = sqrt(r);
//    theta = acos(x2/r);
//    phi   = 0.0;
//    if (theta < 0.0) theta += 2.0*CONST_PI;
//
//    x1 = r; x2 = theta; x3 = phi;
//  }else{
//    print ("! InputDataInterpolate(): invalid or unsupported coordinate transformation.\n");
//    QUIT_PLUTO(1);
//  }
//#elif GEOMETRY == POLAR
//  if (id->geom == GEOMETRY) {
//
//    /* same coordinate system: nothing to do */
//
//  }else if (id->geom == CARTESIAN) {
//    double x, y, z;
//    x = x1*cos(x2);
//    y = x1*sin(x2);
//    z = x3;
//
//    x1 = x; x2 = y; x3 = z;
//  }else{
//    print ("! InputDataInterpolate(): invalid or unsupported coordinate transformation.\n");
//    QUIT_PLUTO(1);
//  }
//#elif GEOMETRY == SPHERICAL
//  if (id->geom == GEOMETRY) {
//
//    /* same coordinate system: nothing to do */
//
//  }else if (id->geom == CARTESIAN) {
//    double x, y, z;
//    x = x1*sin(x2)*cos(x3);
//    y = x1*sin(x2)*sin(x3);
//    z = x1*cos(x2);
//
//    x1 = x; x2 = y; x3 = z;
//  }else{
//    print ("! InputDataInterpolate(): invalid or unsupported coordinate transformation.\n");
//    QUIT_PLUTO(1);
//  }
//#endif


/* --------------------------------------------------------------------- */
/*! - Make sure point (x1,x2,x3) does not fall outside input grid range. 
      Limit to input grid edge otherwise.                                */
/* --------------------------------------------------------------------- */


/* AYW -- DM's periodic input data */
// TODO: Change this to INPUT_DATA_TESSELATE
   #if (CLOUD_REPEAT == YES)
    /*
     DM 23 feb, 2015: Check if PLUTO grid point is outside input grid.
    If yes, then check if PLUTO grid is outside maximum zone of clouds
    defined in definitions_usr.h. If no, then fold PLUTO grid back on
    to appropriate location in input grid, assuming periodic BC.

    XL, XR are domain limits of input grid

    if x < XL
    |-----------|------------|--------------|
          x     XL           XR
          -------------->
    x' = x + (int((x - xL) / Delta) + 1) * Delta , Delta = XR - XL

    if x > XR
    |-----------|------------|--------------|
                XL           XR        x
                        <--------------
    x' = x - (int((x - xR) / Delta) + 1) * Delta

    Rationale: Translate input grid by n steps such that PLUTO grid point
    falls within translated input grid limits
    i.e. XR' = XR + n * Delta, n = int((x - XR) / Delta) + 1
    New x' = XR - (XR' - x) = x - n * Delta

    Same for lower boundary, results in a difference in sign.

    *** To fix: if PLUTO grid is larger than max domain defined in definitions_usr.h
    then PLUTO grid is given value eq to boundary value of input grid.
    */

    /* DM -- */
      D_EXPAND(x1 = shift_idx(x1, g_inputParam[PAR_WX1L], g_inputParam[PAR_WX1H], id->x1[0], id->x1[id_nx1 - 1]);,
               x2 = shift_idx(x2, g_inputParam[PAR_WX2L], g_inputParam[PAR_WX2H], id->x2[0], id->x2[id_nx2 - 1]);,
               x3 = shift_idx(x3, g_inputParam[PAR_WX3L], g_inputParam[PAR_WX3H], id->x3[0], id->x3[id_nx3 - 1]););
    /* -- DM */

#else
  /* Default PLUTO */
  D_EXPAND(if      (x1 < id->x1[0])         x1 = id->x1[0];
           else if (x1 > id->x1[id->nx1-1]) x1 = id->x1[id->nx1-1];  ,
           
           if      (x2 < id->x2[0])         x2 = id->x2[0];
           else if (x2 > id->x2[id->nx2-1]) x2 = id->x2[id->nx2-1];  ,
           
           if      (x3 < id->x3[0])         x3 = id->x3[0];
           else if (x3 > id->x3[id->nx3-1]) x3 = id->x3[id->nx3-1]; )
#endif
/* -- AYW */

/* --------------------------------------------------------------------- */
/*! - Use table lookup by binary search to  find the indices 
      il, jl and kl such that grid points of PLUTO fall between 
      [il, il+1], [jl, jl+1], [kl, kl+1].                                */
/* --------------------------------------------------------------------- */

  il = 0;
  ih = id->nx1 - 1;
  while (il != (ih-1)){
    im = (il+ih)/2;
    if   (x1 <= id->x1[im]) ih = im;   
    else                    il = im;
  }
  
  if (id->nx2 > 1){
    jl = 0;
    jh = id->nx2 - 1;
    while (jl != (jh-1)){
      jm = (jl+jh)/2;
      if (x2 <= id->x2[jm]) jh = jm;   
      else                  jl = jm;
    }
  }

  if (id->nx3 > 1){
    kl = 0;
    kh = id->nx3 - 1;
    while (kl != (kh - 1)){
      km = (kl+kh)/2;
      if (x3 <= id->x3[km]) kh = km;   
      else                  kl = km;
    }
  }

/* --------------------------------------------------------------------- */
/*! - Define normalized coordinates between [0,1]:
      - x[il+1] < x1[i] < x[il+1] ==> 0 < xx < 1
      - y[jl+1] < x2[j] < y[jl+1] ==> 0 < yy < 1
      - z[kl+1] < x3[k] < z[kl+1] ==> 0 < zz < 1                            */
/* --------------------------------------------------------------------- */

  xx = yy = zz = 0.0; /* initialize normalized coordinates */

  if (id->nx1 > 1) xx = (x1 - id->x1[il])/(id->x1[il+1] - id->x1[il]);  
  if (id->nx2 > 1) yy = (x2 - id->x2[jl])/(id->x2[jl+1] - id->x2[jl]);  
  if (id->nx3 > 1) zz = (x3 - id->x3[kl])/(id->x3[kl+1] - id->x3[kl]);

/* --------------------------------------------------------------------- */
/*! - Read data from disk (1 or 2 planes)                                */
/* --------------------------------------------------------------------- */


/*
print ("@@ Interpolation required at %f %f %f\n",x1,x2,x3);
print ("@@ Closest input coord.      %f %f %f\n",
        id->x1[il], id->x2[jl], id->x3[kl]);
print ("@@ Coord indices:            %d %d %d\n",il,jl,kl);
print ("@@ id->klast                 %d\n",id->klast);

static long int ncall=0;

  if ( (kl >= id->klast + ID_NZ_MAX - 1) || (id->klast == -1) ){
print ("@@ kl = %d, klast = %d; ncall = %d\n",kl,id->klast,ncall);
    InputDataReadSlice(indx, kl);
ncall = 0;
  }
ncall++;
*/

/* AYW -- There was a bug here in that if id->nx3 - kl < ID_NZ_MAX the slice is still read,
 * even though the InputDataReadSlice routine would try to read slices at kl > id->nx3,
 * resulting in read errors.
 */
  if ( (kl >= id->klast + ID_NZ_MAX - 1) || (id->klast == -1) ){
    InputDataReadSlice(indx, kl);
  }

/* --------------------------------------------------------------------- */
/*! - Perform bi- or tri-linear interpolation.                           */
/* --------------------------------------------------------------------- */

  Vlo = id->Vin[kl - id->klast];
  Vhi = id->Vin[kl - id->klast+1];

  v =   Vlo[jl][il]*(1.0 - xx)*(1.0 - yy)*(1.0 - zz)  /* [0] is kl */
      + Vlo[jl][il+1]*xx*(1.0 - yy)*(1.0 - zz);
  if (id->nx2 > 1){
    v +=   Vlo[jl+1][il]*(1.0 - xx)*yy*(1.0 - zz)
         + Vlo[jl+1][il+1]*xx*yy*(1.0 - zz);
  }
  if (id->nx3 > 1){
    v +=   Vhi[jl][il]*(1.0 - xx)*(1.0 - yy)*zz
         + Vhi[jl][il+1]*xx*(1.0 - yy)*zz
         + Vhi[jl+1][il]*(1.0 - xx)*yy*zz
         + Vhi[jl+1][il+1]*xx*yy*zz;
  }

  return v;
}

/* ********************************************************************* */
void InputDataReadSlice(int indx, int kslice)
/*! 
 * Read ::ID_NZ_MAX slices from disk starting at the kslice vertical
 * index.
 *
 * \param [in] indx    the structure index (file handle)
 * \param [in] kslice  the starting slice
 *
 *********************************************************************** */
{
  int i,j,k, kmax;
  long int offset;
  double udbl;
  float  uflt;
  inputData *id = id_stack + indx;
  FILE *fp;

/* ----------------------------------------------------
   1. Compute offset.
      Here id->offset (in bytes) is the location of the
      variable in the file, while the second number is
      the vertical slice we want to read.
      Seek from beginning of the file.
   ---------------------------------------------------- */

  offset = id->offset + kslice*id->dsize*id->nx1*id->nx2;

  fp = fopen(id->fname,"rb");
  fseek(fp, offset, SEEK_SET);

/* -----------------------------------------------------
   2. Read binary data at specified position.
   ----------------------------------------------------- */

/* AYW -- There is a bug here in that if id->nx3 - kl < ID_NZ_MAX the slice is still read,
 * even though the InputDataReadSlice routine would try to read slices at kl > id->nx3,
 * resulting in read errors.
 */
//  kmax = MIN(id->nx3,ID_NZ_MAX);
  kmax = MIN(id->nx3 - kslice, ID_NZ_MAX);
  /* -- AYW */
  if (id->dsize == sizeof(double)){

    for (k = 0; k < kmax   ; k++){   /* Read at most 2 planes */  
    for (j = 0; j < id->nx2; j++){
    for (i = 0; i < id->nx1; i++){
      if (fread (&udbl, id->dsize, 1, fp) != 1){  
        print ("! InputDataReadSlice(): error reading data (indx = %d)\n",indx);
        break;
      }
      if (id->swap_endian) SWAP_VAR(udbl);
      id->Vin[k][j][i] = udbl;
    }}}

   }else{

    for (k = 0; k < kmax   ; k++){   /* Read at most 2 planes */  
    for (j = 0; j < id->nx2; j++){
    for (i = 0; i < id->nx1; i++){
      if (fread (&uflt, id->dsize, 1, fp) != 1){
        print ("! InputDataReadSlice(): error reading data (indx = %d)\n",indx);
        break;
      }
      if (id->swap_endian) SWAP_VAR(uflt);
      id->Vin[k][j][i] = (double)uflt;
    }}}
  }

/* -- Update last successfully read slice -- */

  id->klast = kslice;
  fclose(fp);
/*
print ("@@@ Offset = %d\n",offset);
for (k = 0; k < kmax   ; k++){
for (j = 0; j < id->nx2; j++){
for (i = 0; i < id->nx1; i++){
  print ("@@@ Input value (%d %d %d) = %f\n",i,j,k,id->Vin[k][j][i]);
}}}
*/
}


/* AYW --
 * This function shifts a coordinate by a multiple of the input cube width
 * This is used for INPUT_DATA_TESSELATE */
double shift_idx(double x, const double xmin, const double xmax,
                 const double id_x_beg, const double id_x_end) {

    if (x < id_x_beg || x > id_x_end) {
        if (x > xmin && x < xmax) {

            double delta, deltax, intpart, fractpart;
            delta = id_x_end - id_x_beg;
            deltax = x - (x < id_x_beg ? id_x_beg : id_x_end);
            fractpart = modf(fabs(deltax) / delta, &intpart);
            x = (x < id_x_beg ? x + (intpart + 1.0) * delta : x - (intpart + 1.0) * delta);
        }

        if (x < xmin) x = xmin;
        if (x > xmax) x = xmax;
    }

    return x;
}
/* -- AYW */


/* AYW --
 * My version of transforming coordinates. */

/* ********************************************************************* */
void InputDataCoordTransformPoint(int indx, double *x1, double *x2, double *x3){
/*!
 * Converts spatial points from one geometry to another using macros
 *
 * \param [in] indx   input data array element (handle)
 * \param [in] x1     pointer to the coordinate to be transformed
 * \param [in] x2     pointer to the coordinate to be transformed
 * \param [in] x3     pointer to the coordinate to be transformed
 *
 *********************************************************************** */
    inputData *id = id_stack + indx;

// TODO: Actually do this check
    if (id->geom == CARTESIAN) {
        double xc1, xc2, xc3;
        xc1 = CART1(*x1, *x2, *x3);
        xc2 = CART2(*x1, *x2, *x3);
        xc3 = CART3(*x1, *x2, *x3);
        *x1 = xc1;
        *x2 = xc2;
        *x3 = xc3;
    }

    else if (id->geom == CYLINDRICAL) {
        double R, z, phi;
        R   = CYL1(*x1, *x2, *x3);
        phi = CYL2(*x1, *x2, *x3);
        z   = *x3;
        *x1 = R;
        *x2 = z;
        *x3 = phi;
    }

    else if (id->geom == POLAR) {
        double R, phi, z;
        R   = POL1(*x1, *x2, *x3);
        phi = POL2(*x1, *x2, *x3);
        z   = POL3(*x1, *x2, *x3);
        *x1 = R;
        *x2 = phi;
        *x3 = z;
    }

    else if (id->geom == SPHERICAL) {
        double R, theta, phi;
        R     = SPH1(*x1, *x2, *x3);
        theta = SPH2(*x1, *x2, *x3);
        phi   = SPH3(*x1, *x2, *x3);
        *x1 = R;
        *x2 = theta;
        *x3 = phi;
    }
}

/* -- AYW */


/* AYW --
 * Geometric velocity conversions.
 * NOTE: Assumptions on alignment made.
*/

/* ********************************************************************* */
void InputDataCoordTransformVector(int indx, double *x1, double *x2, double *x3, double ***vx1, double ***vx2,
                                   double ***vx3) {
/*!
 * Converts vector components from one geometry to another using macros
 * Unlike coord_transform_point, this is done in one sweep on the domain
 *
 * \param [in] indx   input data array element (handle)
 * \param [in] x1     pointer to the array of coordinates to the vector (usually grid->x[IDIR])
 * \param [in] x2     pointer to the array of coordinates to the vector (usually grid->x[JDIR])
 * \param [in] x3     pointer to the array of coordinates to the vector (usually grid->x[KDIR])
 * \param [in] v1     pointer to the 3D array of vector x1 component
 * \param [in] v2     pointer to the 3D array of vector x2 component
 * \param [in] v3     pointer to the 3D array of vector x3 component
 *
 *********************************************************************** */

    size_t k, j, i;
    double v1, v2, v3;

    inputData *id = id_stack + indx;

    if (id->geom == CARTESIAN) {
        TOT_LOOP(k, j, i) {
                    v1 = vx1[k][j][i];
                    v2 = vx2[k][j][i];
                    v3 = vx3[k][j][i];
                    EXPAND(vx1[k][j][i] = VCART_1(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx2[k][j][i] = VCART_2(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx3[k][j][i] = VCART_3(x1[i], x2[j], x3[k], v1, v2, v3););
                }
    } else if (id->geom == CYLINDRICAL) {
        TOT_LOOP(k, j, i) {
                    v1 = vx1[k][j][i];
                    v2 = vx2[k][j][i];
                    v3 = vx3[k][j][i];
                    EXPAND(vx1[k][j][i] = VCYL_1(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx2[k][j][i] = VCYL_2(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx3[k][j][i] = VCYL_3(x1[i], x2[j], x3[k], v1, v2, v3););
                }
    } else if (id->geom == POLAR) {
        TOT_LOOP(k, j, i) {
                    v1 = vx1[k][j][i];
                    v2 = vx2[k][j][i];
                    v3 = vx3[k][j][i];
                    EXPAND(vx1[k][j][i] = VPOL_1(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx2[k][j][i] = VPOL_2(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx3[k][j][i] = VPOL_3(x1[i], x2[j], x3[k], v1, v2, v3););

                }
    } else if (id->geom == SPHERICAL) {
        TOT_LOOP(k, j, i) {
                    v1 = vx1[k][j][i];
                    v2 = vx2[k][j][i];
                    v3 = vx3[k][j][i];
                    EXPAND(vx1[k][j][i] = VSPH_1(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx2[k][j][i] = VSPH_2(x1[i], x2[j], x3[k], v1, v2, v3);,
                           vx3[k][j][i] = VSPH_3(x1[i], x2[j], x3[k], v1, v2, v3););
                }
    } else {
        print("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
        QUIT_PLUTO(1);
    }
}

/* -- AYW */
