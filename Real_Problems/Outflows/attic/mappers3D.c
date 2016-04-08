/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief 3D wrapper for conservative/primitive conversion.

  Provide 3D wrappers to the standard 1D conversion functions
  ConsToPrim() and PrimToCons().

  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 10, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ConsToPrim3D (Data_Arr U, Data_Arr V, Grid *grid)
/*!
 *  Convert a 3D array of conservative variables \c U to
 *  an array of primitive variables \c V.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]   U      pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [out]  V      pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [in]   grid   pointer to an array of Grid structures
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   current_dir;
  static unsigned char *flag;
  static double **v;
  double prsfix,g,scrh,rhofix,engfix,m2;

  if (v == NULL){
    v    = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);   
  }

  current_dir = g_dir;  /* save current direction */
  g_dir = IDIR;
  KDOM_LOOP(k) { g_k = k;
  JDOM_LOOP(j) { g_j = j;
    ConsToPrim (U[k][j], v, IBEG, IEND, flag);
    IDOM_LOOP(i){ VAR_LOOP(nv){ 
    //-----DM 26feb,2015: fix negative pressure, density, energy----//
	if ((flag[i] & RHO_FAIL) == RHO_FAIL){
	              rhofix=(   V[RHO][k+1][j+1][i-1]+V[RHO][k+1][j+1][i]+V[RHO][k+1][j+1][i+1]
				+V[RHO][k+1][j][i-1]+V[RHO][k+1][j][i]+V[RHO][k+1][j][i-1]
				+V[RHO][k+1][j-1][i-1]+V[RHO][k+1][j-1][i]+V[RHO][k+1][j-1][i+1]
				+V[RHO][k][j+1][i-1]+V[RHO][k][j+1][i]+V[RHO][k][j+1][i+1]
				+V[RHO][k][j][i-1]+V[RHO][k][j][i-1]
				+V[RHO][k][j-1][i-1]+V[RHO][k][j-1][i]+V[RHO][k][j-1][i+1]
				+V[RHO][k-1][j+1][i-1]+V[RHO][k-1][j+1][i]+V[RHO][k-1][j+1][i+1]
				+V[RHO][k-1][j][i-1]+V[RHO][k-1][j][i]+V[RHO][k-1][j][i-1]
				+V[RHO][k-1][j+1][i-1]+V[RHO][k-1][j+1][i]+V[RHO][k-1][j+1][i+1]
				  )/26.0;  
		if (rhofix <= 0.0) { print1("RHO still neg [%d,%d,%d] \n",i,j,k); rhofix=g_smallDensity;}
		if (rhofix != rhofix) {print1("RHO is NAN [%d,%d,%d] \n",i,j,k); rhofix=g_smallDensity;}
		V[RHO][k][j][i]=rhofix;

	}// RHO_FAIL

	if ( (flag[i] & PRS_FAIL) == PRS_FAIL  ){
	            prsfix=(     V[PRS][k+1][j+1][i-1]+V[PRS][k+1][j+1][i]+V[PRS][k+1][j+1][i+1]
				+V[PRS][k+1][j][i-1]+V[PRS][k+1][j][i]+V[PRS][k+1][j][i-1]
				+V[PRS][k+1][j-1][i-1]+V[PRS][k+1][j-1][i]+V[PRS][k+1][j-1][i+1]
				+V[PRS][k][j+1][i-1]+V[PRS][k][j+1][i]+V[PRS][k][j+1][i+1]
				+V[PRS][k][j][i-1]+V[PRS][k][j][i-1]
				+V[PRS][k][j-1][i-1]+V[PRS][k][j-1][i]+V[PRS][k][j-1][i+1]
				+V[PRS][k-1][j+1][i-1]+V[PRS][k-1][j+1][i]+V[PRS][k-1][j+1][i+1]
				+V[PRS][k-1][j][i-1]+V[PRS][k-1][j][i]+V[PRS][k-1][j][i-1]
				+V[PRS][k-1][j+1][i-1]+V[PRS][k-1][j+1][i]+V[PRS][k-1][j+1][i+1]
				  )/26.0;
		if (prsfix <= 0.0) { print1("PRS still neg [%d,%d,%d] \n",i,j,k); prsfix=g_smallPressure;}
		if (prsfix != prsfix) {print1("PRS is NAN [%d,%d,%d] \n",i,j,k); prsfix=g_smallPressure;}
		V[PRS][k][j][i]=prsfix;



          	if ((flag[i] & ENG_FAIL) == ENG_FAIL){
	   		engfix=(         U[k+1][j+1][i-1][ENG]+U[k+1][j+1][i][ENG]+U[k+1][j+1][i+1][ENG]
					+U[k+1][j][i-1][ENG]+U[k+1][j][i][ENG]+U[k+1][j][i-1][ENG]
					+U[k+1][j-1][i-1][ENG]+U[k+1][j-1][i][ENG]+U[k+1][j-1][i+1][ENG]
					+U[k][j+1][i-1][ENG]+U[k][j+1][i][ENG]+U[k][j+1][i+1][ENG]
					+U[k][j][i-1][ENG]+U[k][j][i-1][ENG]
					+U[k][j-1][i-1][ENG]+U[k][j-1][i][ENG]+U[k][j-1][i+1][ENG]
					+U[k-1][j+1][i-1][ENG]+U[k-1][j+1][i][ENG]+U[k-1][j+1][i+1][ENG]
					+U[k-1][j][i-1][ENG]+U[k-1][j][i][ENG]+U[k-1][j][i-1][ENG]
					+U[k-1][j+1][i-1][ENG]+U[k-1][j+1][i][ENG]+U[k-1][j+1][i+1][ENG]
					  )/26.0;
			if (engfix <= 0.0) { 
				print1("ENG still neg [%d,%d,%d] \n",i,j,k); 
				m2=EXPAND(U[k][j][i][MX1]*U[k][j][i][MX1], + U[k][j][i][MX2]*U[k][j][i][MX2], + U[k][j][i][MX3]*U[k][j][i][MX3]);
				engfix=sqrt(1.e-8 + m2 + U[k][j][i][RHO]*U[k][j][i][RHO]);
				}
			if (engfix != engfix) {
				print1("ENG is NAN [%d,%d,%d] \n",i,j,k); 
				m2=EXPAND(U[k][j][i][MX1]*U[k][j][i][MX1], + U[k][j][i][MX2]*U[k][j][i][MX2], + U[k][j][i][MX3]*U[k][j][i][MX3]);
				engfix=sqrt(1.e-8 + m2 + U[k][j][i][RHO]*U[k][j][i][RHO]);
				}
			U[k][j][i][ENG]=engfix;
	 		}//ENG_fail


		scrh   = 1.0/(U[k][j][i][ENG] + prsfix); 
    		EXPAND( V[VX1][k][j][i] = U[k][j][i][MX1]*scrh; ,
       		    	V[VX2][k][j][i] = U[k][j][i][MX2]*scrh; ,
         	  	V[VX3][k][j][i] = U[k][j][i][MX3]*scrh;)
   		g = EXPAND(V[VX1][k][j][i]*V[VX1][k][j][i], + V[VX2][k][j][i]*V[VX2][k][j][i], + V[VX3][k][j][i]*V[VX3][k][j][i]);
    		g = 1.0/sqrt(1.0 - g);
    		#if USE_FOUR_VELOCITY == YES
    		 EXPAND(V[VX1][k][j][i] *= g;  ,
       	 	   	V[VX2][k][j][i] *= g;  ,
          	  	V[VX3][k][j][i] *= g;)
    		#endif

	     /*	V[RHO][k][j][i]=(   V[RHO][k+1][j+1][i-1]+V[RHO][k+1][j+1][i]+V[RHO][k+1][j+1][i+1]
					+V[RHO][k+1][j][i-1]+V[RHO][k+1][j][i]+V[RHO][k+1][j][i-1]
					+V[RHO][k+1][j-1][i-1]+V[RHO][k+1][j-1][i]+V[RHO][k+1][j-1][i+1]
					+V[RHO][k][j+1][i-1]+V[RHO][k][j+1][i]+V[RHO][k][j+1][i+1]
					+V[RHO][k][j][i-1]+V[RHO][k][j][i-1]
					+V[RHO][k][j-1][i-1]+V[RHO][k][j-1][i]+V[RHO][k][j-1][i+1]
					+V[RHO][k-1][j+1][i-1]+V[RHO][k-1][j+1][i]+V[RHO][k-1][j+1][i+1]
					+V[RHO][k-1][j][i-1]+V[RHO][k-1][j][i]+V[RHO][k-1][j][i-1]
					+V[RHO][k-1][j+1][i-1]+V[RHO][k-1][j+1][i]+V[RHO][k-1][j+1][i+1]
				  	)/26.0;  

	     	V[VX1][k][j][i]=(   V[VX1][k+1][j+1][i-1]+V[VX1][k+1][j+1][i]+V[VX1][k+1][j+1][i+1]
					+V[VX1][k+1][j][i-1]+V[VX1][k+1][j][i]+V[VX1][k+1][j][i-1]
					+V[VX1][k+1][j-1][i-1]+V[VX1][k+1][j-1][i]+V[VX1][k+1][j-1][i+1]
					+V[VX1][k][j+1][i-1]+V[VX1][k][j+1][i]+V[VX1][k][j+1][i+1]
					+V[VX1][k][j][i-1]+V[VX1][k][j][i-1]
					+V[VX1][k][j-1][i-1]+V[VX1][k][j-1][i]+V[VX1][k][j-1][i+1]
					+V[VX1][k-1][j+1][i-1]+V[VX1][k-1][j+1][i]+V[VX1][k-1][j+1][i+1]
					+V[VX1][k-1][j][i-1]+V[VX1][k-1][j][i]+V[VX1][k-1][j][i-1]
					+V[VX1][k-1][j+1][i-1]+V[VX1][k-1][j+1][i]+V[VX1][k-1][j+1][i+1]
				  	)/26.0;  	

	     	V[VX2][k][j][i]=(   V[VX2][k+1][j+1][i-1]+V[VX2][k+1][j+1][i]+V[VX2][k+1][j+1][i+1]
					+V[VX2][k+1][j][i-1]+V[VX2][k+1][j][i]+V[VX2][k+1][j][i-1]
					+V[VX2][k+1][j-1][i-1]+V[VX2][k+1][j-1][i]+V[VX2][k+1][j-1][i+1]
					+V[VX2][k][j+1][i-1]+V[VX2][k][j+1][i]+V[VX2][k][j+1][i+1]
					+V[VX2][k][j][i-1]+V[VX2][k][j][i-1]
					+V[VX2][k][j-1][i-1]+V[VX2][k][j-1][i]+V[VX2][k][j-1][i+1]
					+V[VX2][k-1][j+1][i-1]+V[VX2][k-1][j+1][i]+V[VX2][k-1][j+1][i+1]
					+V[VX2][k-1][j][i-1]+V[VX2][k-1][j][i]+V[VX2][k-1][j][i-1]
					+V[VX2][k-1][j+1][i-1]+V[VX2][k-1][j+1][i]+V[VX2][k-1][j+1][i+1]
				  	)/26.0;  				
 
	     	V[VX3][k][j][i]=(   V[VX3][k+1][j+1][i-1]+V[VX3][k+1][j+1][i]+V[VX3][k+1][j+1][i+1]
					+V[VX3][k+1][j][i-1]+V[VX3][k+1][j][i]+V[VX3][k+1][j][i-1]
					+V[VX3][k+1][j-1][i-1]+V[VX3][k+1][j-1][i]+V[VX3][k+1][j-1][i+1]
					+V[VX3][k][j+1][i-1]+V[VX3][k][j+1][i]+V[VX3][k][j+1][i+1]
					+V[VX3][k][j][i-1]+V[VX3][k][j][i-1]
					+V[VX3][k][j-1][i-1]+V[VX3][k][j-1][i]+V[VX3][k][j-1][i+1]
					+V[VX3][k-1][j+1][i-1]+V[VX3][k-1][j+1][i]+V[VX3][k-1][j+1][i+1]
					+V[VX3][k-1][j][i-1]+V[VX3][k-1][j][i]+V[VX3][k-1][j][i-1]
					+V[VX3][k-1][j+1][i-1]+V[VX3][k-1][j+1][i]+V[VX3][k-1][j+1][i+1]
				  	)/26.0;  	
	*/

	
	 		
	}// PRS_FAIL

	V[nv][k][j][i] = v[i][nv];
	}};
  }}
  g_dir = current_dir;  /* restore current direction */
}
/* ********************************************************************* */
void PrimToCons3D (Data_Arr V, Data_Arr U, Grid *grid)
/*!
 *  Convert a 3D array of primitive variables \c V  to
 *  an array of conservative variables \c U.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]    V     pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [out]   U     pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [in]   grid   pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   current_dir;
  static double **v;

  if (v == NULL) v = ARRAY_2D(NMAX_POINT, NVAR, double);

/* ------------------------------------------------------------
     Convert solution vector from primitive to conservative
   ------------------------------------------------------------ */

  current_dir = g_dir; /* save current direction */
  g_dir = IDIR;
  KDOM_LOOP(k) { g_k = k;
  JDOM_LOOP(j) { g_j = j;
    IDOM_LOOP(i)  VAR_LOOP(nv) v[i][nv] = V[nv][k][j][i];
    PrimToCons (v, U[k][j], IBEG, IEND);
  }}
  g_dir = current_dir; /* restore current direction */
}
