# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "rkf45.h"


void mytest(){
# define NEQN 2

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Solve a vector equation using R8_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y'(1) =  Y(2)\n" );
  printf ( "  Y'(2) = -Y(1)\n" );
  printf ( "\n" );

  abserr = sqrt(r8_epsilon());
  relerr = sqrt(r8_epsilon());

  flag = 1;

  t_start = 0.0;
  t_stop = 2.0*3.14159265;

  n_step = 12;

  t = 0.0;
  t_out = 0.0;

  y[0] = 1.0;
  y[1] = 0.0;
  r8_f2 ( t, y,  yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y(1)       Y(2)\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1] );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start 
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start 
            + ( double ) (	   i_step ) * t_stop ) 
            / ( double ) ( n_step );

    flag = r8_rkf45 ( r8_f2, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1] );
  }

  return;
# undef NEQN
}

