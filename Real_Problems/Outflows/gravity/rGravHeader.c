#include "pluto.h"
#include "definitions_usr.h"
#include "rGravTable.h"
#include "ccl.h"

#ifdef GRAV_TABLE

double potential[NPOTVAR];

int rGravHeader(){
  /*
   * This routine reads the header file to the gravity file.
   * The name of the file is GRAVHNAME and is set in definitions_usr.h
   *
   * The structure of the header file is that of a configuration file. 
   * We use the ccl library to parse it. 
   *
   * Note that the ccl library sometimes has trouble compiling. Set
   * CCL_OK in definitions_usr.h accordingly to skip reading the header.
   *
   * The values should be in cgs units.
   */

  struct ccl_t config;
  const struct ccl_pair_t *iter;
  char * pEnd;

  /* Set configuration file details */
  config.comment_char = '#';
  config.sep_char = '=';
  config.str_char = '"';

  /* Parse the file */
  ccl_parse(&config, GRAVHNAME);

  /* Get all key value pairs */

#if GRAV_POTENTIAL == HQNFW
  potential[AHQ]  = strtod ( ccl_get ( &config, "a_hq"     ) , &pEnd ) ;
  potential[DHQ]  = strtod ( ccl_get ( &config, "rho0_hq"  ) , &pEnd ) ;
  potential[ANFW] = strtod ( ccl_get ( &config, "a_nfw"    ) , &pEnd ) ;
  potential[DNFW] = strtod ( ccl_get ( &config, "rho0_nfw" ) , &pEnd ) ;
  potential[NHOT] = strtod ( ccl_get ( &config, "nhot"     ) , &pEnd ) ;
  potential[THOT] = strtod ( ccl_get ( &config, "thot"     ) , &pEnd ) ;

#elif GRAV_POTENTIAL == DOUBLE_ISO
  potential[KAP]  = strtod ( ccl_get ( &config, "kappa"   ) , &pEnd ) ;
  potential[LAM]  = strtod ( ccl_get ( &config, "lambda"  ) , &pEnd ) ;
  potential[RDM]  = strtod ( ccl_get ( &config, "rd"      ) , &pEnd ) ;
  potential[DHOT] = strtod ( ccl_get ( &config, "rho_hot" ) , &pEnd ) ;
  potential[THOT] = strtod ( ccl_get ( &config, "thot"    ) , &pEnd ) ;

#endif

  /* Clean up */
  ccl_release(&config);

  return 0;
}


#endif

