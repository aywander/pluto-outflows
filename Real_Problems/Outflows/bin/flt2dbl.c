#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include"SwapEndian.h"

#define NO 0
#define YES 1

void Print_Usage ();

/* ************************************************ */
int main(int argc, char **argv)
/* 
 *
 * PURPOSE
 *
 *   Convert a binary data file 
 *   from single precision to double
 *   precision.
 *
 *
 *
 *  SYNOPSIS :
 *
 *   flt2dbl file [options] 
 *
 *
 * LAST MODIFIED: 
 *
 *   2012-07-17 09:46 JST AYW Adapted from dbl2flt.c
 *
 *   3 June 2007 by A. Mignone, e-mail: mignone@to.astro.it 
 *
 *************************************************** */
{
  int   arg_is_file[65536], i, swap_flt, swap_dbl;
  int   single_output_file = NO;
  char  *fname_in, fname_out[256];
  double dbl_x;
  float  flt_x;
  FILE  *fin, *fout;

  printf ("flt2dbl - version 1.0\n");
  printf ("Last modified 2012-07-17\n");
  printf ("Original Copyright(C) 2006,2007 by A. Mignone (mignone@to.astro.it)\n\n");
  printf ("Adapted by AYW. \n\n");

  if (argc < 2) Print_Usage();

/* -----------------------------------------
           default values 
   ----------------------------------------- */

  swap_flt = swap_dbl = NO;

/* -------------------------------------------
      parse command line options
   ------------------------------------------- */

  for (i = 1; i < argc; i++){

    arg_is_file[i] = NO;
    if (!strcmp(argv[i], "--help") ||
        !strcmp(argv[i], "-help")){

      Print_Usage();

    }else if (!strcmp(argv[i],"-swap-flt")){

      swap_flt = YES;

    }else if (!strcmp(argv[i],"-swap-dbl")){

      swap_dbl = YES;

    }else if (!strcmp(argv[i],"-o")){
                                                                                                                               
     /* -- open single output file -- */

      sprintf (fname_out,"%s",argv[++i]);
      single_output_file = YES;
      fout = fopen(fname_out,"wb");

    }else{
      arg_is_file[i] = YES;
    }
  }

/* -------------------------------------------------
       open data file for reading and writing 
   ------------------------------------------------- */

  for (i = 1; i < argc; i++){

  /* -- loop only on arguments 
        representing a valid file name -- */ 

    if (!arg_is_file[i]) continue;
    if (!strcmp(argv[i],fname_out)) continue;
    
    fname_in = argv[i];
    if ((fin = fopen (fname_in, "r")) == NULL){
      printf ("! File %s does not exist\n",fname_in);
      break;
    }

  /* -- output name -- */

    if (single_output_file == NO){
      sprintf (fname_out,"%s.dbl",fname_in);
      fout = fopen(fname_out,"wb");
    }

    printf("Converting (float) file %s to (double) %s ...",
            fname_in,fname_out);
    fflush (stdout);

    for (;;){ 
      fread(&flt_x, sizeof(float), 1, fin);
      if (feof(fin)) break; /* check end of file before writing */

      if (swap_flt) SWAP_FLOAT  (flt_x);
      dbl_x = (double) flt_x;
      if (swap_dbl) SWAP_DOUBLE (dbl_x);
      fwrite (&dbl_x, sizeof(double), 1, fout);
    }

    printf("done.\n");
    fclose (fin);
    if (single_output_file == NO) fclose (fout);

  }

  if (single_output_file) fclose (fout);
  return(0);
}

/* ****************************************************************************  */
void Print_Usage ()
/*
 *
 *
 *
 *
 *
 ****************************************************************************** */
{
  printf ("Purpose: convert binary files from single to double precision.\n\n");
  printf ("Usage: flt2dbl file [options]\n\n");
  printf (" file             the name of a valid binary data file(s) in\n");
  printf ("                  double precision;\n");
  printf ("[options] are:\n\n");
  printf (" -o <name>        specify the name of the single precision output\n");
  printf ("                  file. The default is file.dbl;\n");
  printf (" -swap-dbl        swap endianity when writing double precision data.\n"); 
  printf ("                  Use this switch if the byte order of the data and\n");
  printf ("                  the processor differ.\n");
  printf (" -swap-flt        swap endianity when reading single precision data.\n");
  printf ("                  Use this switch if data and processor have the same \n");
  printf ("                  byte order, but you wish to write it on a machine with\n");
  printf ("                  different endianity.\n\n");

  printf ("Example:  Given the files (in single precision) rho.0001, rho.0002, ...,\n");
  printf ("-------   the command\n\n");
  printf ("               flt2dbl rho.*\n\n");
  printf ("          will convert all files beginning with 'rho.' into double precision\n");
  printf ("          files rho.0001.dbl rho.0002.dbl, ...\n\n");

  printf ("Example:  Given the files (in single precision) q.0001, q.0002, ...,\n");
  printf ("-------   the command\n\n");
  printf ("               flt2dbl q.* -o qall.bin\n\n");
  printf ("          will convert and append every file beginning with 'q.' into\n");
  printf ("          to the single file 'qall.bin'.\n\n");

  printf ("Example:       flt2dbl q.* -swap-flt -o qall.bin\n");
  printf ("-------   \n");
  printf ("          same as before, but swap the endianity when reading single\n");
  printf ("          precision data. This is useful, for example, when data is being\n");
  printf ("          read on a machine with byte order different from the processor\n");
  printf ("          where you intend to write the data.\n");

  exit(0);
}

