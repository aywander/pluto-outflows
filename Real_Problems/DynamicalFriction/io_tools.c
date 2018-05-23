//
// Created by Alexander Y. Wagner on 2017/06/29.
//

#include "pluto.h"
#include "io_tools.h"

/*********************************************************************** */
double OutputContextEnter(const char *fname, FILE **fp, double next_output, double output_rate) {
/*!
 * Function called before specifying output in <...>Output fucntions.
 *
 *********************************************************************** */
    /* Open file if first timestep (but not restart).
     * We always write out the first timestep. */
    if (g_stepNumber == 0) {
        *fp = fopen(fname, "w");
        next_output += output_rate + 1;
    }

    /* Prepare for restart or append if this is not step 0  */
    else {

        /* In case of restart, get last timestamp
         * and determine next timestamp */
        if (next_output < 0.0) {
            char sline[512];
            *fp = fopen(fname, "r");
            while (fgets(sline, 512, *fp)) { }
            sscanf(sline, "%lf\n", next_output);
            next_output += output_rate + 1;
            fclose(*fp);
        }

        /* Append if next output step has been reached */
        if (g_time > next_output) *fp = fopen(fname, "a");

    }

    return next_output;
}

/*********************************************************************** */
double OutputContextExit(FILE **fp, double next_output, double output_rate) {
/*!
 * Function called after specifying output in <...>Output fucntions.
 *
 *********************************************************************** */

    if (g_time > next_output) {

        next_output += output_rate;

        fclose(*fp);

    }

    return next_output;
}