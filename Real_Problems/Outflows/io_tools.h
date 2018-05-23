//
// Created by Alexander Y. Wagner on 2017/06/29.
//

#ifndef PLUTO_IO_TOOLS_H
#define PLUTO_IO_TOOLS_H

double OutputContextEnter(const char *fname, FILE **fp, double next_output, double output_rate);

double OutputContextExit(FILE **fp, double next_output, double output_rate);

#endif //PLUTO_IO_TOOLS_H
