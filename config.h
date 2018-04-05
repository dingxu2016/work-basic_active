/*******************************
 type definition of configuration
********************************/

#ifndef CONFIG_H
#define CONFIG_H

#include"sys.h"

typedef struct{
    double x, vx, fx;
    double y, vy, fy;
   // double z, vz, fz;
    double r;
    double theta;
}tpatom;

// particle's properties
extern tpatom *atom;

void alloc_atom(void);
void gen_rand_con(int iseed);

#endif
