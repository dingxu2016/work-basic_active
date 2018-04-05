#ifndef LIST_H
#define LIST_H

#include"sys.h"

#define nlcut 0.3
#define listnum 30

typedef struct{
    double x;
    double y;
} tpold;

extern tpold *old_pos;

extern int *countn;
extern int **nl;


void alloc_list(void);
void make_list(void);
int  check_list(void);
double max(double x, double y);

#endif
