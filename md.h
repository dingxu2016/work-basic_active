#ifndef  MD_H
#define  MD_H

//#define alpha 2.0
#define kspring 10.0
#define mob 1.0

void cal_force(void);
void move(double v0, double dt, double x1, double x2, int iseed);

#endif
