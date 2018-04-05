#include<cstdio>
#include<cstdlib>
#include<cmath>
#include"sys.h"
#include"config.h"
#include"list.h"

tpold *old_pos;
int   *countn;
int   **nl;

void alloc_list(void){
    countn  = (int *)  malloc((sys.natom - 1) * sizeof(int));
    old_pos = (tpold *)malloc(sys.natom * sizeof(tpold));
    nl      = (int **) malloc((sys.natom - 1) * sizeof(int *));
    for (int i = 0; i < sys.natom - 1; i++)
        nl[i] = (int *) malloc(listnum * sizeof(int));
}

void make_list(void){
    double xij, yij, rij, rlist;
    double rx, ry;

    rlist = 1.0 + nlcut;

    for (int i = 0; i<sys.natom - 1; i++)
        countn[i] = - 1;

    for (int i = 0; i < sys.natom - 1; i++){
        rx = atom[i].x;
        ry = atom[i].y;
        for (int j = i+1; j < sys.natom; j++){
            xij = (rx - atom[j].x) * box.xinv;
            xij = (xij - (double)round(xij)) * box.x;
            yij = (ry - atom[j].y) * box.yinv;
            yij = (yij - (double)round(yij)) * box.y;
            rij = sqrt(xij * xij + yij * yij);
            if( rij < rlist ){
                countn[i] += 1;
                nl[i][countn[i]] = j;
            }
            if(countn[i] >= listnum){
                printf("list max number problem");
                exit(0);
            }
        }
    }
    
    for(int i = 0; i < sys.natom; i++){
        old_pos[i].x = atom[i].x;
        old_pos[i].y = atom[i].y;
    }
}

int check_list(void){
    double maxdis = 0.0;

    for(int i=0; i<sys.natom; i++){
        maxdis = max(fabs(atom[i].x - old_pos[i].x), maxdis);
        maxdis = max(fabs(atom[i].y - old_pos[i].y), maxdis);
    }
    maxdis = 2.0 * sqrt(2.0 * maxdis * maxdis);

    if(maxdis > nlcut)
        return 1;
    else
        return 0;
}

double max(double x, double y){
    if(x > y)
        return x;
    else
        return y;
}
