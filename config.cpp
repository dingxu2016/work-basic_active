#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<random>
#include"sys.h"
#include"config.h"
#include"mathe.h"

tpatom *atom;

void alloc_atom(void){
    atom = (tpatom *)malloc(sys.natom * sizeof(tpatom));
}

void gen_rand_con(int iseed){
    //initiate atom's radius
    for(int i=0; i<sys.natom; i++){
        if(i < sys.natom / 2)
            atom[i].r = 0.5;
        else
            atom[i].r = 0.5 * ratio;
    }

    //calculate box's length
    double sdisk = 0.0;
    for(int i=0; i<sys.natom; i++)
        sdisk += PI * atom[i].r * atom[i].r;
    double volsquare = sdisk / sys.phi;
    double temp = sqrt(volsquare);
    box.x = temp;
    box.y = temp;
//    box.z = temp;
    box.xinv = 1.0 / temp;
    box.yinv = 1.0 / temp;
//    box.zinv = 1.0 / temp;

    //initiate randomly atom's position
    /******** C的随机数库写的，不过这个随机数库性能太差**********
    srand(iseed);           
    for(int i=0; i<sys.natom; i++){
        atom[i].x = ( (double)rand()/RAND_MAX - 0.5) * box.x;
        atom[i].y = ( (double)rand()/RAND_MAX - 0.5) * box.y;
    }
    */

    for (int i=0; i<sys.natom; i++)
    {
        atom[i].x = ( random_number(iseed) - 0.5 ) * box.x;
        atom[i].y = ( random_number(iseed) - 0.5 ) * box.y;
    }

    //initiate randomly atom's angle
    for(int i=0; i<sys.natom; i++){
        if(i<sys.natom/2)
            atom[i].theta = random_number(iseed) * PI;
        else
            atom[i].theta = atom[i-sys.natom/2].theta + PI; 
    }
}
