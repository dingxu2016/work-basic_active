#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <cstdbool>
#include "mathe.h"
#include "sys.h"
#include "config.h"

double sigma;
double N[10];
int    *cluster;
bool   check_cal_cluster;

double random_number(int iseed)
{
    static std::default_random_engine random( iseed );
    static std::uniform_real_distribution<double> dis(0.0, 1.0);//半开区间，无1.0 
    return dis(random);
}

double Gaussian_noise(double x1, double x2, int iseed)
{
    static std::default_random_engine gs( iseed+1 );
    static std::normal_distribution<double> n(x1,x2);
    return n(gs);
}

/*
double Gaussian_noise(double x1, double x2){
    double s1, s2;
    double v1, v2;

    //s1 = (((double)rand()) + 1.0) / (((double)RAND_MAX) + 1.0);
    //s2 = (((double)rand()) + 1.0) / (((double)RAND_MAX) + 1.0);
    s1 = random_number();
    while(s1 == 0.0 ){
        s1 = random_number();
    }
    s2 = random_number();
    while(s2 == 0.0 ){
        s2 = random_number();
    }
    v1 = sqrt(- 2.0 * log(s1)) * cos(2.0 * PI * s2);
    v2 = x1 + v1 * x2;
    return v2;
}
*/

double array_mean(int *a, int len){    // len = length
    double sum = 0.0;
    double mean1 = 0.0;
    double N_sub;

    N_sub = (double) sys.natom / len / len;

    for(int i = 0; i < len; i++){
        for(int j = 0; j < len; j++){
            sum += (double)( *(a + i * len + j) - N_sub) * ( *(a + i * len  + j) - N_sub);
        }
    }
    
    mean1 = sqrt( (double) sum / len / len );

    return mean1;
}

void cal_std_deviation(void){
    int js1[2][2]      = {0};
    int js2[3][3]      = {0};
    int js3[4][4]      = {0};
    int js4[6][6]      = {0};
    int js5[10][10]    = {0};
    int js6[11][11]    = {0};
    int js7[16][16]    = {0};
    int js8[32][32]    = {0};
    int js9[45][45]    = {0};
    int js10[100][100] = {0};

    double len[10];
    len[0] = box.x / 2.0;
    len[1] = box.x / 3.0;
    len[2] = box.x / 4.0;
    len[3] = box.x / 6.0;
    len[4] = box.x / 10.0;
    len[5] = box.x / 11.0;
    len[6] = box.x / 16.0;
    len[7] = box.x / 32.0;
    len[8] = box.x / 45.0;
    len[9] = box.x / 100.0;

    int xx[10], yy[10];

    for (int i = 0; i < sys.natom; i++) {
        atom[i].x -= round(atom[i].x * box.xinv) * box.x;
        atom[i].y -= round(atom[i].y * box.yinv) * box.y;
        atom[i].x += box.x / 2.0;
        atom[i].y += box.y / 2.0;
        for (int j = 0; j < 10; j++){
            xx[j] = floor( atom[i].x / len[j] );
            yy[j] = floor( atom[i].y / len[j] );
            switch (j+1){
                case 1 : js1[xx[j]][yy[j]] += 1; break;
                case 2 : js2[xx[j]][yy[j]] += 1; break;
                case 3 : js3[xx[j]][yy[j]] += 1; break;
                case 4 : js4[xx[j]][yy[j]] += 1; break;
                case 5 : js5[xx[j]][yy[j]] += 1; break;
                case 6 : js6[xx[j]][yy[j]] += 1; break;
                case 7 : js7[xx[j]][yy[j]] += 1; break;
                case 8 : js8[xx[j]][yy[j]] += 1; break;
                case 9 : js9[xx[j]][yy[j]] += 1; break;
                case 10 : js10[xx[j]][yy[j]] += 1; break;
            }
        }
    }
    
    N[0] = array_mean(*js1, 2);
    N[1] = array_mean(*js2, 3);
    N[2] = array_mean(*js3, 4);
    N[3] = array_mean(*js4, 6);
    N[4] = array_mean(*js5, 10);
    N[5] = array_mean(*js6, 11);
    N[6] = array_mean(*js7, 16);
    N[7] = array_mean(*js8, 32);
    N[8] = array_mean(*js9, 45);
    N[9] = array_mean(*js10, 100);

}

void cal_cluster_size(void){
    double xij, yij, rij, dij, rx, ry, rr;
    int    **nl;
    bool   *visited;
    int    n = 0;

    nl = (int **) malloc(sizeof(int *) * sys.natom);
    for (int i = 0; i < sys.natom; i++){
        nl[i] = (int *) malloc(sizeof(int) * sys.natom);
    }
    visited = (bool *) malloc(sizeof(bool) * sys.natom);
    if (check_cal_cluster == false)
        cluster = (int *)  malloc(sizeof(int)  * sys.natom);

    check_cal_cluster = true;

    for (int i = 0; i < sys.natom; i++){
        visited[i] = false;
        cluster[i] =   0  ;
        for (int j = 0; j < sys.natom; j++)
            nl[i][j] = 0;
    }

    //make a list
    for (int i = 0; i < sys.natom - 1; i++){
        rx = atom[i].x;
        ry = atom[i].y;
        rr = atom[i].r;
        for (int j = i+1; j < sys.natom; j++){
            xij = (rx - atom[j].x) * box.xinv;
            xij = (xij - (double)round(xij)) * box.x;
            yij = (ry - atom[j].y) * box.yinv;
            yij = (yij - (double)round(yij)) * box.y;
            rij = sqrt(xij * xij + yij * yij);
            dij = rr + atom[j].r;
            if (rij <= dij ){
                nl[i][j] = 1;
                nl[j][i] = 1;
            }
        }
    }

    for(int i = 0; i < sys.natom; i++){
        if( visited[i] == false){
            visited[i] = true;
            cluster[n] =   1 ;
            cal_recursive(i, n, visited, nl);
            n++;
        }
    }
    free(visited);
    for (int i = 0; i < sys.natom; i++)
        free( *(nl+i) );
}

void cal_recursive(int p1, int n, bool *visited, int **nl){
    for (int i = 0; i < sys.natom; i++){
        int flag = 0;
        if( (nl[p1][i] == 1) && (visited[i] == false) ){
            cluster[n] += 1   ;
            visited[i] =  true;
            flag       =  1   ;
            if( flag  == 1)
                cal_recursive(i, n, visited, nl);
            else
                continue;
        }
    }
}
