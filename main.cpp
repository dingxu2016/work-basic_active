#include<cstdio>
#include<cstdlib>
#include<cmath>
#include"sys.h"
#include"config.h"
#include"list.h"
#include"md.h"
#include"mathe.h"

int main(void){
    double dt, v0;
    int iseed;
    long int step;
    char str[80], str1[80], str2[80], str3[80], s[80];
    char str_1[80], str_2[80], str_3[80], str_4[80], str_5[80];
    char str_6[80], str_7[80], str_8[80], str_9[80], str_10[80];

    scanf("%le", &dt);
    scanf("%le", &v0);
    scanf("%d", &iseed);
    scanf("%d", &sys.natom);
    //scanf("%le", &sys.phi);
    scanf("%le", &box.x);
    scanf("%ld", &step);

    box.y    = box.x ;
    box.xinv = 1.0 / box.x ;
    box.yinv = 1.0 / box.y ;

    check_cal_cluster = false;
    sigma = sqrt( 2.0 * Dr );

    alloc_atom();
    gen_rand_con(iseed);
    sprintf(str, "initiation_%.2f_%d_%d.txt", box.x, sys.natom, iseed);
    FILE *fp = fopen(str, "w+");
    for(int i=0; i<sys.natom; i++){
        fprintf(fp, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[i].x, atom[i].y, atom[i].r, atom[i].vx, atom[i].vy);
    }
    fclose(fp);

    int tj = 0; //count the  times about call meke_list


    sprintf(str_1, "particle_%.2f_%d_%d_data.txt", box.x, 1, iseed);
    sprintf(str_2, "particle_%.2f_%d_%d_data.txt", box.x, 2, iseed);
    sprintf(str_3, "particle_%.2f_%d_%d_data.txt", box.x, 3, iseed);
    sprintf(str_4, "particle_%.2f_%d_%d_data.txt", box.x, 4, iseed);
    sprintf(str_5, "particle_%.2f_%d_%d_data.txt", box.x, 5, iseed);
    sprintf(str_6, "particle_%.2f_%d_%d_data.txt", box.x, 6, iseed);
    sprintf(str_7, "particle_%.2f_%d_%d_data.txt", box.x, 7, iseed);
    sprintf(str_8, "particle_%.2f_%d_%d_data.txt", box.x, 8, iseed);
    sprintf(str_9, "particle_%.2f_%d_%d_data.txt", box.x, 9, iseed);
    sprintf(str_10, "particle_%.2f_%d_%d_data.txt", box.x, 10, iseed);
    FILE *fp_1 = fopen(str_1, "w+");
    FILE *fp_2 = fopen(str_2, "w+");
    FILE *fp_3 = fopen(str_3, "w+");
    FILE *fp_4 = fopen(str_4, "w+");
    FILE *fp_5 = fopen(str_5, "w+");
    FILE *fp_6 = fopen(str_6, "w+");
    FILE *fp_7 = fopen(str_7, "w+");
    FILE *fp_8 = fopen(str_8, "w+");
    FILE *fp_9 = fopen(str_9, "w+");
    FILE *fp_10 = fopen(str_10, "w+");

    alloc_list();
    make_list();
    cal_force();
    for (int i = 0; i < step; i++){
        move(v0, dt, mean, sigma, iseed);
        if(check_list()){
            make_list();
            tj += 1;
            printf("tj=%d\ti=%d\n", tj, i);
        }
        cal_force();

        if( i > 0 ){
            if( (i - 0) % 1 == 0){
                //printf("i=%d\n", i);
                //cal_cluster_size();
                //sprintf(str2, "cluster_size_%.2f_%d_%d_%d_data.txt", sys.phi, sys.natom, iseed, i);
                //FILE *fp2 = fopen(str2, "w+");
                for (int j = 0; j < sys.natom; j++){
                    atom[j].x -= round(atom[j].x * box.xinv) * box.x;
                    atom[j].y -= round(atom[j].y * box.yinv) * box.y;
                    if ( j == 0 )
                        fprintf(fp_1, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 1 )
                        fprintf(fp_2, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 2 )
                        fprintf(fp_3, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 3 )
                        fprintf(fp_4, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 4 )
                        fprintf(fp_5, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 5 )
                        fprintf(fp_6, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 6 )
                        fprintf(fp_7, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 7 )
                        fprintf(fp_8, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 8 )
                        fprintf(fp_9, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    if ( j == 9 )
                        fprintf(fp_10, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    //fprintf(fp2, "%d\t%d\n", j, cluster[j]);
                }
                //fclose(fp2);
            }   
        }
    }
    fclose(fp_1);
    fclose(fp_2);
    fclose(fp_3);
    fclose(fp_4);
    
    sprintf(str3, "final_%.2f_%d_%d.txt", box.x, sys.natom, iseed);
    FILE *fp3 = fopen(str3, "w+");
    for(int i=0; i<sys.natom; i++){
        atom[i].x -= round(atom[i].x * box.xinv) * box.x;
        atom[i].y -= round(atom[i].y * box.yinv) * box.y;
        fprintf(fp3, "%26.16e\t%26.16e\t%26.16e\n", atom[i].x, atom[i].y, atom[i].r); 
    }
    fclose(fp3);
    //printf("L=%26.16e", box.x);

    free(atom);
    free(old_pos);
    free(countn);
    free(nl);
    if ( check_cal_cluster == true)
        free(cluster);

    return 0;
}
