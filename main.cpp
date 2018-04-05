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

    scanf("%le", &dt);
    scanf("%le", &v0);
    scanf("%d", &iseed);
    scanf("%d", &sys.natom);
    scanf("%le", &sys.phi);
    scanf("%ld", &step);

    check_cal_cluster = false;
    sigma = sqrt( 2.0 * Dr );

    alloc_atom();
    gen_rand_con(iseed);
    sprintf(str, "initiation_%.2f_%d_%d.txt", sys.phi, sys.natom, iseed);
    FILE *fp = fopen(str, "w+");
    for(int i=0; i<sys.natom; i++){
        fprintf(fp, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[i].x, atom[i].y, atom[i].r, atom[i].vx, atom[i].vy);
    }
    fclose(fp);

    int tj = 0; //count the  times about call meke_list

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

        if( i > 19999 ){
            if( (i - 20000) % 2000 == 0){
                printf("i=%d\n", i);
                cal_cluster_size();
                sprintf(str1, "snapshot_config_%.2f_%d_%d_%d_data.txt", sys.phi, sys.natom, iseed, i);
                sprintf(str2, "cluster_size_%.2f_%d_%d_%d_data.txt", sys.phi, sys.natom, iseed, i);
                FILE *fp1 = fopen(str1, "w+");
                FILE *fp2 = fopen(str2, "w+");
                for (int j = 0; j < sys.natom; j++){
                    atom[j].x -= round(atom[j].x * box.xinv) * box.x;
                    atom[j].y -= round(atom[j].y * box.yinv) * box.y;
                    fprintf(fp1, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    fprintf(fp2, "%d\t%d\n", j, cluster[j]);
                }
                fclose(fp1);
                fclose(fp2);
            }   
        }

    }
    
    sprintf(str3, "final_%.2f_%d_%d.txt", sys.phi, sys.natom, iseed);
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
