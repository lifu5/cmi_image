#include <stdio.h>
#include <math.h>
#include "cmibase.h"

float img2phy_coorx(int x, int width){
    float dx = 1./width;
    return (x+0.5)*dx - 0.5;    
}

float img2phy_coory(int y, int height){
    float dy = 1./height;
    return (y+0.5)*dy - 0.5;
}

void myradon_t(struct cmi_image* img, struct cmi_image* sinog){
    printf("call randon func\n");
    int width = row(img);
    int height = col(img);

    int n_beam = row(sinog);
    int n_angle = col(sinog);

    float *FLdata = sinog->data;
    int* ITdata = img->data;
    double delta_b = 1.414*cmi_max(width, height)/n_beam; //only for square..

    for (int j = 0; j < n_angle; j++){  //loop from [0, pi]
        //printf("%d\n", i);
        double theta = M_PI*j/(double)n_angle;
        for (int x = 0; x<width; x++){
            for (int y = 0; y<height; y++){
                //int a = 1;
                double temp = (x - width/2.)*cos(theta) + (y - height/2.)*sin(theta);
                //x* cos(a) + y*sin(a)
                int b = (int)temp/delta_b + n_beam/2; //get the correpsond b
                FLdata[b*n_angle + j] += ITdata[x*height + y];
                }
            }
        }
}

// use delta 
void get_mask(struct cmi_image* mask){
    int width = row(mask);
    int height = col(mask);
    
    float* FTdata = mask->data;
    for (int x = 0; x < width; x++){
        for (int y = 0; y<height; y++){
            double temp =
            cmi_sqr((x - (width+1)/2.)/width)+cmi_sqr((y - (height+1)/2.)/height);
            if (temp <= 1) FTdata[x*height + y] = 200;
        }
    }
}

//FOV:x^2+y^2<=1
//using delta basis function to support continuous space
//set u_tot as a constant at frist

double calcu_weight_parallel(
        grid_info* grid_points,
        int idx[], int mm,
        double theta,
        double d_bin,
        double pixel_area){

    int m[4];
    int ct_less=0,ct_equ=0,ct_more=0;
    int equ[3] = {-1,-1,-1};// record equ index
    int nequ[3] = {-1, -1,-1}; //
    for(int i=0; i<4; i++)
        m[i] = grid_points[idx[i]].index_m;
    if (m[0]==m[1] &&m[1]==m[2] &&m[2]==m[3] &&m[3]==mm) return pixel_area;
    for (int i=0; i<4; i++){
        if (m[i]<mm){
            nequ[ct_less+ct_more] = idx[i];
            ct_less++;
        }
        else if (m[i]==mm) equ[ct_equ++] =idx[i];
        else{
            nequ[ct_less + ct_more] = idx[i];
            ct_more++;
        }
    }
    double area = 0;
    if (1==ct_equ){   
        double dist_m = grid_points[equ[0]].dist1;
        //printf("%lf,%lf,%lf", dist_m, mm)
        if (ct_less==0){//boundary mm and mm+1
            double bdline = (mm+1)*d_bin-0.5;
            //if (sin(theta)==0 || cos(theta)==0) continue;
            area = cmi_sqr(dist_m - bdline)/sin(theta)/cos(theta)/2;
            //if (area <0.000001) printf("xxx102 %f,%f\n",theta*120/M_PI, dist_m - bdline);
            return cmi_abs(area);
        }
        else if (ct_more==0){//boundary mm-1 and mm
            double bdline = mm*d_bin-0.5;
            //if (sin(theta)==0 || cos(theta)==0) continue;
             area = cmi_sqr(dist_m - bdline)/sin(theta)/cos(theta)/2;
             //printf("%f,%f\n",area,pixel_area);
            //if (area <0.000001) printf("xxx101 %f,%f\n",theta*120/M_PI, dist_m - bdline);
             return cmi_abs(area);
        }else{
            area = pixel_area - calcu_weight_parallel(grid_points, idx, mm-1,theta,d_bin, pixel_area)-\
                   calcu_weight_parallel(grid_points, idx,mm+1,theta,d_bin, pixel_area);
            //FIXME: repeated computational cost
            return cmi_abs(area);
        }
    }
    else if (2 == ct_equ){   
        double dist_m1 = grid_points[equ[0]].dist1;
        double dist_m2 = grid_points[equ[1]].dist1;
        if(ct_less ==0){// boundary line: mm and m+1
            double bdline = (mm+1)*d_bin-0.5;
            if (theta == 0 || cmi_abs(theta-M_PI/2.)<0.02){
                area = cmi_abs(dist_m1 - bdline)*sqrt(pixel_area);
                // printf("%f,%f\n", dist_m1-bdline,dist_m2-bdline);
                // if (area <0.01 && area>0.01) printf("xxxta%f,%f, %f",theta, dist_m2 - bdline, dist_m1 - bdline);
            }//special case
            else{
                area = cmi_sqr(dist_m1 - bdline) - cmi_sqr(dist_m2 -bdline);
                area = cmi_abs(area)/sin(theta)/cos(theta)/2;
                //if (cmi_abs(area) <0.01) printf("x2xx%f,%f, %f\n",theta, dist_m2 - bdline, dist_m1 - bdline);
            }
        }
        else if (ct_more == 0){ //boundary line: mm-1 and mm
            double bdline  = mm*d_bin - 0.5;
            if (theta == 0 || cmi_abs(theta-M_PI/2.)<0.02){
                area = cmi_abs(dist_m1 - bdline)*sqrt(pixel_area);
            }//special case
            else{
                area = cmi_sqr(dist_m1 - bdline) - cmi_sqr(dist_m2 -bdline);
                area = cmi_abs(area)/sin(theta)/cos(theta)/2;
            }
        }else{
            area = pixel_area - calcu_weight_parallel(grid_points, idx, mm-1,theta,d_bin, pixel_area)-\
                   calcu_weight_parallel(grid_points, idx,mm+1,theta,d_bin, pixel_area);
            //FIXME: repeated computational cost
        }
    }
    else if (3==ct_equ){
        double dist_m1 = grid_points[nequ[0]].dist1;
        if (ct_less == 1){ //boundary: m-1,m
            double bdline = (mm)*d_bin - 0.5;
            area = pixel_area - cmi_abs( cmi_sqr(bdline - dist_m1)/sin(theta)/cos(theta)/2 );
            //area = pixel_area - calcu_weight_parallel(grid_points, idx, mm+1, theta, d_bin, pixel_area);
        }
        else if (ct_more == 1){
            double bdline = (mm+1)*d_bin - 0.5;
            //area = pixel_area - calcu_weight_parallel(grid_points, idx, mm-1, theta, d_bin, pixel_area); 
            area = pixel_area - cmi_abs( cmi_sqr(bdline - dist_m1)/sin(theta)/cos(theta)/2 );
            //area = pixel_area - calcu_weight_parallel(grid_points, idx, mm+1, theta, d_bin, pixel_area);
        }   
    }
    else{
        fprintf(stderr,"exceptional case");
        area = 1/2*pixel_area;//XXX bug
    }
    return cmi_abs(area);
}

void get_matrix_A(
    struct cmi_image* atten,
    struct cmi_image* matrix_A,
    grid_info* grid_points,
    int n_bin2,
    int n_angle,
    cd_nodelist* pdistlist){
    
    int width = row(atten);
    int height = col(atten);
    double d_m = 1./n_bin2;
    int N = width*height;
    int M =n_bin2*n_angle;
    double pixel_area = 1./N;
    int Ngrid = (width-1)*(height-1);
    
    printf("-------\n");
    
    //cd_nodelist* pdistlist = (cd_nodelist*)malloc(sizeof(cd_nodelist)*M);
    for (int i =0; i<M; i++) pdistlist[i].len = 0;
    for (int j = 0; j<n_angle; j++){ //loop from 0 to pi
        double theta = M_PI*j/(double)n_angle;
        for (int x =1; x<width-1; x++){ //don't consider the outer layer pixel
            double xx = img2phy_coorx(x, width);
            for (int y = 1; y<height-1; y++){
                double yy= img2phy_coory(y, height);
                double dist2 = -xx*sin(theta) + yy*cos(theta)+0.5;
                int idx[4]; // four corners 
                idx[0] = j*Ngrid + (x-1)*(height-1) +y-1;
                idx[1] = j*Ngrid + (x-1)*(height-1) +y;
                idx[2] = j*Ngrid + x*(height-1) +y-1;
                idx[3] = j*Ngrid + x*(height-1) +y;
                //printf("%d,%d\n", grid_points[idx1].len, idx1);
                if (grid_points[idx[0]].len!=0 &&
                    grid_points[idx[1]].len!=0 &&
                    grid_points[idx[2]].len!=0 &&
                    grid_points[idx[3]].len!=0){
                        int m1,m2,m3,m4;
                        m1 = grid_points[idx[0]].index_m;
                        m2 = grid_points[idx[1]].index_m;
                        m3 = grid_points[idx[2]].index_m;
                        m4 = grid_points[idx[3]].index_m;
                        int min_m = cmi_min(cmi_min(m1,m2), cmi_min(m3,m4));
                        int max_m = cmi_max(cmi_max(m1,m2), cmi_max(m3,m4));
                        //printf("%d\n",max_m-min_m+1);
                        //float temp = 0;
                        for (int mm = min_m; mm<=max_m; mm++){
                            double weight=\
                            calcu_weight_parallel(grid_points,\
                            idx, mm, theta, d_m, pixel_area);
                            //printf("(%f, %d)", weight, mm);
                            //temp += weight;
                            point2d point = {x,y};
                            //if ( j==0 ) printf("%d %d %f %f\n", x,y,yy,cos(theta));
                            cd_node* newnode = initnode(dist2 ,point, weight);
                            listinsert(&pdistlist[j*n_bin2+mm], newnode);
                        }
                        //printf("\n");
                        //if (max_m-min_m+1!=0) printf("%f,%f\n", temp, pixel_area);
                }
            }
        }
    }
    for (int i =0; i<M; i++){
        listsort(&pdistlist[i]);
    }

    float* FTdata = atten->data;
    float* FTdata_A = matrix_A->data;
    for (int i = 0; i<n_angle; i++){
        //double theta = M_PI*i/n_angle;
        for (int j = 0; j<n_bin2; j++){
            double _summ = 0;
            for (cd_node* p = pdistlist[i*n_bin2+j].head; p!=NULL; p = p->next){
                int xx = (int)p->coord.xx;
                int yy = (int)p->coord.yy;

                //printf("%f\n",p->weight);
                FTdata_A[(i*n_bin2+j)*N + (xx*height+yy)] = exp(-_summ);
                _summ += FTdata[xx*height+yy]*(p->weight);
                if (i==30){printf("(%d, %d, %f, %d, %f)", xx, yy, p->weight, j, p->dist); 
                printf("  %f\n", _summ);}
            }
        }
    }
    //printf("%f\n",d_m);
    /* obsolete code
     cd_nodelist* distlist = (cd_nodelist*)malloc(sizeof(cd_nodelist)*M);
    for (int i = 0; i < M; i++)
        distlist[i].len = 0;

    for (int j = 0; j < n_angle; j++){  //loop from [0, pi]
        double theta = M_PI*j/(double)n_angle;
        for (int x = 0; x<width; x++){ //loop 
            double xx = img2phy_coory(x, width);
            for (int y = 0; y<height; y++){
                double yy = img2phy_coory(y, height);
                double dist1 = xx*cos(theta) + yy*sin(theta);
                if (dist1>=-0.5 && dist1<=0.5){ // whether in FOV
                    double dist2 = xx*sin(theta) + yy*cos(theta) + 0.5;
                    int index_m = (int)((dist1+0.5)/d_m);//FIXME: only for even n_bin2
                    if (n_bin2 == index_m)  index_m = n_bin2-1;
                    //printf("%d\n", index_m);
                    point2d point = {x,y};
                    cd_node* newnode = initnode(dist2, point);
                    //printf("%d,%d,%d\n", j*n_bin2+index_m,M, index_m);
                    //printf("%f\n", newnode->dist);
                    listinsert(&distlist[j*n_bin2+index_m], newnode);
                }
            }
        }
    }
    
    for (int i =0; i<M; i++){
        listsort(&distlist[i]);
        //for (int j=0; j<n_angle; j++){
        //    if (i == j*n_bin2+ n_bin2/2) 
        //        printf("%d,%d\n", j,distlist[i].len);
       // }
    }
    //int i = 0;
    //for (mynode* p = distlist[150].head; p!=NULL; p = p->next){
    //    i++;
    //    printf("%d, %f, (%f,%f)\n", i,p->dist, p->coord.xx, p->coord.yy);
    float* FTdata = atten->data;
    float* FTdata_A = matrix_A->data;
    for (int i = 0; i<n_angle; i++){
        double theta = M_PI*i/n_angle;
        double ratio = (double)(distlist[n_bin2/2].len)/distlist[i*n_bin2+n_bin2/2].len;
        double ratio2 = cmi_max(cmi_abs(sin(theta)),cmi_abs(cos(theta)));
        //printf("%f\n",ratio);
        for (int j = 0; j<n_bin2; j++){
            double _summ = 0;
            for (cd_node* p = distlist[i*n_bin2+j].head; p!=NULL; p = p->next){
                int xx = (int)p->coord.xx;
                int yy = (int)p->coord.yy;
                FTdata_A[(i*n_bin2+j)*N + (xx*height+yy)] = exp(-_summ);
                _summ += FTdata[xx*height+yy]*ratio2;
        }
    }
    }*/
}


void get_system_H(
    struct cmi_image* matrix_A,
    cd_nodelist* pdistlist,
    int n_angle,
    int n_bin2,
    int width,
    int height,
    struct cmi_image* system_H){
    
    printf("run compute sinogram func\n");
    
    int N = width*height;
    //int M = n_bin2*n_angle;
    float* FTdataA = matrix_A->data;
    float* FTdataH = system_H->data;
    double d_m = 1./n_bin2;
    for (int j = 0; j<n_angle; j++){        
        double theta = M_PI*j/(double)n_angle;
        for (int x = 0; x<width; x++){ //loop 
            double xx = img2phy_coory(x, width);
            for (int y = 0; y<height; y++){
                double yy = img2phy_coory(y, height);
                double dist1 = xx*cos(theta) + yy*sin(theta);
                if (dist1>=-0.5 && dist1<=0.5){ // whether in FOV
                    int index_m = (int)((dist1+0.5)/d_m);//FIXME: only for even n_bin2
                    // search the weight in struct pdistlist
                    FTdataH[(j*n_bin2+index_m)*N + x*height +y] =
                    FTdataA[(j*n_bin2+index_m)*N + x*height +y];                   
                }
            }
        }
    }
}

void assign_grid(
        grid_info* grid_points,
        geom_paras paras,
        int Ngrid){
       
    int n_angle = paras.n_angle;
    int n_bin = paras.n_bin2;
    int width = paras.width;
    int height = paras.height;
    
    double d_bin = 1./n_bin;
    double Lb = paras.Lb;
    double interv = paras.interv;
    double offsetx = 0.5/width;
    double offsety = 0.5/height;

    for (int i = 0; i<n_angle; i++){
        double theta = M_PI*i/n_angle;
        for (int x = 0; x<width-1; x++){
            double xx = img2phy_coorx(x, width)+offsetx;
            for (int y = 0; y<height-1; y++){
                double yy = img2phy_coory(y, height) + offsety;
                if (cmi_sqr(xx) + cmi_sqr(yy)>0.25) continue; // in FOV
                double dist1 = xx*cos(theta) + yy*sin(theta);
                double dist2 = -xx*sin(theta) + yy*cos(theta) + 0.5;
                if (dist1>=-0.5 && dist1<=0.5){ 
                    int index_m = (int)((dist1+0.5)/d_bin);
                    int temp_idx = i*Ngrid + x*(height-1) +y;
                    grid_points[temp_idx].index_m = index_m;
                    grid_points[temp_idx].dist1 = dist1;
                    double tempd = d_bin*(dist2+interv+Lb)/Lb;
                    int nofcolli =(int)((tempd/d_bin) -1);
 //                   printf("%d\n", nofcolli);
                    for (int m = index_m - nofcolli;
                         m<=index_m + nofcolli; m++){
                        if (m>=0 && m<n_bin){
                           int_node*p = intnode(m);
                           grid_node_insert(
                            &grid_points[temp_idx],p);
                        }
                    }
                }
            }
        }
    }
    //for (int i = 0; i<width-1; i++)
    //    for (int j = 0; j<height-1; j++){
    //        printf("%d,%d,%d\n",grid_points[i*(height-1)+j].len,i,j);          
    //    }
}

void compute_system_mat(
    struct cmi_image* img,
    struct cmi_image* atten,
    struct cmi_image* sinog,
    geom_paras paras){
    
    printf("run compute system func\n");
    int width = row(img);
    int height = col(img);

    //int n_bin1 = row(sinog);
    int n_angle = col(sinog);
    int n_bin2  = row(sinog);
    int N = width*height;
    int M = n_bin2*n_angle;
    int Ngrid  =(width-1)*(height-1);
    
    // assign point on grid belong to which collimator m
    grid_info* grid_points = (grid_info*)malloc(sizeof(grid_info)*Ngrid*n_angle);
    for (int i = 0; i<Ngrid*n_angle; i++)
        grid_points[i].len = 0;
    assign_grid(grid_points, paras, Ngrid);
    // write and test
    /*
    char* fff = "index.DAT";
    FILE* fd = fopen(fff, "w");
    for (int i = 0; i<Ngrid*n_angle; i++){
        double temp = 0;
        if (grid_points[i].len!=0){
            fwrite(&grid_points[i].dist1, sizeof(double) ,1,fd);
        }
        else{
            fwrite(&temp, sizeof(double) ,1,fd);
        }
    }
    fclose(fd);
    */

    printf("creat image matirx\n");   
    struct cmi_image* matrix_A = allocimage2d("sysA", M ,N, JFLOAT);
    struct cmi_image* sys_H = allocimage2d("sysH", M, N , JFLOAT);
    
    cd_nodelist* pdistlist = (cd_nodelist*)malloc(sizeof(cd_nodelist)*M);
    get_matrix_A(atten, matrix_A, grid_points, n_bin2, n_angle, pdistlist);
    //update matrix A and pdistlist
    printf("get matrix A\n");
    get_system_H(matrix_A, pdistlist, n_angle, n_bin2, width, height, sys_H); 
    writerawimage(matrix_A, "H.DAT");
    freeimage(matrix_A);   
    freeimage(sys_H);
}

int main(){
    printf("hello\n");
    //int size = sizeof(void);
    //printf("%d\n", size);
    
    geom_paras paras;
    //set the geometry parameters
    
    paras.n_bin2 = 100;
    paras.n_angle = 120;
    
    paras.width = 128;
    paras.height = 128;

    paras.interv = 0.1;
    paras.Lb = 0.4;
    
    // create two new images
    struct cmi_image* img = allocimage2d("image" , paras.width, paras.height, JINT);
    struct cmi_image* sinogram = allocimage2d("sinogram", paras.n_bin2, paras.n_angle, JFLOAT);
    struct cmi_image* atten = allocimage2d("atten", paras.width, paras.height, JFLOAT);
    get_mask(atten);
    //set response point
    //point2d pos = {300, 200};
    //int* data = img->data;
    //data[(int)(pos.xx*width + pos.yy)] = 1;
    //data[(int)(pos.xx*width +1 + pos.yy)] = 2;
    //data = NULL;// discard this pointer
    compute_system_mat(img, atten, sinogram, paras);   
    
    //myradon_t(img, sinogram);
    //writerawimage(img, "1i.DAT");
    //writerawimage(sinogram, "1r.DAT");

    //resetimage(sinogram);
    
    printf("-----\n");
    freeimage(img);
    freeimage(sinogram);
    return 0;
}
