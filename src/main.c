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
            if (temp <= 1) FTdata[x*height + y] = 0.02;
        }
    }
}

//FOV:x^2+y^2<=1
//using delta basis function to support continuous space
//set u_tot as a constant at frist

void get_matrix_A(
    struct cmi_image* atten,
    struct cmi_image* matrix_A,
    int n_bin2,
    int n_angle){
    
    int width = row(atten);
    int height = col(atten);
    double d_m = 1./n_bin2;
    int M =n_bin2*n_angle;
    int N = width*height;
    //printf("%f\n",d_m);
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
    }
}


void get_system_H(
    struct cmi_image* matrix_A,
    struct cmi_image* system_H,
    int n_angle,
    int n_bin2,
    int width,
    int height){
    
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
                    FTdataH[(j*n_bin2+index_m)*N + x*height +y] =
                    FTdataA[(j*n_bin2+index_m)*N + x*height +y];                   
                }
            }
        }
    }
}

void assign_grid(
        int_nodelist* grid_assign,
        geom_paras paras,
        int Ngird){
    
    
    int n_angle = paras.n_angle;
    int n_bin = paras.n_bin2;
    int width = paras.width;
    int height = paras.height;
    
    double d_bin = 1/n_bin;
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
                double dist1 = xx*cos(theta) + yy*sin(theta);
                double dist2 = xx*sin(theta) + yy*cos(theta) + 0.5;
                if (dist1>=-0.5 && dist1<=0.5){ //in FOV
                    int index_m = (int)((dist1+0.5)/d_bin);
                    double tempd = n_bin*(dist2+interv+Lb)/Lb;
                    int nofcolli =(int)((tempd/d_bin) -1);
                    for (int m = index_m - nofcolli;
                         m<=index_m + nofcolli; m++){
                        if (m>=0 && m<n_bin){
                           int_node*p = intnode(m);
                           listinsert_int(grid_assign, p);
                        }
                    }

                }
                
            }
        }
    }
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
    int_nodelist* grid_assign = (int_nodelist*)malloc(sizeof(int_nodelist)*Ngrid*n_angle);
    for (int i = 0; i<Ngrid*n_angle; i++)
        grid_assign[i].len = 0;
    assign_grid(grid_assign, paras, Ngrid);
    printf("creat image matirx\n");   
    struct cmi_image* matrix_A = allocimage2d("sysA", M ,N, JFLOAT);
    struct cmi_image* sys_H = allocimage2d("sysH", M, N , JFLOAT);
    get_matrix_A(atten, matrix_A, n_bin2, n_angle);
    printf("get matrix A\n");
    //get_system_H(matrix_A, sys_H, n_angle, n_bin2, width, height); 
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
    
    paras.width = 144;
    paras.height = 144;

    paras.interv = 0.25;
    paras.Lb = 0.25;
    
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
