#include <stdio.h>
#include <math.h>
#include "cmibase.h"

float img2phy_coorx(int x, int width){
    float dx = 1./width;
    return (x-0.5)*dx - 0.5;
    
}

float img2phy_coory(int y, int height){
    float dy = 1./height;
    return (y-0.5)*dy - 0.5;
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
    
    int* ITdata = mask->data;
    for (int x = 0; x < width; x++){
        for (int y = 0; y<height; y++){
            double temp =
            cmi_sqr((x - (width+1)/2.)/width)+cmi_sqr((y - (height+1)/2.)/height);
            if (temp <= 1) ITdata[x*height + y] = 1;
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
    //printf("%f\n",d_m);
    cmilist* distlist = (cmilist*)malloc(sizeof(cmilist)*M);
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
                    double dist2 = xx*sin(theta) + yy*cos(theta) + 1;
                    int index_m = (int)(dist1/d_m) + n_bin2/2;//FIXME: only for even n_bin2
                    point2d point = {x,y};
                    mynode* newnode = initnode(dist2, point);
                    //printf("%d,%d,%d\n", j*n_bin2+index_m,M, index_m);
                    listinsert(&distlist[j*n_bin2+index_m], newnode);
                }

            }
        }
    }
    for (int i =0; i<M; i++)
        listsort(&distlist[i]);
}

void compute_systemM(
    struct cmi_image* img,
    struct cmi_image* mask,
    struct cmi_image* sinog){
    
    printf("run compute system func\n");
    int width = row(img);
    int height = col(img);

    //int n_bin1 = row(sinog);
    int n_angle = col(sinog);
    int n_bin2  = 100;
    int N = width*height;
    int M = n_bin2*n_angle;

    struct cmi_image* matrix_A = allocimage2d("sysA", M ,N, JFLOAT);
    
    get_matrix_A(mask, matrix_A, n_bin2, n_angle);
    printf("get matrix A\n");
    freeimage(matrix_A);   
}


int main(){
    printf("hello\n");
    //int size = sizeof(void);
    //printf("%d\n", size);
    
    //set the parameters
    int n_bin1 = 100;
    int n_angle = 80;

    int width = 256;
    int height = 256;
    
    // create two new images
    struct cmi_image* img = allocimage2d("image" , width, height, JINT);
    struct cmi_image* sinogram = allocimage2d("sinogram", n_bin1, n_angle, JFLOAT);
    struct cmi_image* mask = allocimage2d("mask", width, height, JINT);
    get_mask(mask);
    //set response point
    point2d pos = {300, 200};
    int* data = img->data;
    data[(int)(pos.xx*width + pos.yy)] = 1;
    data[(int)(pos.xx*width +1 + pos.yy)] = 2;
    data = NULL;// discard this pointer
    compute_systemM(img, mask, sinogram);   
    
    //myradon_t(img, sinogram);
    //writerawimage(img, "1i.DAT");
    //writerawimage(sinogram, "1r.DAT");

    //resetimage(sinogram);
    
    printf("-----\n");
    freeimage(img);
    freeimage(sinogram);
    return 0;
}
