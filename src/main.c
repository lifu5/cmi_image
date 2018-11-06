#include <stdio.h>
#include <math.h>
#include "cmibase.h"

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

int main(){
    printf("hello\n");
    //int size = sizeof(void);
    //printf("%d\n", size);
    
    //set the parameters
    int n_beam = 300;
    int n_angle = 200;

    int width = 400;
    int height = 400;
    
    // create two new images
    struct cmi_image* img = allocimage2d("image" , width, height, JINT);
    struct cmi_image* sinogram = allocimage2d("sinogram", n_beam, n_angle, JFLOAT);
    
    //set response point
    point2d pos = {200, 200};
    int* data = img->data;
    data[(int)(pos.xx*width + pos.yy)] = 1;
    
    myradon_t(img, sinogram);
    writerawimage(img, "1i.DAT");
    writerawimage(sinogram, "1r.DAT");
    printf("-----\n");
    freeimage(img);
    freeimage(sinogram);
    return 0;
}
