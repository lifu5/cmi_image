#include <stdio.h>
#include "cmibase.h"


void myradon_t(struct cmi_image* img, struct cmi_image* sinog){
    int width = row(img);
    int height = col(img);

    int n_beam = row(sinog);
    int n_angle = col(sinog);
    //printf("%d\n", n_angle);
    for (int x = 0; x<width; x++)
        for (int y = 0; y<height; y++){
            
        }
}

int main(){
    printf("hello\n");
    int size = sizeof(void);
    printf("%d\n", size);

    int n_beam = 300;
    int n_angle = 240;

    int width = 400;
    int height = 400;

    struct cmi_image* img = allocimage2d("image" , width, height, JINT);
    struct cmi_image* sinogram = allocimage2d("sinogram", n_beam, n_angle, JFLOAT);

    myradon_t(img, sinogram);

    freeimage(img);
    freeimage(sinogram);

    return 0;
}
