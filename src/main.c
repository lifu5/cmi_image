#include <stdio.h>
#include "cmibase.h"

int main(){
    printf("hello\n");
    struct cmi_image* img;
    int size = sizeof(void);
    printf("%d\n", size);
    size = sizeof(int16_t);
    printf("%d\n", size);
    size = sizeof(double);
    printf("%d\n", size);
    return 0;
}
