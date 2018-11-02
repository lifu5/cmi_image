#include<stdio.h>
#include<stdint.h>

#ifdef  _BASE_H
#define _BASE_H

#ifdef __cpluscplus
extern "C"{
#endif

# define INT32_MAX  2147483547
# define INT32_MIN (-INT32_MAX-1)
# define M_E    2.7182818284590452354   /* e */
# define M_PI   3.1415926535897932384  /* pi */

#define cmi_abs(X)  ((X)>=0?(X):-(X))
#define cmi_max(X,Y) ((X)>=(Y)?(X):(Y))
#define cmi_min(X,Y) ((X)<=(Y)?(X):(Y))
#define cmi_odd(X) ((X)&1)
#define cmi_even(X) (((X)&1)==0)
#define cmi_sqr(x) ((x)*(x))

struct cmi_image {

    char *name;
    int32_t n_row;
    int32_t n_col;
    int32_t n_depth;

    int32_t data_type;    //data storage type
    void* data;           //pointer for raw data
};

extern void writeimage(
    struct cmi_image* image,
    char* filename
);
// binary data

extern cmi_image* readimage(
    char* filename
);
// binary data

extern void freeimage(cmi_image* image);

extern cmi_image* zeros_image_2d(
    int32_t row,
    int32_t col,
    int32_t type
);

#ifdef __cplusplus
}
#endif
#endif
