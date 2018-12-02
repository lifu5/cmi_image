#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmibase.h"
#include "cmiutils.h"

inline int row(struct cmi_image* f){
    return f->n_row;
}

inline int col(struct cmi_image* f){
    return f->n_col;
}

inline int depth(struct cmi_image* f){
    return f->n_depth;
}

int pixel_num(struct cmi_image* f){
    int32_t _row = row(f);
    int32_t _col = col(f);
    int32_t _depth = depth(f);
    int N = _row * _col * _depth;
    return N;
}

struct cmi_image* allocimage(
    char* name,
    int32_t rr,
    int32_t cc,
    int32_t dp,
    int32_t dt){
/* allocat space  return a struct pointer
 */
    int N = rr * cc * dp;
    struct cmi_image* output;
    int d_size;
    switch (dt){
        case JBOOL:     d_size = sizeof(_Bool); break;
        case JINT:      d_size = sizeof(int);   break;
        case JFLOAT:    d_size = sizeof(float); break;
        case JDOUBLE:   d_size = sizeof(double);break;
        default: fprintf(stderr,"invaild data_type %d\n", dt);
                 return NULL;
    } // select data type

    output = (struct cmi_image*) malloc(sizeof(struct cmi_image));
    if (output == NULL){
        fprintf(stderr,"malloc failed\n");
        return NULL;
    }

    output->data = (void *)calloc(1, N*d_size);
    if (output->data == NULL){
        fprintf(stderr, "calloc failed\n");
        return NULL;
    }
    if (name != NULL){
        output->name = (char *)calloc(1, strlen(name)+1);
        strcpy(output->name, name);
    }
    else
        output->name = NULL;
    //assign value    
    output->n_row = rr;
    output->n_col = cc;
    output->n_depth = dp;
    output->data_type = dt;

    return output;
}

struct cmi_image* allocimage2d(
    char* name,
    int32_t rr,
    int32_t cc,
    int32_t dt){
 
    int32_t depth = 1;
    return allocimage(name, rr, cc, depth, dt);
}

void freeimage(struct cmi_image* image){
    if (image->name != NULL) free(image->name);
    free(image->data);
    free(image);
}

void writerawimage(struct cmi_image* image, char* filename){
    FILE* fd = NULL;
    int32_t wrt;
    int _row = row(image);
    int _col = col(image);
    int _depth = depth(image);
    int N = _row * _col * _depth;

    fd = fopen(filename, "w");
    if (!fd){
        fprintf(stderr, "cannot open file: %s\n", filename);
        exit(0);
    }
    fwrite(&_row, sizeof(int), 1, fd);
    fwrite(&_col, sizeof(int), 1, fd);
    fwrite(&_depth, sizeof(int), 1, fd);
    if (image->data_type == JBOOL){
        wrt = fwrite((_Bool*)image->data, sizeof(_Bool), N, fd);
        if (wrt != N){
            fprintf(stderr, "write error!\n");
            exit(0);
        }
    }
    if (image->data_type == JINT){
        wrt = fwrite((int*)image->data, sizeof(int), N, fd);
        if (wrt != N){
            fprintf(stderr, "write error!\n");
            exit(0);
        }
    }
    if (image->data_type == JFLOAT){
        wrt = fwrite((float*)image->data, sizeof(float), N, fd);
        if (wrt != N){
            fprintf(stderr, "write error!\n");
            exit(0);
        }
    }
    if (image->data_type == JDOUBLE){
        wrt = fwrite((double*)image->data, sizeof(double), N, fd);
        if (wrt != N){
            fprintf(stderr, "write error!\n");
            exit(0);
        }
    }
    fclose(fd);
}

struct cmi_image* copyimage(struct cmi_image* src){
    int32_t row = src->n_row;
    int32_t col = src->n_col;
    int32_t dpt = src->n_depth;
    int32_t dt = src->data_type;

    struct cmi_image* dst;

    int N = row * col * dpt;
    dst = allocimage(NULL, row, col, dpt, dt);
    if (dst == NULL){
        fprintf(stderr, "allocate image failed\n");
        return NULL;
    }

    switch(dt){
        case JBOOL: memcpy(src->data, dst->data, (N*sizeof(_Bool))); break;
        case JINT: memcpy(src->data, dst->data, (N*sizeof(int32_t))); break;
        case JFLOAT: memcpy(src->data, dst->data, (N*sizeof(float))); break;
        case JDOUBLE: memcpy(src->data, dst->data, (N*sizeof(double))); break;
        default:
            fprintf(stderr, "bad data type\n");
            return NULL;
    }

    if (src->name != NULL){
        dst->name = (char*)calloc(1, strlen(src->name) +1);
        strcpy(dst->name, src->name);
    }
    return dst;
}

void doubleimage(struct cmi_image* g){
    int _row = row(g);
    int _col = col(g);
    int _depth = depth(g);

    int N = _row * _col * _depth;
    void* buffer = (void* )calloc(1, N*sizeof(double));
    double* DB = (double*)(buffer);
    if (JBOOL == g->data_type){
        _Bool* data_ptr = g->data;
        for (int i = 0; i < N; i++) DB[i] = (double)data_ptr[i];
    }
    else if (JINT == g->data_type){
        int* data_ptr = g->data;
        for (int i = 0; i < N; i++) DB[i] = (double)data_ptr[i];
    }
    else if (JFLOAT == g->data_type){
        float* data_ptr = g->data;
        for (int i = 0; i < N; i++) DB[i] = (double)data_ptr[i];
    }
    else if (JDOUBLE == g->data_type) return;
    else{
        fprintf(stderr, "bad data type\n");
        return;
    }
    g->data_type = JDOUBLE;
    free(g->data); 
    g->data = NULL;
    g->data = buffer;
}

void floatimage(struct cmi_image* g){
    int _row = row(g);
    int _col = col(g);
    int _depth = depth(g);

    int N = _row * _col * _depth;
    void* buffer = (void* )calloc(1, N*sizeof(float));
    float* FL = (float*)(buffer);
    if (JBOOL == g->data_type){
        _Bool* data_ptr = g->data;
        for (int i = 0; i < N; i++) FL[i] = (float)data_ptr[i];
    }
    else if (JINT == g->data_type){
        int* data_ptr = g->data;
        for (int i = 0; i < N; i++) FL[i] = (float)data_ptr[i];
    }
    else if (JFLOAT == g->data_type){
        return; 
    }
    else if (JDOUBLE == g->data_type){
        double* data_ptr = g->data;
        for (int i = 0; i < N; i++) FL[i] = (float)data_ptr[i];
    }
    else{
        fprintf(stderr, "bad data type\n");
        return;
    }
    g->data_type = JFLOAT;
    free(g->data); 
    g->data = NULL;
    g->data = buffer;
}

int sizeof_unit(int32_t type){
    if (JBOOL ==  type) return sizeof(_Bool);
    else if (JINT == type) return sizeof(int);
    else if (JFLOAT == type) return sizeof(float);
    else if (JDOUBLE == type) return sizeof(double);
    else{
        fprintf(stderr, "bad data type\n");
        return -1;
    }
}

void resetimage(struct cmi_image* g){
    int N = pixel_num(g);
    int size_u = sizeof_unit(g->data_type);

    memset(g->data, 0, N*size_u);
}   

