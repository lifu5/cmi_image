#include <stdio.h>
#ifndef _UTILS_H
#define _UTILS_H


#include "rbtree.h"
#include <stdlib.h>
#include <unistd.h>
#include <string.h>


struct point2d{
    float xx, yy;
};
typedef struct point2d point2d;

struct point3d{
    float xx, yy, zz;
};
typedef struct point3d point3d;

int get_filesize(char* filename);

struct mynode{
    double dist;
    point2d coord;
    struct mynode* next;
};

typedef struct mynode mynode;
mynode* initnode(double dist, point2d coord);


struct cmilist{
    mynode* head;
    //mynode* curr;
    mynode* tail;
    int len;
};
typedef struct cmilist cmilist;

cmilist* initlist(mynode*);
// provide limit function
void listinsert(cmilist* list, mynode* newnode);

void listsort(cmilist* list);  

//map part is underconstucted
struct map{
    struct rb_node node;
    float dist; //key: sort by dist
    point2d* coord; // value
};
typedef struct map map_t;
typedef struct rb_root root_t;
typedef struct rb_node rb_node_t;

map_t* get_coord(root_t* root, float dist);

int insert(root_t* root, float key, point2d* coord);

#endif
