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

struct coord_dist_node{
    double dist;
    point2d coord;
    double weight;
    struct coord_dist_node* next;
};

typedef struct coord_dist_node cd_node;


struct cd_nodelist{
    cd_node* head;
    //mynode* curr;
    cd_node* tail;
    int len;
};

struct int_node{
    int val;
    struct int_node* next;
};
typedef struct int_node int_node; 

struct grid_info{
    int_node* head;
    int_node* tail;
    int index_m; //parallel
    double dist1; //projection on detect plane
    int len;
};
typedef struct grid_info grid_info;
typedef struct cd_nodelist cd_nodelist;

cd_nodelist* initlist(cd_node*);


// only provide limit function for list structure
cd_node* initnode(double dist, point2d coord, double weight);
int_node* intnode(int val);
void grid_node_insert(grid_info* list, int_node* newnode);
void listinsert(cd_nodelist* list, cd_node* newnode);

void listsort(cd_nodelist* list);  

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
