#include <stdlib.h>
#include "cmiutils.h"

int get_filesize(char* filename){
    FILE* fp =fopen("filename", "rb");
    int size;
    if(fp == NULL){ //open files failed 
        perror("can not open the file\n");
        return -1;
    }
    fseek(fp, 0, SEEK_END);
    size=ftell(fp);
    fclose(fp);
    return size;
}

cd_nodelist* initlist(cd_node *node){
    cd_nodelist* p;
    p = (cd_nodelist*)malloc(sizeof(cd_nodelist));
    p->head = node;
    p->tail = node;
    p->tail->next = NULL;
    p->len = 1;
    return p;
}

cd_node* initnode(double dist, point2d coord, double weight){
    cd_node* p;
    p = (cd_node*)malloc(sizeof(cd_node));
    p->dist = dist;
    p->weight = weight;
    p->coord = coord;
    p->next = NULL;
    return p;
}


int_node* intnode(int val){
    int_node* p;
    p = (int_node*)malloc(sizeof(int_node));
    p->val = val;
    p->next = NULL;
    return p;
}

void grid_node_insert(grid_info* gridp, int_node* newnode){
    if (0 == gridp->len){
        gridp->head= newnode;
        gridp->tail = newnode;
        gridp->tail->next = NULL;
        gridp->len = 1;
        return;
    }
    else{
    gridp->tail->next = newnode;
    gridp->tail = gridp->tail->next;
    gridp->tail->next = NULL;
    gridp->len++;
    }
}


void listinsert(cd_nodelist* list, cd_node* newnode){
    if (0 == list->len){
        list->head= newnode;
        list->tail = newnode;
        list->tail->next = NULL;
        list->len = 1;
        return;
    }
    else{
    list->tail->next = newnode;
    list->tail = list->tail->next;
    list->tail->next = NULL;
    list->len++;
    }
}

void listsort(cd_nodelist* list){
    //bubble sort
    cd_node *p,*q;
    double tmp_d;
    point2d tmp_p;
    for (p = list->head; p!=NULL; p = p->next){
        for (q = p->next; q!=NULL; q = q->next){
            if ( p->dist> q->dist){
                tmp_d = q->dist; tmp_p = q->coord;
                q->dist = p->dist; q->coord = p->coord;
                p->dist = tmp_d; p->coord = tmp_p;
            }
        }
    }
}


//underconstruction
map_t* get_coord(root_t* root, float dist){
    rb_node_t* node = root->rb_node;
    while(node){
        map_t* data = container_of(node, map_t, node);
        if (data->dist < dist){
            node = node->rb_left;
        }else if(data->dist > dist){
            node = node->rb_right;
        }else{
            return data;
        }
    }
    return NULL;
}

int insert(root_t* root, float key, point2d* coord){
    map_t *data = (map_t*)malloc(sizeof(map_t));
    data->dist = key;
    data->coord = (point2d*)malloc(sizeof(point2d));
    data->coord = coord;

    //rb_node_t **new_node = &(root->rb_node), *parent = NULL;
    //while(*new_node){
    //    map_t *this_node = container_of(*new_node, map_t, node);
    //}
    return 0;
}
