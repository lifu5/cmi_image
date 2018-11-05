#include <stdlib.h>
#include "cmiutils.h"

int get_filesize(char* filename){
    FILE* fp =fopen("filename", "rb");
    int size;
    if(fp == NULL){ //open files failed 
        perror("can not open the file");
        return -1;
    }
    fseek(fp, 0, SEEK_END);
    size=ftell(fp);
    fclose(fp);
    return size;
}




