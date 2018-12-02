#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _CONFIG_H
#define _CONFIG_H

#define MAX_BUF_LEN 1024
#define MAX_KY


void trim(char* strIn, char* strOut);

void getValue(char * keyAndValue, char * key, char * value);

void readCFG(
        const char *filename/*in*/,
        const char *key/*in*/,
        const char **value/*out*/);


#endif
