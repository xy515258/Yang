
#ifdef __cplusplus 
extern "C" {
#endif

#ifndef KD_ALLOC_H
#define KD_ALLOC_H

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

void * kd_alloc(int type_size, int kdim, int *dims);
void * kd_alloc2( int typesize, int kdim, ... );
void debug(char ** ptr, int kdim, int * dims);

#endif

#ifdef __cplusplus
}
#endif

