
#include "kd_alloc.h"

void * kd_alloc(int type_size, int kdim, int *dims)
{
    int i,j;
    char * ptr;
    int dimv, size;
    int memchunk;
    
    char ** current_address;
    char * forward_address;

    memchunk = 0;
    dimv = 1;
    size = sizeof(char *);
    for (i=0; i<kdim; i++)
    {
        if (i == kdim-1) size = type_size;
        dimv *= dims[i];
        memchunk += dimv*size;
    }

    ptr = (char *) malloc(memchunk);

    dimv = 1;
    current_address = (char **) ptr;
    size = sizeof(char *);

    for (i=0; i<kdim-1; i++)
    {
        
        dimv *= dims[i];
        forward_address = (char *) (current_address + dimv);

        if (i >= kdim-2) size = type_size;

        for (j=0; j<dimv; j++)
            current_address[j] = forward_address + j*dims[i+1]*size;

        current_address = (char **) forward_address;

    }

    return ptr;

}


void * kd_alloc2( int typesize, int kdim, ... )
{
    int i;
    char * ptr;
    int dims[kdim];
    va_list dim_list;

    va_start( dim_list, kdim );

    for (i=0; i<kdim; i++)
        dims[i] = va_arg( dim_list, int );

    va_end(dim_list);
    ptr = (char *) kd_alloc( typesize, kdim, dims );
    return ptr;
}

#include <stdio.h>

void debug(char ** ptr, int kdim, int * dims)
{
    int i,j;
    int dimv = 1;
    int offset = 0;

    for (i=0; i<kdim-1; i++)
    {
        dimv *= dims[i];
        printf("\n dim = %d dims = %d\n", i+1, dims[i]);

        for (j=0; j<dimv; j++)
            printf("%d: \t%p -> %p\n", j, 
                                       (void *)(ptr+j+offset), 
                                       (void *) *(ptr+j+offset));

        offset += dimv;

    }

}

/*
int main()
{

    double ** array;
    int dims[2];
    dims[0] = 4;
    dims[1] = 3;

    array = (double **) kd_alloc2(sizeof(double), 2, dims[0], dims[1]);

    debug((char **) array, 2, dims);

    int i,j;
    for (i=0; i<dims[0]; i++)
    for (j=0; j<dims[1]; j++)
        array[i][j] = i*dims[1] + j;

    for (i=0; i<dims[0]; i++)
    for (j=0; j<dims[1]; j++)
        printf("%f\n", (float) array[i][j]);

    free(array);

    return 0;

}
*/
