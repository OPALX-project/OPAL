/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
  ========================================================================
  NSGA2

  Implements data exchange trough files.

  file: nsga2_io.c
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  last change: $date$
  ========================================================================
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "nsga2.h"


int read_pop(char *filename, pop *pp, int size, int dim)
/* Reads individuals from file into pop */
{
    int i, j;
    int entries = 0;
    char tag[4];
    FILE *fp;
    int result;

    assert(dim >= 0);
    assert(pp != NULL);

    fp = fopen(filename, "r");
    assert(fp != NULL);

    result = fscanf(fp, "%d", &entries);
    if (entries == 0) /* file has not been written yet */
    {
       return (1);  /* signalling that reading failed */
    }
    assert(entries == size * (dim + 1));

    for (j = 0; j < size; j++)
    {
	/* reading index of individual */
	result = fscanf(fp, "%d", &(pp->ind_array[j]->index));

	for (i = 0; i < dim; i++)
	{
	    /* reading objective values of ind */
	    result = fscanf(fp, "%le", &(pp->ind_array[j]->f[i]));
	    assert(pp->ind_array[j]->f[i] >= 0);
	    if (result == EOF) /* file not completely written */
	    {
		fclose(fp);
		return (1); /* signalling that reading failed */
	    }
	}
    }

    /* after all data elements: "END" expected */
    result = fscanf(fp, "%s", tag);
    if (strcmp(tag, "END") != 0)
    {
	fclose(fp);
	return (1);  /* signalling that reading failed */
    }
    else /* "END" ok */
    {
	fclose(fp);

        /* delete file content if reading successful */
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "0");
        fclose(fp);

        return (0);  /* signalling that reading was successful */
    }
}



void write_pop(char* filename, pop* pp, int size)
/* Writes a pop or PISA_to a given filename. */
{
     int i;
     FILE *fp;

     assert(0 <= size && size <= pp->size);

     fp = fopen(filename, "w");
     assert(fp != NULL);

     fprintf(fp, "%d\n", size); /* number of elements */

     for (i = 0; i < size; i++)
     {
          fprintf(fp, "%d\n", pp->ind_array[i]->index);
     }

     fprintf(fp, "END");
     fclose(fp);
}



int check_file(char* filename)
{
     int control_element = 1;
     FILE *fp;

     fp = fopen(filename, "r");
     assert(fp != NULL);
     if (fscanf(fp, "%d", &control_element) != 1) {
         fprintf(stderr, "Couldn't read file '%s'.", filename);
     }

     fclose(fp);

     if(0 == control_element)
          return (0); /* file is ready for writing */
     else
          return (1); /* file is not ready for writing */
}
