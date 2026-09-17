#include "wfg.h"
#include "avl.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

static void trimLine(char line[])
{
	int i = 0;

	while(line[i] != '\0')
	{
		if (line[i] == '\r' || line[i] == '\n')
		{
			line[i] = '\0';
			break;
		}
		i++;
	}
}

void printContents(FILECONTENTS *f)
{
	int i;
	for ( i = 0; i < f->nFronts; i++)
	{
		printf("Front %d:\n", i+1);
		int j;
		for ( j = 0; j < f->fronts[i].nPoints; j++)
		{
			printf("\t");
			int k;
			for (k = 0; k < f->fronts[i].n; k++)
			{
				printf("%f ", f->fronts[i].points[j].objectives[k]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

FILECONTENTS *readFile(const char filename[])
{
	FILE *fp;
	char line[BUFSIZ];
	int front = 0, point = 0, objective = 0;

	FILECONTENTS *fc = (FILECONTENTS *) std::malloc(sizeof(FILECONTENTS));
	fc->nFronts = 0;
	fc->fronts = NULL;

	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "File %s could not be opened\n", filename);
		exit(EXIT_FAILURE);
	}

	while(fgets(line, sizeof line, fp) != NULL)
	{
		trimLine(line);
		if (std::strcmp(line, "#") == 0)
		{
			front = fc->nFronts;
			fc->nFronts++;
			fc->fronts = (FRONT*)realloc(fc->fronts, sizeof(FRONT) * fc->nFronts);
			fc->fronts[front].nPoints = 0;
			fc->fronts[front].points = nullptr;
		}
		else
		{
			FRONT *f = &fc->fronts[front];
			point = f->nPoints;
			f->nPoints++;
			f->points = (POINT*)std::realloc(f->points, sizeof(POINT) * f->nPoints);
			f->n = 0;
			f->points[point].objectives = NULL;
			f->points[point].tnode = (avl_node_t*) std::malloc(sizeof(avl_node_t));
			char *tok = std::strtok(line, " \t\n");
			do
			{
				POINT *p = &f->points[point];
				objective = f->n;
				f->n++;
				p->objectives = (OBJECTIVE*) std::realloc(p->objectives, sizeof(OBJECTIVE) * f->n);
				p->objectives[objective] = std::atof(tok);
			} while ((tok = std::strtok(NULL, " \t\n")) != NULL);
		}
	}

	fc->nFronts--;
	// for (int i = 0; i < fc->nFronts; i++) fc->fronts[i].n = fc->fronts[i].points[0].nObjectives;
        fclose(fp);
	/* printf("Read %d fronts\n", fc->nFronts);
	   printContents(fc); */
	return fc;
}
