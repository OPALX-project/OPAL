#ifndef _WFG_H_
#define _WFG_H_

#include <cstdio>
#include <cstdlib>
#include <string>

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
      struct avl_node_t * tnode;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
} FRONT;

typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

FILECONTENTS *readFile(const char[]);

extern void printContents(FILECONTENTS *);

#endif
