/* vector.h -- methods for vector object */

#ifndef VECTOR_INC
#define VECTOR_INC

#include "stdio.h"

// stringWrap Object Definition 
typedef struct{
  int length;			// length of string 
  double *dubs;			// true doubles underneath
  int size;				// current alloc'd Size
} vector;

vector *vectorNew(int startingSize);
vector *vectorClone(vector *v);
void vectorFree(vector *v);
void vectorInit(vector *v);
void vectorSetIndexTo(vector *v, int index, double val);
void vectorAppend(vector *v, double val);
double vectorGetIndex(vector *v, int index);










#endif
