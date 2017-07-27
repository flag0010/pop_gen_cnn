/* stringWrap.c -- methods for stringWrap object
which eases the use of C strings */

#ifndef STRINGWRAP_INC
#define STRINGWRAP_INC

#include "stdio.h"

//#define BUFFER 300000000 /*enough for normal lines */
#define BUFFER 10000
#define SMALL	100
#define LARGE	10000
#define max(a, b) (a > b ? a : b)

// stringWrap Object Definition 
typedef struct{
  int length;			// length of string 
  char *cString;		// true cString underneath
  int size;				// current alloc'd Size
} stringWrap;

stringWrap *stringWrapNew(int startingSize);
void stringWrapFree(stringWrap *aStringWrap);
stringWrap *stringWrapNewFromChar(char *aCString);
void stringWrapRealloc(stringWrap *aStringWrap, int newSize);
void stringWrapAppendChar(stringWrap *destWrap,char *aCString);
void stringWrapAppend(stringWrap *dst, stringWrap *src);
void stringWrapPrintf(stringWrap *aStringWrap);
void stringWrapInit(stringWrap *aStringWrap);
int stringWrapReadLine(stringWrap *aStringWrap, FILE *aFile);

stringWrap *stringWrapSlice(stringWrap *aStringWrap, int start, int end);
int stringWrapLength(stringWrap *aStringWrap);
int stringWrapHasChar(stringWrap *aStringWrap, char aChar);
int stringWrapCountChar(stringWrap *aStringWrap, char aChar);
void stringWrapAppendSingle(stringWrap *destWrap, char aChar);
void stringWrapChomp(stringWrap *dst);
stringWrap *stringWrapFindOther(stringWrap *aStringWrap, char aChar);


void itoa(int n, char s[]);
void reverse(char s[]);
#endif

