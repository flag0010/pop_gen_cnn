

#include "vector.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//initialization, alloc, free, etc
vector *vectorNew(int startingSize){
	vector *new;
	
	new = malloc(sizeof(vector));
	new->dubs = malloc(startingSize * sizeof(double));
	new->size = startingSize;
	vectorInit(new);
	return(new);
}

vector *vectorClone(vector *v){
	vector *new;
	int i;

	new = malloc(sizeof(vector));
	new->dubs = malloc(v->size * sizeof(double));
	new->size = v->size;
	for(i = 0; i < new->size; i++){
		new->dubs[i] = v->dubs[i];
	}
	return(new);
}

void vectorInit(vector *v){
	int i;
	v->length = 0;
	for(i = 0; i < v->size;i++){
		v->dubs[i] = 0.0;
	}
}

void vectorFree(vector *v){
	free(v->dubs);
	free(v);
}

void vectorRealloc(vector *v, int newSize){
	v->dubs = realloc(v->dubs, newSize * sizeof(double));
	v->size = newSize;
}


//setters / getters
void vectorSetIndexTo(vector *v, int index, double val){
	assert(index < v->length);
	v->dubs[index] = val;
}

void vectorAppend(vector *v, double val){

	if(v->length + 1 >= v->size){
			vectorRealloc(v, v->size * 2) ; // tune?
	}
	v->dubs[v->length] = val;
	v->length += 1;	
}

double vectorGetIndex(vector *v, int index){
	return(v->dubs[index]);
}

