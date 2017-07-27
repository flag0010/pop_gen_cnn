/* stringWrap.c -- methods for stringWrap object
which eases the use of C strings */
#include "stringWrap.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

//initialization, alloc, free, etc
stringWrap *stringWrapNew(int startingSize){
	stringWrap *newWrap;
	int i;
	
	newWrap = malloc(sizeof(stringWrap));
	assert(newWrap);
	newWrap->cString = malloc(startingSize * sizeof(char) + 1); //plus 1 for terminator
	assert(newWrap->cString);
	for(i=0;i<startingSize;i++){newWrap->cString[i] ='\0';}
	newWrap->size = startingSize;
	stringWrapInit(newWrap);
	return(newWrap);
}

void stringWrapFree(stringWrap *aStringWrap){
	free(aStringWrap->cString);
	free(aStringWrap);
}

stringWrap *stringWrapNewFromChar(char *aCString){
	stringWrap *newWrap;
	
	newWrap = stringWrapNew(strlen(aCString));
	stringWrapAppendChar(newWrap, aCString);
	return(newWrap);
}

void stringWrapRealloc(stringWrap *aStringWrap, int newSize){
	aStringWrap->cString = realloc(aStringWrap->cString, newSize * sizeof(char));
	assert(aStringWrap->cString);
	aStringWrap->size = newSize;
}

//stringWrapInit-- sets length to zero, and init.s cString
void stringWrapInit(stringWrap *aStringWrap){
	aStringWrap->length = 0;
	aStringWrap->cString[0] = '\0';
}

int stringWrapLength(stringWrap *aStringWrap){
	int i=0;
	while(aStringWrap->cString[i] != '\0' && aStringWrap->cString[i] != '\n'){
		i++;
	}
	return(i);
}
		

//this is the key function that allows for auto-realloc of strings
// might need tuning
void stringWrapAppendChar(stringWrap *destWrap,char *aCString){
	int i, len;

	len = strlen(aCString);
	if(destWrap->length + len + 1 >= destWrap->size){
		stringWrapRealloc(destWrap, max(destWrap->length + len, destWrap->size * 2) + 1); // tune?
	}
	for (i = 0; i < len && aCString[i] != '\0'; i++){
		//beware that this upper cases the characters
		destWrap->cString[destWrap->length + i] = toupper(aCString[i]);
	}
	//add termination
	destWrap->length += i;
	destWrap->cString[destWrap->length] = '\0'; 

}

void stringWrapAppendSingle(stringWrap *destWrap, char aChar){
	int i, len;

	len = 1;
	if(destWrap->length + len + 1 >= destWrap->size){
		stringWrapRealloc(destWrap, max(destWrap->length + len, destWrap->size * 2) + 1); // tune?
	}
	for (i = 0; i < len && aChar != '\0'; i++){
		destWrap->cString[destWrap->length + i] = aChar;
	}
	//add termination
	destWrap->length += i;
	destWrap->cString[destWrap->length] = '\0';
}

void stringWrapAppend(stringWrap *dst, stringWrap *src){
	stringWrapAppendChar(dst, src->cString);
}


void stringWrapChomp(stringWrap *dst){
	if(dst->cString[dst->length - 1] == '\n')
		dst->cString[dst->length - 1] = '\0';
}

//
// I/O -- this will need work
//

void stringWrapPrintf(stringWrap *aStringWrap){
	printf("%s",aStringWrap->cString);
}

//stringWrapReadline puts the next line of a file into a stringWrap
int stringWrapReadLine(stringWrap *aStringWrap, FILE *aFile){
	char buf[BUFFER];
		
	stringWrapInit(aStringWrap);

	do
	{
		if (fgets(buf, BUFFER, aFile) == NULL)
			return(EOF);
		else{
			//buf[BUFFER-2] ='\n';
			stringWrapAppendChar(aStringWrap,buf);
		}
	}
	while(strchr(buf, '\n') == NULL);

	return(0);	
}

// string utilities

//slice returns a new stringWrap with the substring
stringWrap *stringWrapSlice(stringWrap *aStringWrap, int start, int end){
	char tmpChar[end-start+2];
	int i;
		
	for(i=0;i<=end-start;i++){
		tmpChar[i] = aStringWrap->cString[i+start];
	}
	tmpChar[i]='\0';
	return(stringWrapNewFromChar(tmpChar));
}

int stringWrapHasChar(stringWrap *aStringWrap, char aChar){
	int i,flag = 0;
	i = 0;
	while(i < aStringWrap->length && flag == 0){
		if (aStringWrap->cString[i] == aChar){
			flag = 1;
		}
		i++;
	}
	return(flag);
}

int stringWrapCountChar(stringWrap *aStringWrap, char aChar){
	int i, count;
	count = 0;
	for(i = 0; i < aStringWrap->length; i++){
		if (aStringWrap->cString[i] == aChar){
			count+=1;
		}
	}
	return(count);
}

//returns first char that is not aChar
stringWrap *stringWrapFindOther(stringWrap *aStringWrap, char aChar){
	char tmpChar[2];
	int i;
		
	for(i=0;i< aStringWrap->length;i++){
		if (aStringWrap->cString[i] != aChar){
		tmpChar[0] = aStringWrap->cString[i];
		}
	}
	tmpChar[1]='\0';
	return(stringWrapNewFromChar(tmpChar));
}

// Basics... from K&R
/* reverse:  reverse string s in place */
void reverse(char s[]){
    int i, j;
    char c;

    for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

/* itoa:  convert n to characters in s */
void itoa(int n, char s[]){
    int i, sign;

    if ((sign = n) < 0)  /* record sign */
        n = -n;          /* make n positive */
    i = 0;
    do {       /* generate digits in reverse order */
        s[i++] = n % 10 + '0';   /* get next digit */
    } while ((n /= 10) > 0);     /* delete it */
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
}

