//CSTR Simulator

#include <stdio.h>
#include <stdlib.h>
#include "runkut.h"


int main(int argc, char *argv[]){

	printf("\n\n");
	printf("****************************************\n");
	printf("************ CSTR Simulator ************\n");
	printf("****************************************\n");
	
	procQuimRK(argc, argv);

	return 0;
}
