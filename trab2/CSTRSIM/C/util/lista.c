
#include <stdio.h>
#include <stdlib.h>
#include "lista.h"

fault *newNode(int c, double nv, double it, double t){
	fault *new;
	new = (fault*)malloc(sizeof(fault));
	
	new->fcode = c;
	new->newValue = nv;
	new->inst = it;
	new->tau = t;
	new->next = NULL;
	
	return new;
}

fault *addNode(fault *root, fault *new){
	fault *aux = root;
	
	if(root == NULL)
		return new;
	
	while(aux->next != NULL)
		aux = aux->next;
	
	aux->next = new;
	
	return root;
}

fault *addFault(fault *root, int c, double nv, double it, double t){
	fault *aux = root;
	
	if(root == NULL)
		return newNode(c, nv, it, t);
	
	while(aux->next != NULL)
		aux = aux->next;
	
	aux->next = newNode(c, nv, it, t);
	
	return root;
}

fault *delNode(fault *root, int c){
	fault *aux = root;
	fault *prv = NULL;
	
	if(root == NULL)
		fprintf(stderr, "Nao existem falhas\n");
	
	while(aux != NULL && aux->fcode != c){
		prv = aux;
		aux = aux->next;
	}
	
	if(aux == NULL){
		fprintf(stderr, "Falha nao encontrada\n");
		return root;
	}
	
	if(prv == NULL)
		root = aux->next;
	else
		prv->next = aux->next;
	
	free(aux);
	
	return root;
}

fault *delLast(fault *root){
	fault *aux = root;
	fault *prv = NULL;
	
	while(aux->next != NULL){
		prv = aux;
		aux = aux->next;
	}
	
	prv->next = NULL;
	free(aux);
	
	return root;
}

void printFaultList(fault *lf){
	fault *aux = lf;
	
	while(aux != NULL){
		fprintf(stdout, "%d\t%f\t%f\t%f\t\n", aux->fcode, aux->newValue, aux->inst, aux->tau);
		aux = aux->next;
	}
}
