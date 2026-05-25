//Funcoes para manipulacao do arquivo de configuracao

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lista.h"
#define M 100

void leComentario(int n, char *str, FILE *arq){
	int i;
	str[0] = '\0';
	for(i=1; i<=n; i++)
		fgets(str, M, arq);
}

fault *transformaLinha(fault *ft, char *str){
	int i=0, j;
	char str1[M];
	str1[0] = '\0';
	
	j=0;
	while(str[i] != 32 && str[i] != 9){
		
		if(str[i] >= '0' && str[i] <= '9' || str[i] == 46){
			str1[j] = str[i];
		}
		i++;
		j++;
	}
	str1[j] = '\0';
	int code = -1;
	code = atoi(str1);

	i++;
	j=0;
	str1[0] = '\0';
	while(str[i] != 32 && str[i] != 9){
		
		if(str[i] >= '0' && str[i] <= '9' || str[i] == 46 || str[i] == 45){
			str1[j] = str[i];
		}
		i++;
		j++;
	}
	str1[j] = '\0';
	double nv = -1;
	nv = atof(str1);
	
	i++;
	j=0;
	str1[0] = '\0';
	while(str[i] != 32 && str[i] != 9){
		
		if(str[i] >= '0' && str[i] <= '9' || str[i] == 46){
			str1[j] = str[i];
		}
		i++;
		j++;
	}
	str1[j] = '\0';
	double inst = -1;
	inst = atof(str1);
	
	i++;
	j=0;
	str1[0] = '\0';
	while(str[i] != '\0'){
		
		if(str[i] >= '0' && str[i] <= '9' || str[i] == 46){
			str1[j] = str[i];
		}
		i++;
		j++;
	}
	str1[j] = '\0';
	double tau = -1;
	tau = atof(str1);
	
	return addFault(ft, code, nv, inst, tau);
}

int transformaListaFalhas(fault *listaFalhas, int *falhas, double *nv, double *it, double *tu){
	
	fault *aux = listaFalhas;
	int i = 0;

	while(aux != NULL){
		falhas[aux->fcode] = aux->fcode;
		nv[aux->fcode] = aux->newValue;
		it[aux->fcode] = aux->inst;
		tu[aux->fcode] = aux->tau;
		//printf("transformaListaFalhas> i=%3d fault=%3d  tau=%.2e\n",i,aux->fcode,aux->tau);
		aux = aux->next;
		i++;
	}

	return i;
}

unsigned int leInteiro(FILE *arq, char *str){
	do{
		fgets(str, M, arq);
	}while(str[0] == 35 || str[0] == 32 || str[0] == 10 || str[0] == 9);
	
	return atoi(str);
}

double leReal(FILE *arq, char *str){
	do{
		fgets(str, M, arq);
	}while(str[0] == 35 || str[0] == 32 || str[0] == 10 || str[0] == 9);

	return atof(str);
}

fault *leFalhas(FILE *arq, fault *list, char *str){
	int ln = 0, esc = 0;

	do{
		fgets(str, M, arq);
		//if(str[0] == 35 && str[1] == 35)
		//	break;
		if(ln++ > 99){
			//fprintf(stderr,"ERROR: Could not find the end of file [##]\n\n");
			//return NULL;
			esc = 1;
			break;
		}
	}while(str[0] == 35 || str[0] == 32 || str[0] == 10 || str[0] == 9);
	
	if(esc == 1)
		return NULL;

	ln = strlen(str);
	str[0] = '\0';
	
	fseek(arq, -ln, SEEK_CUR);

	while(!feof(arq) && esc == 0){
		str[0] = '\0';
		if(str[0] == 35 || str[0] == 32 || str[0] == 10 || str[0] == 9){
			fgets(str, 256, arq);
		}
		else{
			fgets(str, 256, arq);
			list = transformaLinha(list, str);
		}
	}
	return list;
}

int leArquivoFalhas(FILE *arq, unsigned int *sd, double *its, int *git,
		int *falhas, double *nv, double *inst, double *tau){

	int i;
	
	fault *listaFalhas = NULL;
	char str[M];
	
	*sd = leInteiro(arq, str);
	*its = leReal(arq, str);
	*git = leInteiro(arq, str);//printf("\n\n%u \n\n", *git);
	if(*git < 1) *git = 1;
	str[0] = '\0';

	listaFalhas = leFalhas(arq, listaFalhas, str);
	fclose(arq);

	if(listaFalhas == NULL){
		return 0;
	}

	//printFaultList(listaFalhas);	
	int qtf = transformaListaFalhas(listaFalhas, falhas, nv, inst, tau);
	//printf("Quantidade de falhas no arquivo: %d\n", qtf);

    /*for(i=0;i<49;i++)
        printf("leArquivoFalhas> nf=%2d i=%2d fl=%2d nv=%7.2lf inst=%7.2lf tau=%7.2lf\n",
	            qtf,i,falhas[i],nv[i],inst[i],tau[i]);       /**/

	while( listaFalhas != NULL )	{
		fault *ant = listaFalhas;
		listaFalhas = listaFalhas->next;
		free( ant );
	}
	return qtf;
}

void initVets(int *fl, double *nv, double *it, double *tu, int n){
	int i;
	
	for(i=0; i<n; i++){
		fl[i] = 0;
		nv[i] = 0;
		it[i] = 0;
		tu[i] = 0;
	}
}

void printMsg(FILE *arq, int tp){
	
	if(tp == 1){
		fprintf(arq, "#\n# 1. Title: CSTR Simulator output\n#\n");
		fprintf(arq, "# 2. Source:\n#      CSTR Simulator\n#\n# 3. Number of Instances: See the file used to configure the simulator\n#\n");
		fprintf(arq,"# 4. Number of Attributes: 38 numeric and the class\n#\n");
		fprintf(arq, "# 5. Attribute Information:\n#    1. Tank head\n#    2. Outlet tank head\n#    3. Head after pipe 1\n");
		fprintf(arq, "#    4. Head after valve\n#    5. Head after pipe 2\n#    6. Leak head\n#    7. Tank Inlet flow rate\n");
		fprintf(arq, "#    8. Flow rate in pipe 1\n#    9. Flow rate in pipe 2\n#    10. Leak flow rate\n#    11. Reactor temperature\n");
		fprintf(arq, "#    12. Cooling system inlet head\n#    13. Head after pipe 1\n#    14. Head after valve\n#    15. Head after jacket obstruction\n");
		fprintf(arq, "#    16. Head after pipe 2\n#    17. Leak to environment head\n#    18. Leak to tank head\n");
		fprintf(arq, "#    19. Cooling system inlet flow rate (pipe 1)\n#    20. Flow rate in pipe 2\n#    21. Leak to environment flow rate\n");
		fprintf(arq, "#    22. Leak to tank flow rate\n#    23. Coolant inlet temperature\n#    24. Jacket temperature\n#    25. Generation rate of product B\n");
		fprintf(arq, "#    26. Generation rate of product C\n#    27. Feed temperature\n#    28. Feed concentration\n#    29. Concentration of reactant A\n#    30. Concentration of product B\n");
		fprintf(arq, "#    31. Concentration of product C\n#    32. Level controller output\n#    33. Coolant controller output\n#    34. Coolant setpoint\n");
		fprintf(arq, "#    35. Inventory\n#    36. Mol balance\n#    37. Cooling water head loss\n#    38. Effluent head loss\n#    39. Class\n#\n");
		fprintf(arq, "# 6. Missing Attribute Values: None\n#\n#\n#\n");
		fprintf(arq, "38\n");
	}

	if(tp == 2){
		fprintf(arq, "#\n# 1. Title: CSTR Simulator output\n#\n");
		fprintf(arq, "# 2. Source:\n#      CSTR Simulator\n#\n# 3. Number of Instances: See the file used to configure the simulator\n#\n");
		fprintf(arq,"# 4. Number of Attributes: 18 numeric and the class\n#\n");
		fprintf(arq, "# 5. Attribute Information:\n#    1. Feed concentration\n#    2. Feed flowrate\n#    3. Feed temperature\n");
		fprintf(arq, "#    4. Reactor level\n#    5. Product A concentration\n#    6. Product B concentration\n#    7. Reactor temperature\n");
		fprintf(arq, "#    8. Coolant flowrate\n#    9. Product flowrate\n#    10. Coolant inlet temperature\n#    11. Coolant inlet pressure\n");
		fprintf(arq, "#    12. Level controller output\n#    13. Coolant controller output\n#    14. Coolant setpoint\n#    15. Inventory\n");
		fprintf(arq, "#    16. Mol balance\n#    17. Cooling water head loss\n#    18. Effluent head loss\n#    19. Class\n#\n");
		fprintf(arq, "# 6. Missing Attribute Values: None\n#\n#\n#\n");
		fprintf(arq, "18\n");	
	}

}
