//Tallison 	Carvalho Moraes
#include <stdio.h>
#include "params.h"
#include "gaussrand.h"
#include "nrutil.h"
#include "ludcmp.h"

//Metodo de Newton-Raphson para calculo das variaves.
//Tem como parametros os vetores de fluxos e de cabecais, as aberturas das
//valvulas, as aberturas dos vazamentos e a quantidade de variaves(fluxos + cabecais)
double *newtRap(double *q, double *h, double *k, double mvalv, double *mvo, double ktb1, 
	            double mvalvc, int n){	
	
	//func - vetor de funções
	//jac - Jacobiana, com as derivadas das funções
	double *func, **jac, **auxjac, *r, *in, d;
	
	//n - tamanho do vetor de raizes iniciais
	//indx - ponteiro utilizado 
	int *indx;
	
	//Parametros do problema
	//double k[13];
	double kefl, kefv, keflc;
	double kefv1c, kefv2c;
	k[1] = ktb1;
	
	//Calculo das constantes de efluente, valvula e vazamento do reator
	kefl = 1/(2*g*A*A);
	k[3] = 1/(mvalv*mvalv);
	k[2] = 1/(mvo[0]*mvo[0]);
	kefv = 1/(2*g*mvo[0]*mvo[0]);
	
	//Calculo das constantes de efluente, valvula e vazamento do resfriamento
	k[6] = 1/(mvalvc*mvalvc);
	k[8] = 1/(mvo[1]*mvo[1]);
	k[7] = 1/(mvo[2]*mvo[2]);
	keflc = 1/(2*g*Ac*Ac);	// Ac = Cooling system cross-section area
	k[9] = (mvo[3] >= 1)? 0.0 : (1/(mvo[3]*mvo[3]));
	kefv1c = 1/(2*g*mvo[1]*mvo[1]);
	kefv2c = 1/(2*g*mvo[2]*mvo[2]);
	
	//Inicializacao dos vetores e da matriz utilizando a biblioteca
	//Numerical Recipes in C
	func = dvector(1,n);
	jac = dmatrix(1,n,1,n);
	r = dvector(1,n);
	in = dvector(1,n);
	indx = ivector(1,n);
	
	//Loop principal do metodo de Newton-Raphson
	int it = 0;
	while(it < 10){
		
		double q2sqr = q[2]*q[2]; double q3sqr = q[3]*q[3]; double q4sqr = q[4]*q[4];
		//Matriz de funcoes	
		//Saida do reator
		func[1] = h[2]-h[1]-h[0]+k[11]*q[2]+k[12]*q2sqr;	//Bomba
		func[2] = h[3]-h[2]+ktb1*q2sqr; 		        	//Tubo 1
		func[3] = h[4]-h[3]+k[3]*q4sqr;			        	//Valvula
		func[4] = h[6]-h[3]+k[2]*q3sqr;			        	//Vazamento
		func[5] = h[5]-h[4]+k[4]*q4sqr;			        	//Tubo 2
		func[6] = -h[5]+kefl*q4sqr;			            	//Efluente
		func[7] = q[4]+q[3]-q[2];			            	//Soma dos fluxos
		func[8] = -h[6]+kefv*q3sqr;				            //Efluente do vazamento
		
		double q5sqr = q[5]*q[5]; double q6sqr = q[6]*q[6]; double q7sqr = q[7]*q[7]; double q8sqr = q[8]*q[8];
		//Sistema de resfriamento
		func[17] = -q[5]+q[8]+q[6]+q[7];			//Soma dos fluxos
		func[9]  = -h[7]+h[8]+k[5]*q5sqr;			//Tubo 1
		func[10] = -h[8]+h[9]+k[6]*q5sqr;			//Valvula
		func[13] = -h[9]+h[13]+k[7]*q6sqr;			//Vazamento para o reator
		func[11] = -h[9]+h[12]+k[8]*q7sqr;			//Vazamento para o ambiente
		func[12] = -h[10]+h[11]+k[10]*q8sqr;		//Tubo 2
		func[14] = -h[11]+keflc*q8sqr;				//Efluente do Tubo2
		func[15] = -h[12]+kefv1c*q7sqr;				//Efluente do vazamento para o ambiente
		func[16] = -h[13]+kefv2c*q6sqr;				//Efluente do vazamento para o reator
		func[18] = -h[9]+h[10]+k[9]*q8sqr;			//Obstrucao na jaqueta de resfriamento
		
		//Zera os elementos da matriz
		initmatrix(jac, n);
		
		//Jacobiana
		//Saida do reator
		jac[1][1] = 1;
		jac[1][5] = 1*k[11]+2*k[12]*q[2];
		
		jac[2][1] = -1;
		jac[2][2] = 1;
		jac[2][5] = 2*ktb1*q[2];
		
		jac[3][2] = -1;
		jac[3][3] = 1;
		jac[3][6] = 2*k[3]*q[4];
		
		jac[4][2] = -1;
		jac[4][7] = 2*k[2]*q[3];
		jac[4][8] = 1;
		
		jac[5][3] = -1;
		jac[5][4] = 1;
		jac[5][6] = 2*k[4]*q[4];
		
		jac[6][4] = -1;
		jac[6][6] = 2*kefl*q[4];
		
		jac[7][5] = -1;
		jac[7][6] = 1;
		jac[7][7] = 1;
		
		jac[8][7] = 2*kefv*q[3];
		jac[8][8] = -1;
		
		//Sistema de resfriamento
		jac[9][9] = 1;
		jac[9][15] = 2*k[5]*q[5];
		
		jac[10][9] = -1;
		jac[10][10] = 1;
		jac[10][15] = 2*k[6]*q[5];
		
		jac[11][10] = -1;
		jac[11][13] = 1;
		jac[11][17] = 2*k[8]*q[7];
		
		jac[12][11] = -1;
		jac[12][12] = 1;
		jac[12][16] = 2*k[10]*q[8];
		
		jac[13][10] = -1;
		jac[13][14] = 1;
		jac[13][18] = 2*k[7]*q[6];
		
		jac[14][12] = -1;
		jac[14][16] = 2*keflc*q[8];
		
		jac[15][13] = -1;
		jac[15][17] = 2*kefv1c*q[7];
		
		jac[16][14] = -1;
		jac[16][18] = 2*kefv2c*q[6];
		
		jac[17][15] = -1;
		jac[17][16] = 1;
		jac[17][17] = 1;
		jac[17][18] = 1;
		
		jac[18][10] = -1;
		jac[18][11] = 1;
		jac[18][16] = 2*k[9]*q[8];
		
		//Inicialização do vetor de raizes iniciais
		r[1] = h[2];
		r[2] = h[3];
		r[3] = h[4];
		r[4] = h[5];
		r[5] = q[2];
		r[6] = q[4];
		r[7] = q[3];
		r[8] = h[6];
		r[9] = h[8];
		r[10] = h[9];
		r[11] = h[10];
		r[12] = h[11];
		r[13] = h[12];
		r[14] = h[13];
		r[15] = q[5];
		r[16] = q[8];
		r[17] = q[7];
		r[18] = q[6];
		
		//Utilizando a resolucao pela formula Ax=b sugerida pela documentacao
		//Multipliacao da Jacobiana pelo vetor de raizes. (Funcao se encontra em ludcmp.c)
		//Resultado é retornado em in
		matrixvector(jac,r,n,in);
		
		//Subtração do vetor resultado da operação anterior pelo vetor (Funcao se encontra em ludcmp.c)
		//de funções com a raizes atuais
		subvector(in,func,n,in);
		
		//Variavel necessaria pois as funcoes abaixo destroem a matriz
		//original
		auxjac = jac;
		
		//Computa o resultado do sistema linear utilizando 
		//decomposicao LU da matriz (Funcoes se encontram em ludcmp.c)
		ludcmp(auxjac,n,indx,&d);
		lubksb(auxjac,n,indx,in);
		
		//Atualização das raizes
		h[2] = in[1];
		h[3] = in[2];
		h[4] = in[3];
		h[5] = in[4];
		q[2] = in[5];
		// fprintf(stdout, "Q4: Before=%.5e After=%.5e\n", q[4], in[6]);
		q[4] = in[6];
		q[3] = in[7];
		h[6] = in[8];
		h[8] = in[9];
		h[9] = in[10];
		h[10] = in[11];
		h[11] = in[12];
		h[12] = in[13];
		h[13] = in[14];
		q[5] = in[15];
		q[8] = in[16];
		q[7] = in[17];
		q[6] = in[18];	
		
		it++;
	}
	
	free_dvector(func, 1, n);
	free_dmatrix(jac, 1, n, 1, n);
	free_dvector(r, 1, n);
	free_dvector(in, 1, n);
	free_ivector(indx, 1, n);

	return k;
}

