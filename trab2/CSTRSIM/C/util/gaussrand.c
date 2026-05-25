//Gerador de numeros aleatorios com distribuicao normal
//Metodo discutido por Knuth, criado por Marsaglia

#include <stdlib.h>
#include <math.h>

double gauss_rand(){
	
	double U1, U2, W, mult, mu, sigma;
	static double X1, X2;
	static int call = 0;
	mu = 0;
	sigma = 1;

	if (call == 1)
	{
	  call = !call;
	  return (mu + sigma * (double) X2);
	}

	do
	{
	  U1 = -1 + ((double) rand () / RAND_MAX) * 2;
	  U2 = -1 + ((double) rand () / RAND_MAX) * 2;
	  W = pow (U1, 2) + pow (U2, 2);
	}
	while (W >= 1 || W == 0);

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return (mu + sigma * (double) X1);
}

void insereRuido(double *q, double *h, double *r, double *c, double *t, double *cnt, double *ms, double *sdev, double* cnst, double *out){

	//Sensor reads
	int i;
	for(i=0; i<14; i++)
		ms[i] = ms[i] + gauss_rand()*sdev[i];

	out[0] = ms[3];						//#4 Reactor level
	out[1] = h[2] + gauss_rand()*0.01;			//Head after pump
	out[2] = h[3] + gauss_rand()*0.01;			//Head after pipe 1
	out[3] = h[4] + gauss_rand()*0.01;			//Head after valve
	out[4] = h[5] + gauss_rand()*0.01;			//Head after pipe 2
	out[5] = h[6] + gauss_rand()*0.001;			//Leak head
	out[6] = ms[1];						//#2 Feed flowrate
	out[7] = q[2] + gauss_rand()*0.002;			//Reactor outlet flowrate
	out[8] = ms[8];						//#9 Product flowrate
	out[9] = q[3] + gauss_rand()*0.002;			//Leak flowrate
	out[10] = ms[6];					//#7 Reactor temperature
	out[11] = ms[10];					//#11 Coolant inlet pressure
	out[12] = h[8] + gauss_rand()*0.01;			//Head after pipe 1 (cooling system)
	out[13] = h[9] + gauss_rand()*0.01;			//Head afer valve
	out[14] = h[10] + gauss_rand()*0.01;			//Head after jacket obstruction
	out[15] = h[11] + gauss_rand()*0.01;			//Head after pipe 2
	out[16] = h[12] + gauss_rand()*0.001;			//Leak to environment head
	out[17] = h[13] + gauss_rand()*0.001;			//Leak to tank head
	out[18] = ms[7];					//#8 Coolant flowrate
	out[19] = q[8] + gauss_rand()*0.002;			//Flow rate in pipe 2
	out[20] = q[7] + gauss_rand()*0.002;			//Leak to environment flowrate
	out[21] = q[6] + gauss_rand()*0.002;			//Leak to tank flowrate
	out[22] = ms[9];					//#10 Coolant inlet temperature
	out[23] = t[4] + gauss_rand();				//Jacket temperature
	out[24] = r[0] + gauss_rand()*0.001;			//Generation rate of product B
	out[25] = r[1] + gauss_rand()*0.001;			//Generation rate of product C
	out[26] = ms[2];					//#3 Feed temperature
	out[27] = ms[0];					//#1 Feed concentration
	out[28] = ms[4];					//#5 Product A concentration
	out[29] = ms[5];					//#6 Product B concentration
	out[30] = c[3] + gauss_rand()*0.015;			//Product C concentration
	out[31] = ms[11];					//#12 Level controller output
	out[32] = ms[12];					//#13 Coolant controller output
	out[33] = ms[13];					//#14 Coolant setpoint
	out[34] = cnst[0] + gauss_rand()*0.000005; 
	out[35] = cnst[1] + gauss_rand()*0.0005;
	out[36] = cnst[2] + gauss_rand()*0.000005;
	out[37] = cnst[3] + gauss_rand()*0.000005;
}
