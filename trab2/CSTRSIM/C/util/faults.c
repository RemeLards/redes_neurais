//Manipula as falhas que foram passadas pelo arquivo
#include <math.h>
#include <stdio.h>

void initSensorFaults(double *h, double *q, double *t, double *c, double *cnt, double *oldVal){
	
	//Sensor faults - fixed bias
	oldVal[23] = c[0];			//#1 Feed concentration
	oldVal[24] = q[1];			//#2 Feed flowrate
	oldVal[25] = t[1];			//#3 Feed temperature
	oldVal[26] = h[1];			//#4 Reactor level
	oldVal[27] = c[1];			//#5 Product A concentration
	oldVal[28] = c[2];			//#6 Product B concentration
	oldVal[29] = t[2];			//#7 Reactor temperature
	oldVal[30] = q[5];			//#8 Coolant flowrate
	oldVal[31] = q[4];			//#9 Product flowrate
	oldVal[32] = t[3];			//#10 Coolant inlet temperature
	oldVal[33] = h[7];			//#11 Coolant inlet pressure
	oldVal[34] = cnt[1];		//#12 Level controller output
	oldVal[35] = cnt[3];		//#13 Coolant controller output
	oldVal[36] = cnt[2];		//#14 Coolant setpoint

	//Sensor faults - fixed value
	oldVal[37] = c[0];			//#1 Feed concentration
	oldVal[38] = q[1];			//#2 Feed flowrate
	oldVal[39] = t[1];			//#3 Feed temperature
	oldVal[40] = h[1];			//#4 Reactor level
	oldVal[41] = c[1];			//#5 Product A concentration
	oldVal[42] = c[2];			//#6 Product B concentration
	oldVal[43] = t[2];			//#7 Reactor temperature
	oldVal[44] = q[5];			//#8 Coolant flowrate
	oldVal[45] = q[4];			//#9 Product flowrate
	oldVal[46] = t[3];			//#10 Coolant inlet temperature
	oldVal[47] = h[7];			//#11 Coolant inlet pressure
	oldVal[48] = cnt[1];		//#12 Level controller output
	oldVal[49] = cnt[3];		//#13 Coolant controller output
	oldVal[50] = cnt[2];		//#14 Coolant setpoint
}

/*
	simTime = simulation time
	fl = fault flag, 1=present 0=absent
	inst = starting instant
	nv = new value
	tau = time constant
	oldVal = nominal value
	ms = measured value (if applicable)
	newVal = new value

	FAULT RANGES (add one in the tables in the paper):
	0	Normal
	1-21	Process Faults (Table 4 in paper)
	22-35	Sensor Faults with fixed bias (Table 5, first part, in paper)
	36-50	Sensor Faults with fixed value (Table 5, second part, in paper)
*/
void faultHandler(double simTime, int *fl, double *inst, double *nv, double *tau, double *oldVal, double *ms, double *newVal)
{
	//Faults
	int i;
	for(i=2; i <= 22; i++){
		if(simTime >= inst[i] && inst[i] != 0){
			fl[i] = 1;
			//printf("faultHandler> simTime=%8.2f i=%2d inst=%8.2f fl=%3d tau=%6.2lf\n",simTime,i,inst[i],fl[i],tau[i]);
			newVal[i] = nv[i]-(nv[i]-oldVal[i])*exp(tau[i]*(inst[i]-simTime));
		}else{
			fl[i] = -1;
			newVal[i] = oldVal[i];
		}
	}

	for(i=2; i<= 15; i++){
		if(simTime >= inst[i+21] && inst[i+21] != 0){
			fl[i+21] = 1;
			//printf("faultHandler> simTime=%8.2f i=%2d inst=%8.2f fl=%3d tau=%6.2lf\n",simTime,i,inst[i],fl[i],tau[i]);
			ms[i-2] = nv[i+21]-(nv[i+21]-0.0)*exp(tau[i+21]*(inst[i+21]-simTime)) + oldVal[i+21];
		}else if(simTime >= inst[i+35] && inst[i+35] != 0){
			fl[i+35] = 1;
			//printf("faultHandler> simTime=%8.2f i=%2d inst=%8.2f fl=%3d tau=%6.2lf\n",simTime,i,inst[i],fl[i],tau[i]);
			ms[i-2] = nv[i+35]-(nv[i+35]-oldVal[i+35])*exp(tau[i+35]*(inst[i+35]-simTime));
		}else{
			fl[i+21] = -1;
			fl[i+35] = -1;
			ms[i-2] = oldVal[i+35];
		}
	}
	for(i=0; i < 51; i++)	{
		if( fl[i]==1 )	// fault active
		{
			// printf("faultHandler> simTime=%8.2f i=%2d inst=%8.2f fl=%3d tau=%6.2lf oldVal=%6.2lf newVal=%6.2lf",
			//	simTime,i,inst[i],fl[i],tau[i],oldVal[i],newVal[i]);
			//if( i < 15 )	// Measured variables
			//	printf("  ms[%3d]=%6.2lf", i, ms[i] );
			//printf("\n");
		}
	}	/**/
}
