//CSTR Simulator

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util/params.h"
#include "util/gaussrand.h"
#include "util/nrutil.h"
#include "util/ludcmp.h"
#include "util/newtrap.h"
#include "util/controle.h"
#include "util/lista.h"
#include "util/arquivo.h"
#include "util/faults.h"

#define BOOLSTR(b)      (b ? "TRUE " : "FALSE" )

int procQuim(int argc, char *argv[]){
	
	//Output files
	FILE *cstr38, *cstr18;
	FILE *matlab34, *matlab14;
		
	cstr38 = fopen("output/cstr38.txt", "w");
	if(!cstr38){
		fprintf(stderr, "ERROR: could not open output file\n");
		return -1;
	}
	printMsg(cstr38,1);
	
	cstr18 = fopen("output/cstr18.txt", "w");
	if(!cstr18){
		fprintf(stderr, "ERROR: could not open output file\n");
		return -1;
	}
	printMsg(cstr18,2);

	matlab34 = fopen("matlab/cstr38.txt", "w");
	if(!matlab34){
		fprintf(stderr, "ERROR: could not open output file\n");
		return -1;
	}
	
	matlab14 = fopen("matlab/cstr18.txt", "w");
	if(!matlab14){
		fprintf(stderr, "ERROR: could not open output file\n");
		return -1;
	}

	FILE *conf = NULL;
	if(argc < 2){
		conf = fopen("default.cfg", "r");
		if(!conf){
			fprintf(stderr, "ERROR: could not open fault configuration file\n");
			return -1;
		}
		fprintf(stdout, "File: default.cfg\n");
	}else{
		conf = fopen(argv[1], "r");
		if(!conf){
			fprintf(stderr, "ERROR: could not open fault configuration file\n");
			return -1;
		}
		fprintf(stdout, "File: %s\n", argv[1]);
	}
	
	int i, s0, s1, nf = 51;
	int it, ij, git, fl[nf];
	double simTime, maxTime, inst[nf];
	double nv[nf], tau[nf], oldVal[nf], newVal[nf];
	unsigned int seed = 0;

	initVets(fl, nv, inst, tau, nf); // Set all 'fl, 'nv, 'inst', 'tau', 'nf' arrays to zero
	leArquivoFalhas(conf, &seed, &maxTime, &git, fl, nv, inst, tau);

	printf("Iterations: %f\n", maxTime);
	for(i=1;i<nf;i++)
		if( fl[i] != 0 )
	 		printf("Total # faults=%2d fault# =%2d extent=%7.3lf instant=%6.3f tau=%7.3lf\n",
				nf-1,i,nv[i],inst[i],tau[i]);	/**/
	
	simTime = 0.0;
	it = 0;
	s0 = 0;
	s1 = 0;

	//Seed for random number generator
	if(seed != -1){
		srand(seed);
		printf("Seed for random number generator: %d\n", seed);
	}
	else{
		seed = (unsigned)time(NULL);
		srand(seed);
		printf("Seed for random number generator: %d\n", seed);
	}
	
	//Initial root vector for Newton-Raphson
	double h[14], q[9], out[38], ms[14];
	
	//Sensor measurements
	ms[0] = 20.0;		//#1 Feed concentration
	ms[1] = 0.25;		//#2 Feed flowrate
	ms[2] = 30.0;		//#3 Feed temperature
	ms[3] = 2.0;		//#4 Reactor level
	ms[4] = 2.85;		//#5 Product A concentration
	ms[5] = 17.11;		//#6 Product B concentration
	ms[6] = 80.0;		//#7 Reactor temperature
	ms[7] = 0.9;		//#8 Coolant flowrate
	ms[8] = 0.25;		//#9 Product flowrate
	ms[9] = 20.0;		//#10 Coolant inlet temperature
	ms[10] = 10.0;		//#11 Coolant inlet pressure
	ms[11] = 0.04;		//#12 Level controller output
	ms[12] = 0.60;		//#13 Coolant controller output
	ms[13] = 0.9;		//#14 Coolant setpoint

	q[0] = 0;
	//Feed flowrate
	q[1] = 0.25;

	//Reactor
	h[0] = 47.0;	//Pressure head of the pump
	h[1] = 2.0;		//Reactor head (level)
	h[2] = 5.0;		//Head after pump
	h[3] = 2.0;		//Head after pipe 1
	h[4] = 1.0;		//Head after valve
	h[5] = 1.0;		//Head after pipe 2
	h[6] = 0.0;		//Leak head
	q[2] = 0.25;	//Flow in pipe 1
	q[3] = 0.25;	//Flow in pipe 2
	q[4] = 0.0;		//Leak flow

	//Cooling system
	h[7] = 10.0;	//Pressure head of the cooling system
	h[8] = 5.0;		//Head after pipe 1
	h[9] = 2.0;		//Head after valve
	h[10] = 2.0;	//Head after jacket obstruction
	h[11] = 2.0;	//Head after pipe 2
	h[12] = 0.0;	//Leak to environmet head
	h[13] = 0.0;	//Leak to reactor head
	q[5] = 0.9;		//Flow in pipe 1
	q[6] = 0.0;		//Leak to reactor flow
	q[7] = 0.0;		//Leak to environment flow
	q[8] = 0.9;		//Flow in pipe 2

	double sdev[14] = {0.15,0.0002,0.15,0.01,0.02,0.14,0.15,0.0003,0.000002,0.15,0.01,0.00005,0.00005,0.00005};

	//Additional variables for effluent faults
	double h5fault = h[5], h11fault = h[11];
	
	//Total number os variables (flows + heads). Except reactor inlet head,
	//the reactor head and cooling system intake head
	int n = 18;
	
	//Temperature vector
	double t[5], c[4], r[2], cc0, dn, dt, ktb1;
	double v, v0, vold, rho, rhold, cpc, cpcold;
	double eb, ec, qrxn, qext, qc;
	double flowinoutdiff, sumconc, mol, molt0, molin, molout, dmolin, dmolout, dmolout_OLD;
	double k[13], cnst[4];
	double k0[5], t2_old, h1_old;

	molt0 = 60;	// moles at t=0, needed for constraint 2 (Mol Balance)
	//CSTR volume
	v = v0 = 3.0;

	t[0] = 0;
	//Initial temperatures
	t[1] = 30.0;		//Feed temperature
	t[2] = 80.0;		//Reactor temperature
	t[3] = 20.0;		//Coolant inlet temperature
	t[4] = 40.0;		//Jacket temperature
	
	//Constants
	ua = 1901.0;
	rho = 1000.0;
	cpc = 4.2;
	tarea = 1.5;
	ktb1 = 10.0;
	dt = 0.02;
	
	//Frequency factor for Arrhenius equation
	k0b = 2500.0;
	k0c = 3000.0;
	
	//Activation energy
	eb = 25000.0;
	ec = 45000.0;
	
	//Leaks and obstruction coefficients
	double mvo[4];
	
	mvo[0] = 1.0e-7;	// ==> k[2]=1e-14	Leak in CSTR exit
	mvo[1] = 1.0e-7;	// ==> k[8]=1e-14	Leak from jacket to environment
	mvo[2] = 1.0e-7;	// ==> k[6]=1e-14	Leak from jacket to reactor
	mvo[3] = 1.0;		// ==> k[9]=1		Jacket obstruction

	//Constraints variables
	flowinoutdiff = 0.0;
	molin = 0.0;
	molout = 0.0;

	//Nominal values
	k0[0] = 10.0;	// k[1] nominal value
	k0[1] = 250;	// k[4] nominal value
	k0[2] = 5.18;	// k[5] nominal value
	k0[3] = 0.0;	// k[9] nominal value
	k0[4] = 4.22;	// k[10] nominal value

	//Controller variables
	double sp[3], cnt[4], cntold[4], error[3][2], sp_error[3];
	
	sp[0] = 0;
	sp[1] = 2.0;
	sp[2] = 80.0;
	sp[3] = 0.9;

	cnt[0] = 0;
	cnt[1] = 0.101;
	cnt[2] = 0.9;
	cnt[3] = 0.60;

	cntold[0] = 0;
	cntold[1] = 0.101;
	cntold[2] = 0.9;
	cntold[3] = 0.60;

	error[0][0] = 0.0;
	error[0][1] = 0.0;
	
	error[1][0] = 0.0;
	error[1][1] = 0.0;
	
	error[2][0] = 0.0;
	error[2][1] = 0.0;
	
	//Minor loss coefficients
	k[0] = 0;			        // dummy -- not used
	k[1] = 10.0;			    // Pipe:  Reactor exit (fixed)
	k[2] = 1/(mvo[0]*mvo[0]);	// Leak:  Pump (depends on leak 'valve')
	k[3] = 1/(cnt[1]*cnt[1]);	// Valve: Product flow
	k[4] = 250.0;		    	// Pipe:  Effluent pipe (fixed)
	k[5] = 5.18;		    	// Pipe:  Cooling inlet (fixed)
	k[6] = 1/(cnt[3]*cnt[3]);	// Valve: Cooling flow
	k[7] = 1/(mvo[2]*mvo[2]);	// Leak: Jacket to tank
	k[8] = 1/(mvo[1]*mvo[1]);	// Leak: Jacket to environment
	k[9] = (mvo[3] >= 1)? 0.0 : (1/(mvo[3]*mvo[3]));	// Obstruction: Intra-jacket
	k[10] = 4.22;			    // Pipe: Cooling outlet
	k[11] = 1.0;	    		// Pump: Flow coefficient (fixed)
	k[12] = 2.0;		    	// Pump: Major loss coefficient

	//Heat of reaction
	dhb = 30000.0;
	dhc = -10000.0;
	
	qrxn = 0.0;
	qext = 0.0;
	qc = 0.0;

	//Initial values for concentrations
	c[0] = 20.0;		//Feed concentration
	c[1] = 2.850;		//Product A concentration
	c[2] = 17.114;		//Product B concentration
	c[3] = 0.0226;		//Product C concentration
	cc0 = c[3];

	//Stores initial values
	oldVal[2] = ktb1;	oldVal[3] = mvo[3];	oldVal[4] = mvo[1];	
	oldVal[5] = mvo[2];	oldVal[6] = mvo[0];	oldVal[7] = h[0];
	oldVal[8] = ua; 	oldVal[9] = qext; 	oldVal[10] = eb;
	oldVal[11] = ec; 	oldVal[12] = q[1]; 	oldVal[13] = t[1];
	oldVal[14] = c[0];	oldVal[15] = t[3]; 	oldVal[16] = h[7]; 	
	oldVal[17] = h[11]; 	oldVal[18] = 26.0;	oldVal[19] = sp[1];
	oldVal[20] = sp[2];	oldVal[21] = cnt[1]; 	oldVal[22] = cnt[3];

	//Main loop
	printf("Running simulation...\n\n");
	while(simTime < maxTime) {

		printf("Iteration=%8d -- Time [min] %8.2f of total %8.2f\r", it, simTime, maxTime );
		
		//Handles all faults
		//=============================================
		initSensorFaults(h, q, t, c, cnt, oldVal);
		faultHandler(simTime, fl, inst, nv, tau, oldVal, ms, newVal); // eventually change oldVal to newVal, if fault
		//=============================================

		ktb1 = newVal[2];		mvo[3] = newVal[3];		mvo[1] = newVal[4];
		mvo[2] = newVal[5];		mvo[0] = newVal[6];		h[0] = newVal[7];
		ua = newVal[8]; 		qext = newVal[9];		eb = newVal[10];
		ec = newVal[11]; 		q[1] = newVal[12]; 		t[1] = newVal[13];
		c[0] = newVal[14]; 		t[3] = newVal[15]; 		h[7] = newVal[16];
		h11fault = newVal[17];		h5fault = newVal[18];		sp[1] = newVal[19];
		sp[2] = newVal[20];		cnt[1] = newVal[21];		cnt[3] = newVal[22];

		//PID controller call (reactor outlet valve)
		if(fl[21] != 1){
			sp_error[0] = sp[1] - ms[3];
			cnt[1] = ctrl(sp_error[0], error[0][0], error[0][1], cntold[1], dt);
			error[0][0] = sp_error[0];
			error[0][1] = error[0][0];
			cntold[1] = cnt[1];
		}
		
		//Cascade PID controller call (cooling system valve)
		sp_error[1] = sp[2] - ms[6];
		cnt[2] = master_ctrl(sp_error[1], error[1][0], error[1][1], cntold[2], dt);
		error[1][0] = sp_error[1];
		error[1][1] = error[1][0];
		cntold[2] = cnt[2];
		
		//Cascade PID controller call (cooling system valve)
		if(fl[22] != 1){
			sp_error[2] = ms[7] - sp[3];
			cnt[3] = slave_ctrl(sp_error[2], error[2][0], error[2][1], cntold[3], dt);
			error[2][0] = sp_error[2];
			error[2][1] = error[2][0];
			cntold[3] = cnt[3];
			sp[3] = cnt[2];
		}
		
		// Peep some values
		//fprintf(stdout,"PEEP BEFORE> h[7]=%10.4e h[8]=%10.4e h[9]=%10.4e h[10]=%10.4e h[11]=%10.4e h[10]-h[11]=%e\n", h[7],h[8],h[9],h[10],h[11],h[10]-h[11]);
		// fprintf(stdout,"PEEP BEFORE> cA=%10.4e cB=%10.4e cC=%10.4e cC0=%10.4e\n", ms[4],ms[5],ms[8],cc0);

		//Newton-Raphson
		//=============================================
		newtRap(q, h, k, cnt[1], mvo, ktb1, cnt[3], n);
		//=============================================
		//Check for fault in effluent
		if(fl[17] == 1) h[11] = h11fault;
		if(fl[18] == 1) h[5] = h5fault;

		// fprintf(stdout,"PEEP AFTER > h[7]=%10.4e h[8]=%10.4e h[9]=%10.4e h[10]=%10.4e h[11]=%10.4e h[10]-h[11]=%e\n", h[7],h[8],h[9],h[10],h[11],h[10]-h[11]); getchar();
		// fprintf(stdout," PEEP AFTER> cA=%10.4e cB=%10.4e cC=%10.4e cC0=%10.4e\n", ms[4],ms[5],ms[8],cc0);

		//Cooling system temperature
		t[4] = (ua*t[2] + rho*cpc*q[8]*t[3])/(rho*cpc*q[8] + ua);
		
		//Heat update
		qc = ua*(t[2]-t[4]);
		
		//Reactor volume
		//Head h[1]
		vold = v;
		v = vold + dt*(q[1] + q[6] - q[2]);
		h1_old = h[1];
		h[1] = v/tarea;
		
		//p(rho)
		rhold = rho;
		rho = (1/v)*(vold*rho)+(1/v)*(dt*(q[1]*rho+q[7]*rho-q[2]*rho));
		
		//cpc
		cpcold = cpc;
		cpc = (1/v)*(vold*cpc)+(1/v)*(dt*(q[1]*cpc+q[7]*cpc-q[2]*cpc));
		
		//Generation rates of product B and C
		r[0] = c[1]*k0b*exp(-eb/(R*(273.5+t[2])));
		r[1] = c[1]*k0c*exp(-ec/(R*(273.5+t[2])));
		
		//Euler
		c[1] = (vold/v)*c[1] + (dt/v)*(c[0]*q[1]-c[1]*q[2]-(r[0]+r[1])*vold);
		c[2] = (vold/v)*c[2] + (dt/v)*(-c[2]*q[2]+r[0]*vold);
		c[3] = (vold/v)*c[3] + (dt/v)*(-c[3]*q[2]+r[1]*vold);
		//fprintf(stdout, "R - %lf %lf %lf\n", c[1],c[2],c[3]);
		
		//Heat update
		qrxn = (dhb*r[0]+dhc*r[1])*vold;
		dn = v*rho*cpc;
		
		//Euler
		t2_old = t[2];
		t[2] = (vold/v)*t[2] +(dt/vold)*((1/(rho*cpc))*(qrxn+qext-qc)+q[1]*t[1]+q[6]*t[4]-q[2]*t[2]);
		//fprintf(stdout, "%lf\n", t[2]);
		
		//Evaluate safety systems
		if((h1_old >= 2.75 || t2_old >= 130.0) && s0 == 0){
			oldVal[12] = 0.0;
			printf("\nEmergency shut down initiated at iteration %d at time %f...\n", it, simTime);
			printf("Reason h1=%f >= 2.75 = %s --- T2=%.2f > 130=%s\n",
				h1_old,BOOLSTR(h1_old >= 2.75),t2_old,BOOLSTR(t2_old >= 130.0));
			s0 = 1;
		}
		if(h1_old <= 1.2 && s1 == 0){
			oldVal[7] = 0.0;
			printf("\nLow level forces pump shut down with h1=%f at iteration %d at time %f...\n", h1_old, it, simTime);
			s1 = 1;
		}
		if(h1_old <= 0.0){
			printf("\nFailure due to low level with h1=%f at iteration %d at time %f...\n\n", h1_old, it, simTime);
			return 1;
		}

		//Gaussian noise
		insereRuido(q, h, r, c, t, cnt, ms, sdev, cnst, out);


		// Noise must be added before constraint test

		// CONSTRAINTS
		// =====================================================
		// ### INVENTORY ###
		flowinoutdiff = flowinoutdiff + (q[1] - q[4])*dt; // <<<<<<<<<<<< OLD : Q_1 - Q_4
        // ==========================================================
        // TO DO: Correction for the C implementation not yet working
		// flowinoutdiff = flowinoutdiff + (q[1] - q[4] + q[6] - q[3])*dt; // <<<<<<<<<<<< CORRECTION : Q_1 - Q_4 + Q_6 - Q_3
        // ==========================================================
		cnst[0] = v - v0 - flowinoutdiff;

		// ### MOL BALANCE ###
		sumconc = c[1] + c[2] + cc0;
		mol = sumconc*v;
		// // fprintf(stdout,"sumconc=%.10e V=%.10e mol=%.10e\n", sumconc, v, mol);
		// // getchar();

		dmolin = c[0]*q[1]*dt;
		molin = molin + dmolin;
		dmolout = sumconc*ms[8]*dt; // <<<<<<<<<<<< OLD : Q_4 = ms[8]
        // ==========================================================
        // TO DO: Correction for the C implementation not yet working
		//dmolout = sumconc*(ms[8]+q[3])*dt; // <<<<<<<<<<<< CORRECTION : Q_4 + Q_3
        // ==========================================================
		molout = molout + dmolout;
		cnst[1] = mol - molt0 - molin + molout;	//Mol balance, cA(t=0)*V(t=0)=20*3=60

		// double z2old, auxdiff;
		// auxdiff = z2old - cnst[1]; z2old = cnst[1];
		// fprintf(stdout,"z2=%.10e dmolout_OLD=%.10e dmolout=%.10e auxdiff=%.10e\n",
		// 		cnst[1], dmolout_OLD, dmolout, auxdiff);
		// fprintf(stdout,"z2=%.10e sumc=%50.42e cA0=%.2e cA=%.6e cB=%.6e cC=%.6e cCt0=%.2e Q1=%.2e Q4=%.2e(meas=%.2e) V=%.2e mol=%.20e molin=%.5e molout=%.5e dmol=%.20e diffio=%.8e\n",
		// 	cnst[1], sumconc, c[0], c[1], c[2], c[3], cc0, q[1], q[4], ms[8], v, mol, molin, molout, dmolin-dmolout, molin-molout ); /**/
		// getchar();

		// fprintf(stdout,"z2=%.10e sumc=%50.42e Q1-Q4=%.20e Q1=%.20e Q2=%.20e Q3=%.20e Q4=%.20e Q6=%.20e dmol=%.20e diffio=%.8e\n",
		// 	cnst[1], sumconc, q[1]-q[4], q[1], q[2], q[3], q[4], q[6], dmolin-dmolout, molin-molout ); /**/

		// ### COOLING CIRCUIT HEAD LOSS ###
		cnst[2] = (h[7] - h[11]) - (k[6] + k0[2] + k0[3] + k0[4])*q[5]*q[5];
		/* fprintf(stdout,"z3=%.30e h7=%e h11=%e Q5=%e k6=%e kN2=%e kN3=%e kN4=%e\n", 
			cnst[2], h[7], h[11], q[5], k[6], k0[2], k0[3], k0[4] );    /**/

		// ### EFFLUENT HEAD LOSS ###
		cnst[3] = (h[1] + h[0] - h[5]) - (k[11]*q[4]) -(k[12] + k[3] + k0[0] + k0[1])*q[4]*q[4];
		/* fprintf(stdout,"z4=%e h1=%e h0=%e h5=%e  Q4=%e k11=%e k12=%e k3=%e kN0=%e kN1=%e\n", 
			cnst[3], h[1], h[0], h[5], q[4], k[11], k[12], k[3], k0[0], k0[1] );    /**/
			
		//printf("Loss factors: K1=%.2e K2=%.2e K3=%.2e K4=%.2e K5=%.2e K6=%.2e K7=%.2e K8=%.2e K9=%.2e K10=%.2e K11=%.2e K12=%.2e\n",
		//      k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[9],k[10],k[11],k[12]);
			
		//Print output
		if(( (it+1) % (int)(git/dt) ) == 0){

			//Print without gaussian noise
			/*fprintf(cstr38, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
				    h[1], h[2], h[3], h[4], h[5], h[6], q[1], q[2], q[4], q[3], t[2]);
			fprintf(matlab34, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
				    h[1], h[2], h[3], h[4], h[5], h[6], q[1], q[2], q[4], q[3], t[2]);
			
			fprintf(cstr38, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
				    h[7], h[8], h[9], h[10], h[11], h[12], h[13], q[5], q[8], q[7], q[6], t[4]);
			fprintf(matlab34, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
				    h[7], h[8], h[9], h[10], h[11], h[12], h[13], q[5], q[8], q[7], q[6], t[4]);

			fprintf(cstr38, "%lf %lf %lf %lf %lf %lf %lf %lf ", t[3], r[0], r[1], t[1], c[0], c[1], c[2], c[3]);
			fprintf(matlab34, "%lf %lf %lf %lf %lf %lf %lf %lf ", t[3], r[0], r[1], t[1], c[0], c[1], c[2], c[3]);

			fprintf(cstr38, "%lf %lf %lf ", cnt[1], cnt[3], cnt[2], cnst[0], cnst[1], cnst[2], cnst[3]);
			fprintf(matlab34, "%lf %lf %lf ", cnt[1], cnt[3], cnt[2], cnst[0], cnst[1], cnst[2], cnst[3]);
			
			fprintf(cstr18, "%lf %lf %lf %lf %lf %lf %lf ", c[0],q[1],t[1],h[1],c[1],c[2],t[2]);
			fprintf(cstr18, "%lf %lf %lf %lf %lf %lf %lf ", q[5],q[4],t[3],h[7],cnt[1],cnt[3],cnt[2]);
			fprintf(cstr18, "%lf %lf %lf %lf ", cnst[0], cnst[1], cnst[2], cnst[3]);
			
			fprintf(matlab14, "%lf %lf %lf %lf %lf %lf %lf ", c[0],q[1],t[1],h[1],c[1],c[2],t[2]);
			fprintf(matlab14, "%lf %lf %lf %lf %lf %lf %lf ", q[5],q[4],t[3],h[7],cnt[1],cnt[3],cnt[2]);
			fprintf(matlab14, "%lf %lf %lf %lf ", cnst[0], cnst[1], cnst[2], cnst[3]);
			/******************************************************************************************************/

			//Print with gaussian noise
			for(ij=0; ij<38; ij++)
				fprintf(cstr38, "%lf ", out[ij]);

			for(ij=0; ij<38; ij++)
				fprintf(matlab34, "%lf ", out[ij]);

			fprintf(cstr18, "%lf %lf %lf %lf %lf %lf %lf ", out[27],out[6],out[26],out[0],out[28],out[29],out[10]);
			fprintf(cstr18, "%lf %lf %lf %lf %lf %lf %lf ", out[18],out[8],out[22],out[11],out[31],out[32],out[33]);
			fprintf(cstr18, "%lf %lf %lf %lf ", out[34], out[35], out[36], out[37]);

			fprintf(matlab14, "%lf %lf %lf %lf %lf %lf %lf ", out[27],out[6],out[26],out[0],out[28],out[29],out[10]);
			fprintf(matlab14, "%lf %lf %lf %lf %lf %lf %lf ", out[18],out[8],out[22],out[11],out[31],out[32],out[33]);
			fprintf(matlab14, "%lf %lf %lf %lf ", out[34], out[35], out[36], out[37]);
			/**********************************************************************************************************/
			//Print fault code
			int hasfl=0;
			for(ij=2; ij<nf; ij++){
				if(fl[ij] == 1){
					hasfl = 1;
					fprintf(cstr38, "%d_", ij);
					fprintf(cstr18, "%d_", ij);
					fprintf(matlab34, "%d_", ij);
					fprintf(matlab14, "%d_", ij);
				}
			}

			if(hasfl == 0){
				fprintf(cstr38, "1");
				fprintf(cstr18, "1");
				fprintf(matlab34, "1");
				fprintf(matlab14, "1");
			}

			fprintf(cstr38, "\n");
			fprintf(cstr18, "\n");
			fprintf(matlab34, "\n");
			fprintf(matlab14, "\n");
		}

		it = it+1;
		simTime = simTime + dt;
	}
	
	printf("\nSimulation complete, check the output folder(s)\n\n");
	
	fclose(cstr38);
	fclose(cstr18);
	fclose(matlab34);
	fclose(matlab14);
		
	return 0;
}

