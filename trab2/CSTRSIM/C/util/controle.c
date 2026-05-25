//Controle PID
#include "params.h"

double ctrl(double e, double e_1, double e_2, double pre_out, double dt){

	double t1, t2, t3, out;
	
	t1 = (1 + (dt/Ti) + (Td/dt)) * e;
	t2 = (-1 - (2 * Td/dt)) * e_1;
	t3 = (Td/dt) * e_2;
	out =  pre_out - (Kp * (t1 + t2 + t3));
	
	if(out > 1.0) return 1.0;
	if(out <= 0.0) return 0.000001;
	
	return out;
}

double master_ctrl(double e, double e_1, double e_2, double pre_out, double dt){

	double t1, t2, t3, out;
	
	t1 = (1 + (dt/Ti_master) + (Td_master/dt)) * e;
	t2 = (-1 - (2 * Td_master/dt)) * e_1;
	t3 = (Td_master/dt) * e_2;
	out =  pre_out - (Kp_master * (t1 + t2 + t3));
	
	if(out > 1.0) return 1.0;
	if(out <= 0.0) return 0.0;
	
	return out;
}

double slave_ctrl(double e, double e_1, double e_2, double pre_out, double dt){

	double t1, t2, t3, out;
	
	t1 = (1 + (dt/Ti_slave) + (Td_slave/dt)) * e;
	t2 = (-1 - (2 * Td_slave/dt)) * e_1;
	t3 = (Td_slave/dt) * e_2;
	out =  pre_out - (Kp_slave * (t1 + t2 + t3));
	
	if(out > 1.0) return 1.0;
	if(out <= 0.0) return 0.000001;
	
	return out;
}