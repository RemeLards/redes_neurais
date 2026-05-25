//CSTR parameters

//Gravity acceleration
static double g = 9.81;

//Gas constant
static double R = 8.3144621;

//Heat exchange constant
static double ua = 1901.0;

//Floor area of CSTR
static double tarea = 1.5;

//Pipe cross-section area
static double A = 0.011;

//Cooling system cross-section area
static double Ac = 0.8;

//Primary heat of reaction
static double dhb = 30000.0;

//Secondary heat of reaction
static double dhc = -10000.0;

//Primary rate constant pre-exponential factor
static double k0b = 2500.0;

//Secondary rate constant pre-exponential factor
static double k0c = 3000.0;

//Controller parameters
static double Kp = 0.3;
static double Ti = 0.1;
static double Td = 0.15;

//Controller parameters
static double Kp_master = 0.5;
static double Ti_master = 2.0;
static double Td_master = 0.25;

//Controller parameters
static double Kp_slave = 0.15;
static double Ti_slave = 0.01;
static double Td_slave = 0.035;