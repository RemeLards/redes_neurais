
typedef struct l{
	int fcode;
	double newValue;
	double inst;
	double tau;
	struct l *next;
} fault;

fault *newNode(int c, double nv, double it, double t);
fault *addNode(fault *root, fault *new);
fault *addFault(fault *root, int c, double nv, double it, double t);
fault *delNode(fault *root, int c);
fault *delLast(fault *root);
void printFaultList(fault *lf);
