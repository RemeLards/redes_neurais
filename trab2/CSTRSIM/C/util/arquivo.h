void printMsg(FILE *arq, int tp);
void leComentario(int n, char *str, FILE *arq);
fault *transformaLinha(fault *ft, char *str);
void initVets(int *fl, double *nv, double *it, double *tu, int n);
void transLinha(char *str, int *fl, double *nv, double *it, double *tau, int k);
int leArquivoFalhas(FILE *arq, unsigned int *sd, double *its, int *git, int *falhas, double *nv, double *inst, double *tau);
void initFiles(FILE *cstr34, FILE *cstr14, FILE *matlab34, FILE *matlab14);