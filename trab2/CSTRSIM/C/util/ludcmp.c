//Biblioteca utilizada. Retirada de Numerical Recipes in C 2ed

#include <math.h>
#include "nrutil.h"
#define TINY 1.0e-30

//Realiza a decomposicao LU de uma matriz
void ludcmp(double **a, int n, int *indx, double *d)
/*Given a matrix a[1..n][1..n] , this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation effected by the partial
pivoting; d is output as ±1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.*/
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;	//vv stores the implicit scaling of each row.
	vv=dvector(1,n);
	*d=1.0;		//No row interchanges yet.
	for (i=1;i<=n;i++) {	//Loop over rows to get the implicit scaling information.
		big=0.0;
		for (j=1;j<=n;j++)
		if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		//No nonzero largest element.
		vv[i]=1.0/big;		//Save the scaling.
	}
	for (j=1;j<=n;j++) {	//This is the loop over columns of Crout’s method.
		for (i=1;i<j;i++) {	//This is equation (2.3.12) except for i = j.
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;		//Initialize for the search for largest pivot element.
		for (i=j;i<=n;i++) {	//This is i = j of equation (2.3.12) and i = j + 1 . . . N of equation (2.3.13).
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
			//Is the figure of merit for the pivot better than the best so far?
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {	//Do we need to interchange rows?
			for (k=1;k<=n;k++) {	//Yes, do so...
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);		//...and change the parity of d.
			vv[imax]=vv[j];		//Also interchange the scale factor.
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		/*If the pivot element is zero the matrix is singular (at least to the precision of the
		algorithm). For some applications on singular matrices, it is desirable to substitute
		TINY for zero.*/
		if (j != n) {		//Now, finally, divide by the pivot element.
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}		//Go back for the next column in the reduction.
	free_dvector(vv,1,n);
}

//Resolve um sistema linear com base nos resultados da funcao anterior
void lubksb(double **a, int n, int *indx, double b[])
/*Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp . indx[1..n] is input
as the permutation vector returned by ludcmp . b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a , n , and indx are not modified by this routine
and can be left in place for successive calls with different right-hand sides b . This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion.*/
{
	int i,ii=0,ip,j;
	double sum;
	for (i=1;i<=n;i++) {		//When ii is set to a positive value, it will become the
		ip=indx[i];			//index of the first nonvanishing element of b. We now
		sum=b[ip];			//do the forward substitution, equation (2.3.6). The
		b[ip]=b[i];			//only new wrinkle is to unscramble the permutation
		if (ii)				//as we go.
		    for (j=ii;j<=i-1;j++)
			sum -= a[i][j]*b[j];
		else if (sum) ii=i;		//A nonzero element was encountered, so from now on we
		b[i]=sum;			//will have to do the sums in the loop above.
	}
	for (i=n;i>=1;i--) {		//Now we do the backsubstitution, equation (2.3.7).
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];	//Store a component of the solution vector X.
	}				//All done!
}

//Realiza a multiplicacao de duas matrizes
void multmatrix(double **a, double **b, double **c, int n){
	int i, j, k;
	double sum;
	for(i=1; i<=n; i++){
		for(j=1; j<=n; j++){
			sum = 0;
			for(k=1;k<=n;k++){
				sum += a[i][k]*b[k][j];
				c[i][j] = sum;
			}
		}
	}
}

//Realiza a multiplicacao de uma matriz por um vetor
void matrixvector(double **a, double *b, int n, double *r){
		int i, j, k;
		double sum;
		for(i=1; i<=n; i++){
			sum = 0;
			for(j=1; j<=n; j++){
				sum += a[i][j]*b[j];
				r[i] = sum;
			}
		}
}

//Realiza a subtracao de dois vetores
void subvector(double *a, double *b, int n, double *r){
	int i;
	for(i=1; i<=n; i++){
		r[i] = a[i]-b[i];
	}
}

//Inicializa uma matriz com zero
void initmatrix(double **a, int n){
	int i, j;
	for(i=1; i<=n; i++){
		for(j=1; j<=n; j++){
			a[i][j] = 0;
		}
	}
}
