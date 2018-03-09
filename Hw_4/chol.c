#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *load_matrix(char *filename, int *n, int *m); // carga una matriz
void *banachiewicz(float *matrix_A,float *L,int rows); //funcion para descomposicion cholesky-banachiewicz 
void *diag(int j,float *L,float *matrix_A,int rows); //funcion auxiliar 1 de banachiewicz para los elementos diagonales
void *nondiag(int i,int j,float *L,float *matrix_A,int rows);//funcion auxiliar 2 de banachiewicz para los elementos no diagonales
void *transponer(float *L,float *T, int rows); // transpone una matriz en otra
void *iniciocero(float *L, int N); // llena una matriz de ceros(no la inicializa ni le aloca memoria)
void multimatri(float *matrix_A,float *matrix_B,int n_row,int n_cols); // multiplica dos matrices
void *soltrianinf(float *matrix_A,float *matrix_B,float *matrix_X,int n_rows,int n_cols);// resuelve una matriz inferior
void *soltriansup(float *matrix_A,float *matrix_B,float *matrix_X,int n_rows,int n_cols);// resuelve una matriz superior 

int main(int argc, char **argv){

float *matrix_b,*matrix_X,*matrix_Y,*matrix_A, *L, *T;
int rows,cols,N,j,i;
matrix_b = load_matrix(argv[2], &rows, &cols);
matrix_A = load_matrix(argv[1], &rows, &cols);
N = rows*cols;
L = malloc(N * sizeof(float));
T = malloc(N * sizeof(float));
matrix_X = malloc(rows * sizeof(float));
matrix_Y = malloc(rows * sizeof(float));


banachiewicz(matrix_A,L,rows);
transponer(L,T,rows);
soltrianinf(L,matrix_b,matrix_Y,rows,cols);
soltriansup(T,matrix_Y,matrix_X,rows,cols);

    for(j=0;j<rows;j++){
      printf("%f\n", matrix_X[j]);
    }
    printf("\n");
 

/*iniciocero(L,N);
banachiewicz(matrix_A,L,rows);
transponer(L,T,rows);

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      printf(" %f ", L[i*cols + j]);
    }
    printf("\n");
  }


for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      printf(" %f ", T[i*cols + j]);
    }
    printf("\n");
  }
*/
}

//llena una matriz de elementos nulos
void *iniciocero(float *L, int N){

int i;
	for(i = 0;i < N; i++){
	L[i] = 0.0;
	}
}

// genera la traspuesta de una matriz en otra matriz
void *transponer(float *L,float *T, int rows){

int i,j;
	for(i = 0;i < rows; i++){
		for( j = 0; j < rows; j++){
		T[i*rows + j] = L[j*rows + i];
		}
	}
}

//hace la descomposicion cholesky-banachiewicz
void *banachiewicz(float *matrix_A,float *L,int rows){

int i,j;

	for(j = 0;j < rows;j++){
		for(i = j;i < rows;i++){
			if(i==j){
			diag(j,L,matrix_A,rows);
			}
			else{
			nondiag(i,j,L,matrix_A,rows);
			}
		}
	}
} 

//calcula un elemento diagonal de la descomposicion cholesky-banachiewicz
void *diag(int j,float *L,float *matrix_A,int rows){

int k;
float l=0.0;

		for(k = 0;k < j;k++){
		l = l + pow(L[j*rows + k],2);
		}

L[j*rows + j] = sqrt(matrix_A[j*rows + j] - l);
	
}

//calcula un elemento no diagonal de la descomposicion cholesky-banachiewicz
void *nondiag(int i,int j,float *L,float *matrix_A,int rows){

int k;
float l=0.0;

	for(k = 0;k < j;k++){
	l = l + (L[i*rows + k])*(L[j*rows + k]);
	}
L[i*rows + j] = (1.0/L[j*rows + j])*(matrix_A[i*rows + j] - l);
}

//tomado del repositorio de Jaime Forero
float *load_matrix(char *filename, int *n, int *m){
  float *matrix;
  FILE *in;
  int n_row, n_cols;
  int i;
  int j;

  if(!(in=fopen(filename, "r"))){
    printf("Problem opening file %s\n", filename);
    exit(1);
  }

  fscanf(in, "%d %d\n", &n_row, &n_cols);
  //printf("%d %d\n", n_row, n_cols);

  matrix = malloc(n_row * n_cols * sizeof(float));

  for(i=0;i<n_row;i++){
    for(j=0;j<n_cols;j++){
      fscanf(in, "%f", &matrix[i*n_cols + j]);
    }
  }    
  *n = n_row;
  *m = n_cols;
  return matrix;
}

void multimatri(float *matrix_A,float *matrix_B,int n_row,int n_cols){

int i,j,k;
k=0;
float *matrixtemp;
matrixtemp = malloc(n_row * n_cols * sizeof(float));

	for(i=0;i<n_row;i++){

		for(k=0;k<n_row;k++){  

		matrixtemp[i*n_cols + k] = 0.0;
  
			for(j=0;j<n_cols;j++){

			matrixtemp[i*n_cols + k] = matrixtemp[i*n_cols + k] + (matrix_A[i*n_cols + j])*(matrix_B[j*n_row + k]);

			}
		}
   	}
	for(i=0;i<n_row*n_cols;i++){

	matrix_A[i] = matrixtemp[i];
	}

}


void *soltrianinf(float *matrix_A,float *matrix_B,float *matrix_X,int n_rows,int n_cols){

float propunt;//contador productopunto
int i,j,cont,N;
N = n_rows*n_cols; 
	for(i=0;i<n_rows;i++){
	matrix_X[i]= 0.0;
	}

matrix_X[0] = matrix_B[0]/matrix_A[0];

	for(i=1;i< n_rows;i++){
	propunt = 0.0;
		for(j=0;j<i;j++){
		propunt = propunt + matrix_X[j]*matrix_A[i*n_cols + j];
		}

	matrix_X[i] = (matrix_B[i] - propunt)/matrix_A[i*n_cols + i];
	}
}

void *soltriansup(float *matrix_A,float *matrix_B,float *matrix_X,int n_rows,int n_cols){

float propunt;//contador productopunto
int i,j,cont,N;
N = n_rows*n_cols; 
	for(i=0;i<n_rows;i++){
	matrix_X[i]= 0.0;
	}

matrix_X[n_rows-1] = matrix_B[n_rows-1]/matrix_A[N-1];

	for(i=1;i< n_rows;i++){
	propunt = 0.0;
		for(j=0;j<i;j++){
		propunt = propunt + matrix_X[n_rows-1-j]*matrix_A[(n_rows-i-1)*n_cols + (n_cols-j-1)];
		}

	matrix_X[(n_cols-i-1)] = (matrix_B[(n_cols-i-1)] - propunt)/matrix_A[(n_rows-i-1)*n_cols + (n_cols-i-1)];
	}
}
