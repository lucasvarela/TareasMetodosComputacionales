#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *f_fn(float *u);
float energy(float p, float rho, float v);
float *copy(float *in, int N);
float *add_vector(float *a, float *b);
float *scalar_mult(float lambda, float *a);
float pressure(float *u);

float g = 1.4;
float c = 3E8;

int main(int argc, char **argv){
  float p_I = 100000.0, p_D = 10000.0,rho_I = 1.0,rho_D = 0.125,gamma = 1.4;
  char *t_string = argv[1];
  float t = atof(argv[1]);
  int N_t,N_x,i,j;
  float dt,dx;
  N_x = 2000;
  N_t = 2000;
  dx = 20.0/N_x;
  dt=t/N_t;

  float **u = malloc(N_x*sizeof(float *));
  float **u_old = malloc(N_x*sizeof(float *));
  float **u_half = malloc(N_x*sizeof(float *));
  float **f = malloc(N_x*sizeof(float *));
  float **f_old = malloc(N_x*sizeof(float *));
  float **f_half = malloc(N_x*sizeof(float *));

  for(i=0; i<N_x; i++){
    u[i]=malloc(3*sizeof(float));
    u_old[i]=malloc(3*sizeof(float));
    u_half[i]=malloc(3*sizeof(float));

    f[i]=malloc(3*sizeof(float));
    f_old[i]=malloc(3*sizeof(float));
    f_half[i]=malloc(3*sizeof(float));

    u[i][1]=0;
    float rho = i<N_x/2.0 ? rho_I: rho_D;
    float p = i<N_x/2.0 ? p_I: p_D;
    u[i][0]= rho;
    u[i][2] = rho*energy(p,rho,0);
  }
  for (i=0; i<N_t;i++)
    u_old[i] = copy(u[i],3);

  for(j=0;j<N_t;j++){
    for(i=0;i<N_x;i++)
      f_old[i] = f_fn(u_old[i]);
    for(i=0;i<N_x-1;i++)
      u_half[i] = add_vector(scalar_mult(0.5,add_vector(u_old[i+1],u_old[i])),scalar_mult(-0.5*dt/dx,add_vector(f_old[i+1],scalar_mult(-1,f_old[i]))));
    for(i=0; i<N_x-1; i++)
      f_half[i] = f_fn(u_half[i]);  
    for (i=1; i<N_x-1;i++){
      u[i]=add_vector(u_old[i],scalar_mult(-dt/dx,add_vector(f_half[i],scalar_mult(-1.0,f_half[i-1]))));
    }
    for (i=0; i<N_x;i++){
      u_old[i] = copy(u[i],3);
    }
  }

  FILE *out;

  char filename[100];
  sprintf(filename, "estado_%s.dat", t_string);

  out = fopen(filename, "w");

  for(i=0;i<N_x;i++){
    fprintf(out, "%f %f %f %f\n", dx*i-10.0,  u[i][1]/u[i][0], pressure(u[i]), u[i][0]);
  }

  fclose(out);
}

float pressure(float *u){
  float out = (g-1)*(u[2]-0.5*pow(u[1],2)/u[0]);
  return out;
}

float *f_fn(float *u){
  float *out = malloc(3*sizeof(float));
  out[0]=u[1];
  out[1]=pow(u[1],2)/u[0]+(g-1)*(u[2]-0.5*pow(u[1],2)/u[0]);
  out[2]=(u[2]+(g-1)*(u[2]-0.5*pow(u[1],2)/u[0]))*u[1]/u[0];
  return out;
}


float energy(float p, float rho, float v){
  return p/((g-1)*rho)+0.5*pow(v,2);
}

float *copy(float *in, int N){
  float *out = malloc(N*sizeof(float));
  int i;
  for (i=0; i<N; i++)
    out[i]=in[i];
  return out;
}

// Vector addition
float *add_vector(float *a, float *b){
  float *c = malloc(3*sizeof(float));
  int i;
  for(i=0; i<3; i++)
    c[i]=a[i]+b[i];
  return c;
}

//Scalar times vector multiplication
float *scalar_mult(float lambda, float *a){
  float *new = malloc(3*sizeof(float));
  int i;
  for(i=0; i<3; i++)
    new[i]=lambda*a[i];
  return new;
}
