#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float T = 40; // N
float L = 100; // m
float time = 120; //s

float *copy(float *in, int n);
void initial_conditions(float *u, int points_x, float dx);

int main(int argc, char **argv){
  float rho = atof(argv[1]);
  float c = sqrt(T/rho);
  int points_x = 101*10;
  float dx = L/(points_x-1);
  float dt = 0.05 * dx / c;

  float t = 0;
  float *u = malloc(points_x*sizeof(float));
  float *u_old = malloc(points_x*sizeof(float));
  float *u_old2 =  malloc(points_x*sizeof(float));
  char filename[50];
  sprintf(filename, "string_%.2f.dat", rho);
  FILE *out = fopen(filename, "w");
  float c_prime = dx/dt;

  int print = 0;

  initial_conditions(u, points_x, dx);
  int i;
  if (t>= print){
    for (i = 0; i <points_x/10; i++){
      fprintf(out, "%f ", u[i*10]);
    }
    fprintf(out, "\n");
    print += 1;
  }

  u_old = copy(u,points_x);

  t=t+dt;
  for (i=1; i<points_x-1; i++){
    u[i]= u_old[i] + 0.5*pow(c/c_prime,2)*(u_old[i+1]+u_old[i-1]-2*u_old[i]);
  }
  if (t>= print){
    for (i = 0; i <points_x/10; i++){
      fprintf(out, "%f ", u[i*10]);
    }
    fprintf(out, "\n");
    print += 1;
  }

  u_old2 = copy(u_old, points_x);
  u_old = copy(u,points_x);

  while (t<time){

    t=t+dt;
    for (i=1; i<points_x-1; i++){
      u[i]=2*u_old[i]-u_old2[i]+pow(c/c_prime,2)*(u_old[i+1] + u_old[i-1]-2*u_old[i]);
    }
  if (t>= print){
    for (i = 0; i <points_x/10; i++){
      fprintf(out, "%f ", u[i*10]);
    }
    fprintf(out, "\n");
    print += 1;
  }
    u_old2 = copy(u_old, points_x);
    u_old = copy(u,points_x);
  }


  return 0;
}

void initial_conditions(float *u, int points_x, float dx){
  int i;
  u[0]=0;
  u[points_x-1]=0;
  for (i=1; i<points_x-1; i++){
    if(i*dx <= 0.8*L)
      u[i]=1.25*i*dx/L;
    else if(i*dx>0.8*L)
      u[i]=5-5*i*dx/L;
  }
}

float *copy(float *in, int n){
  int i; 
  float *out = malloc(n*sizeof(float));
  for (i=0; i<n; i++)
    out[i]=in[i];
  return out;
}
