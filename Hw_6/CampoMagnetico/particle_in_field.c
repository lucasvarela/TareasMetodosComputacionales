#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double componente_vxB_z(double x, double y,double z,double v_y,double v_x, double cte1);
double componente_vxB_y(double x, double y,double z,double v_z,double v_x, double cte1);
double componente_vxB_x(double x, double y,double z, double v_y, double v_z, double cte1);
double Gamma(double v_x,double v_y,double v_z);


//Funcion principal
int main(int argc, char **argv){
FILE *fileout;
char filename[100000];
int i;
int j = 1;
double v_x_old,v_y_old,v_z_old,x_old,y_old,z_old,t;
double q_m,alpha_grados,alpha,v,dt,K_i,N,cte1,gamma;

cte1 = -3.0/( pow(10,5) ); // Se define esta constante por conveniencia -B_0* (RT a la tres) unidades en radios de la tierra y Teslas
q_m = 9.579*(pow(10,7) );  //cargamasa q/m en 1/(Teslas segundos)



//Posicion inicial
x_old = 2.0;
y_old = 0.0;
z_old = 0.0;

//Condiciones iniciales dadas de forma inconveniente
K_i = (1.065789/( pow(10,9) ) )*( atof(argv[1]) );
alpha_grados = atof(argv[2]);

//Se crea archivo que va a almacenar los resultados y se le da un nombre que depende de las condicione iniciales.
sprintf(filename, "trayectoria_%d_%d.dat", (int)atof(argv[1]) ,(int)alpha_grados);
fileout= fopen(filename, "w");
    
    
    
//Volviendo utiles esas condiciones iniciales
v = 47.0* sqrt( 1.0 - 1.0/( pow((K_i + 1.0 ),2) )  ) ; //magnitud de velocidad inicial en radios tierra sobre segundo
alpha = alpha_grados *( 3.1416 )/( 180.0 );
N = 10000000;
dt = 100.0/N;  // en segundos
int particion = N/1000;

//velocidad inicial
v_x_old = 0.0;
v_y_old = v*sin(alpha);
v_z_old = v*cos(alpha);

 
gamma = Gamma(v_x_old,v_y_old,v_z_old);

double modulo;
double v_ref  = sqrt(pow(v_x_old,2)+pow(v_y_old,2)+pow(v_z_old,2)); 
//modulo = sqrt(pow(v_x_old,2)+pow(v_y_old,2)+pow(v_z_old,2));
//printf("v=%f vx=%f vy=%f vz=%f\n",v, v_x_old, v_y_old, v_z_old);    

//RungeKuttaFourthOrderStep

double average_k_1_x,average_k_2_x,average_k_1_y,average_k_2_y,average_k_1_z,average_k_2_z,k_1_prime2_x,k_1_prime2_y,k_1_prime2_z;
double k_1_prime1_x,k_1_prime1_y,k_1_prime1_z,x_1,x_2,x_3,v_x_1,v_x_2,v_x_3,y_1,y_2,y_3,v_y_1,v_y_2,v_y_3,z_1,z_2,z_3,v_z_1,v_z_2,v_z_3;
double k_2_prime2_x,k_2_prime2_y,k_2_prime2_z,k_3_prime1_x,k_3_prime1_y,k_3_prime1_z,k_3_prime2_x,k_3_prime2_y,k_3_prime2_z;
double k_4_prime1_x,k_4_prime1_y,k_4_prime1_z,k_4_prime2_x,k_4_prime2_y,k_4_prime2_z,k_2_prime1_z,k_2_prime1_y,k_2_prime1_x;

for(i=1;i<N;i++){


// RungeKuttaFourthOrderStep
  
  if( ( pow(x_old,2) + pow(y_old,2) + pow(z_old,2) ) <= 1.0){
	printf("El proton colisiono con la tierra.");
	i = N + 42;	
	}
  
  
    k_1_prime2_x = ( (q_m)/gamma )*(componente_vxB_x(x_old,y_old,z_old,v_y_old,v_z_old,cte1) );
    k_1_prime2_y = ( (q_m)/gamma )*(componente_vxB_y(x_old,y_old,z_old,v_z_old,v_x_old,cte1) );
    k_1_prime2_z = ( (q_m)/gamma )*(componente_vxB_z(x_old,y_old,z_old,v_y_old,v_x_old,cte1) );

    k_1_prime1_x = v_x_old;
    k_1_prime1_y = v_y_old;
    k_1_prime1_z = v_z_old;

    //first step

    x_1 = x_old + 0.5 *dt * k_1_prime1_x;
    v_x_1 = v_x_old + 0.5 * dt * k_1_prime2_x;
    k_2_prime1_x = v_x_1;
    y_1 = y_old + 0.5 *dt * k_1_prime1_y;
    v_y_1 = v_y_old + 0.5 * dt * k_1_prime2_y;
    k_2_prime1_y = v_y_1;
    z_1 = z_old + 0.5 *dt * k_1_prime1_z;
    v_z_1 = v_z_old + 0.5 * dt * k_1_prime2_z;
    k_2_prime1_z = v_z_1;

    k_2_prime2_x = ( (q_m)/gamma )*(componente_vxB_x(x_1,y_1,z_1,v_y_1,v_z_1,cte1) );
    k_2_prime2_y = ( (q_m)/gamma )*(componente_vxB_y(x_1,y_1,z_1,v_z_1,v_x_1,cte1) );
    k_2_prime2_z = ( (q_m)/gamma )*(componente_vxB_z(x_1,y_1,z_1,v_y_1,v_x_1,cte1) );

    //second step
 
    x_2 = x_old + 0.5 *dt * k_2_prime1_x;
    v_x_2 = v_x_old + 0.5 * dt * k_2_prime2_x;
    k_3_prime1_x = v_x_2;
    y_2 = y_old + 0.5 *dt * k_2_prime1_y;
    v_y_2 = v_y_old + 0.5 * dt * k_2_prime2_y;
    k_3_prime1_y = v_y_2;
    z_2 = z_old + 0.5 *dt * k_2_prime1_z;
    v_z_2 = v_z_old + 0.5 * dt * k_2_prime2_z;
    k_3_prime1_z = v_z_2;

    k_3_prime2_x = ( (q_m)/gamma )*(componente_vxB_x(x_2,y_2,z_2,v_y_2,v_z_2,cte1) );
    k_3_prime2_y = ( (q_m)/gamma )*(componente_vxB_y(x_2,y_2,z_2,v_z_2,v_x_2,cte1) );
    k_3_prime2_z = ( (q_m)/gamma )*(componente_vxB_z(x_2,y_2,z_2,v_y_2,v_x_2,cte1) );
    
    //third
    
    x_3 = x_old + dt * k_3_prime1_x;
    v_x_3 = v_x_old + dt * k_3_prime2_x;
    k_4_prime1_x = v_x_3;
    y_3 = y_old + dt * k_3_prime1_y;
    v_y_3 = v_y_old + dt * k_3_prime2_y;
    k_4_prime1_y = v_y_3;
    z_3 = z_old + dt * k_3_prime1_z;
    v_z_3 = v_z_old + dt * k_3_prime2_z;
    k_4_prime1_z = v_z_3;

    k_4_prime2_x = ( (q_m)/gamma )*(componente_vxB_x(x_3,y_3,z_3,v_y_3,v_z_3,cte1) );
    k_4_prime2_y = ( (q_m)/gamma )*(componente_vxB_y(x_3,y_3,z_3,v_z_3,v_x_3,cte1) );
    k_4_prime2_z = ( (q_m)/gamma )*(componente_vxB_z(x_3,y_3,z_3,v_y_3,v_x_3,cte1) );

 
    //fourth step

    average_k_1_x = (1.0/6.0)*(k_1_prime1_x + 2.0*k_2_prime1_x + 2.0*k_3_prime1_x + k_4_prime1_x);
    average_k_2_x = (1.0/6.0)*(k_1_prime2_x + 2.0*k_2_prime2_x + 2.0*k_3_prime2_x + k_4_prime2_x);
    average_k_1_y = (1.0/6.0)*(k_1_prime1_y + 2.0*k_2_prime1_y + 2.0*k_3_prime1_y + k_4_prime1_y);
    average_k_2_y = (1.0/6.0)*(k_1_prime2_y + 2.0*k_2_prime2_y + 2.0*k_3_prime2_y + k_4_prime2_y);
    average_k_1_z = (1.0/6.0)*(k_1_prime1_z + 2.0*k_2_prime1_z + 2.0*k_3_prime1_z + k_4_prime1_z);
    average_k_2_z = (1.0/6.0)*(k_1_prime2_z + 2.0*k_2_prime2_z + 2.0*k_3_prime2_z + k_4_prime2_z);

    t = t + dt;
    x_old = x_old + dt * average_k_1_x;
    v_x_old= v_x_old + dt * average_k_2_x;
    y_old = y_old + dt * average_k_1_y;
    v_y_old= v_y_old + dt * average_k_2_y;
    z_old = z_old + dt * average_k_1_z;
    v_z_old= v_z_old + dt * average_k_2_z;



   if( i == (j*particion - 1)  ){
    modulo = sqrt(pow(v_x_old,2)+pow(v_y_old,2)+pow(v_z_old,2));
    fprintf(fileout,"%f %f %f\n",x_old, y_old, z_old);
	j = j + 1;
	}

}


//modulo = sqrt(pow(v_x_old,2)+pow(v_y_old,2)+pow(v_z_old,2));
    //printf("%f %f %f\n",modulo, y_old, z_old);  

//printf("%f\n", v_ref);
fclose(fileout);
return 0;

}

//funcion que calcula gamma
double Gamma(double v_x,double v_y,double v_z){

double respuesta, beta_cuad;

beta_cuad = ( pow(v_x,2) + pow(v_y,2) +pow(v_z,2) )/(pow(47,2));     //c definido en radios de la tierra / s cuad
respuesta = 1.0/( sqrt(1.0 - beta_cuad) );

return respuesta;
}


// funcion que calcula la componente z del producto cruz vxB
double componente_vxB_z(double x, double y,double z,double v_y,double v_x, double cte1){

double respuesta, r;

r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
respuesta = (1.0/( pow(r,5) ) )*cte1*z*3.0*( y*v_x - x*v_y );

return respuesta;

}


// funcion que calcula la componente y del producto cruz vxB
double componente_vxB_y(double x, double y,double z,double v_z,double v_x, double cte1){

double respuesta, r;

r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
respuesta = (1.0/( pow(r,5) ) )*cte1*( 3.0*x*z*v_z - ( 2.0*pow(z,2) - pow(x,2) - pow(y,2) )*v_x );

return respuesta;

}

// funcion que calcula la componente x del producto cruz vxB
double componente_vxB_x(double x, double y,double z, double v_y, double v_z, double cte1){

double respuesta, r;

r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
respuesta = cte1*( ( 2.0*pow(z,2) - pow(x,2) - pow(y,2) )*v_y - 3.0*y*z*v_z )/( pow(r,5) );

return respuesta;

}
