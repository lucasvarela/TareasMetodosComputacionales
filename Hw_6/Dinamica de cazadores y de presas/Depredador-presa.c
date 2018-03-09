

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float funcion_prima_x(float x, float y);
float funcion_prima_y(float x, float y);
void Runge_kutta(float x,float y,float t,float paso,float *siguiente_iteracion);
void imprimir (float *x ,float *y,float *t,int total);



int main(int argc, char *argv[])

{/*punteros que almacenan posiciones y tiempo*/
    float *x;
    float*y;
    float*t;
    /*puntero que contiene informacion de cada iteracion*/
    float*siguiente_iteracion;
    /*posicion inicial y final y distancia de los pasos*/
    float minimo = 0.0;
    float maximo = 1.0;
    float paso=0.0001;
    int i;
    float total=(int)((maximo-minimo)/paso);
    x=malloc(total*sizeof(float));
    y=malloc(total*sizeof(float));
    t=malloc(total*sizeof(float));
    siguiente_iteracion=malloc(3*sizeof(float));
    x[0]=atof(argv[1]);
    t[0]=0.0;
    y[0]=atof(argv[2]);

    for (i=0; i<total; i++) {
        Runge_kutta(x[i], y[i], t[i],paso,siguiente_iteracion);
        x[i+1]=siguiente_iteracion[0];
        y[i+1]=siguiente_iteracion[1];
        t[i+1]=siguiente_iteracion[2];
    }


    imprimir(x,y,t,total);
    
}

/*pendientes para x y y*/

float funcion_prima_x(float x, float y)
{
    return (20.0)*x-x*y;


}

float funcion_prima_y(float x, float y){
    return -(30.0)*y+x*y;
    
    
    
}
/*se realiza un runge kutta de 4 ordern en x y y simultaneamente ya que las dos son funciones de t. Se */
void Runge_kutta(float x,float y,float t,float paso,float *siguiente_iteracion)
{   float k1,l1,k2,l2,k3,l3,k4,l4;
 

    /*primer paso*/
    k1=funcion_prima_x(x,y)*paso;
    l1=funcion_prima_y(x,y)*paso;
    
    /*segundo paso*/
    k2=funcion_prima_x(x+k1/2,y+l1/2)*paso;
    l2=funcion_prima_y(x+k1/2,y+l1/2)*paso;

    /*tercer paso*/
    k3=funcion_prima_x(x+k2/2,y+l2/2)*paso;
    l3=funcion_prima_y(x+k2/2,y+l2/2)*paso;
    
    /*cuarto paso*/
    k4=funcion_prima_x(x+k3,y+l3)*paso;
    l4=funcion_prima_y(x+k3,y+l3)*paso;
    
    /*promedio para x y y*/
    float x_promedio= (k1+(2.0)*k2+(2.0)*k3+k4)/(6.0);
    float y_promedio= (l1+(2.0)*l2+(2.0)*l3+l4)/(6.0);

    /*se suma el  promedio para avanzar*/
    float x_final=x+x_promedio;
    float y_final=y+y_promedio;
    float t_final=t+paso;
    /*la nueva posicion para x,y, y t se almacena en un puntero el cual pasa la informacion al puntero que contiene las posiciones de y,x, y t*/
    siguiente_iteracion[0]=x_final;
    siguiente_iteracion[1]=y_final;
    siguiente_iteracion[2]=t_final;


}

/*se imprime información en un archivo*/
void imprimir (float *x ,float *y,float *t,int total)

{
    FILE *fileout;
    int i;
    char filename[100000];
    sprintf(filename, "poblaciones_%d_%d.dat", (int)x[0],(int)y[0]);

    
    fileout= fopen(filename, "w");
    for (i=0;i <total;i++)
    {
        fprintf(fileout,"%f %f\n",x[i],y[i]);
    }
    fclose(fileout);
}






