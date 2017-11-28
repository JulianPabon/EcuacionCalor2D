#include<iostream>
#include<cmath>
#include<armadillo>

using namespace std;
using namespace arma;

//Llenar la matriz A de acuerdo a la formula
void llenarMatrizA(mat &A, int Nx, int Nz,double sx, double sy){
    double tmp = 1 + (2*sy) + (2*sx);
    int nodo,superior,inferior;
    //Empieza en puntos z = x = 1, porque la cero es una condicion
    //de borde, si se inicia en cero la formula daria indices -1.
    //Lo anterior se aplica tambien para z = Nz -1 y  x = Nx -1
    for(int z=1; z<Nz-1; z++){        // iteterar sobre filas
        for(int x=1; x<Nx-1; x++){    // iterar sobre columnas
            //PROBLEMA CON SUPERIOR E INFERIOR INDICES NEGATIVOS
            nodo = ( (z-1)* Nx ) + (x-1);
            superior = nodo + Nx;
            inferior = nodo - Nx;
            A(nodo, nodo) =  tmp;
            if(z==1){
                A(nodo, superior) = -sy;
                if(x==1){                  //bottom-left
                    A(nodo, nodo+1) = -sx;
                }
                else if(x==Nx-2){         //bottom-right
                    A(nodo, nodo-1) =  -sx;
                }else{                    //bottom-middle
                    A(nodo, nodo-1) = -sx;
                    A(nodo, nodo+1) = -sx;
                }
            }else  if(z==Nz-2){
                A(nodo, inferior) =  -sy;
                if(x==1){                 //top-left
                    A(nodo, nodo+1) = -sx;
                }
                else if(x==Nx-2){        //top-right
                    A(nodo, nodo-1) = -sx;
                }else{                   //top-middle
                    A(nodo, nodo-1) = -sx;
                    A(nodo, nodo+1) = -sx;
                }
            }else{
                A(nodo, inferior) = -sy;
                A(nodo, superior) = -sy;
                if(x==1){                //left
                    A(nodo, nodo+1) = -sx;
                }
                else if(x==Nx-2){       //right
                    A(nodo, nodo-1) = -sx;
                }
                else{    //complete equation central points
                    A(nodo, nodo-1) = -sx;
                    A(nodo, nodo+1) = -sx;
                }
            }
        }
    }
}

int main(){
    
    /* 
        Hacer pruebas con variaciones de deltaX, deltaY
        El calculo y las matrices son de los tamaños interiores
        a los bordes.
    */
    //Variables
    /*  
        temp0: temperatura  a propagar
        x, z: puntos donde empieza la propagacion
    */
    double deltaX,deltaZ,deltaT,k, tem0;
    int Nx,Nz,T,nodos, x, z;
    cin>>deltaX>>deltaZ>>deltaT>>Nx>>Nz>>T>>k>>tem0>>x>>z;

    //Cantidad de puntos de la malla menos las filas y columnas
    //que componen las condiciones de borde y cuya temperatura es 0
    nodos = Nx * Nz;

    // ARMADILLO CREA LAS MATRICES Y VECTORES POR DEFECTO CON CEROS
    //Matriz tridiagonal A
    mat A = mat(nodos, nodos);
    //Vector B, puntos de la malla con temperaturas conocidas
    vec B = vec(nodos); 
    //AX=B;
    vec X = vec(nodos);
    // Condición inicial, lugar(x,z) de la malla donde se pone la temperatura tem0
    // x -> derecha - izquierda : horizontal
    // z -> arriba - abajo : vertial
    X((Nx * z + x)) = temp0;

    double sx = (k*deltaT)/(deltaX*deltaX);
    double sy = (k*deltaT)/(deltaZ*deltaZ);
    llenarMatrizA(A, Nx, Nz, sx, sy);

    //Calculo de temperatura de la malla para cada tiempo
    for(int t = 0; t < T-1; t++){
        //En el tiempo 0, X contiene la malla con la condicion inical
        B = X;
        //Solucion del sistema de ecuaciones usando armadillo
        //X contiene la temperatura en el tiempo t   
        X = solve(A,B);
    }
    return 0;
}
