#include<iostream>
#include<cmath>
#include "arrayfire.h"

using namespace std;
using namespace arma;

//Llenar la matriz A de acuerdo a la formula
//LA MATRIZ A SOLO SE CONTRUYE UNA VEZ POR ESO NO ES 
//NECESARIO PARALELIZARLA
//Llenar la matriz A de acuerdo a la formula
//Como los valores de los bordes ya se conocen, 
//la matriz A solo se construye con los nodos interiores
void llenarMatrizA(mat &A, int Nx, int Nz,double sx, double sz){
    double alfa = 1 + (2*sz) + (2*sx);
    int nodo, superior, inferior, derecha, izquierda;
    //Empieza en puntos z = x = 1, porque la cero es una condicion
    //de borde, si se inicia en cero la formula daria indices -1.
    //Lo anterior se aplica tambien para z = Nz -1 y  x = Nx -1

    //EMPEZAR ITERADORES EN 1, 1 PARA PODER SSABER CUANDO SE VA A 
    //ESTAR EN LOS NODOS DE FRONTERA
    for(int j = 1; j <= Nz; j++){        // iteterar sobre filas
        for(int i = 1; i <= Nx; i++){    // iterar sobre columnas
            nodo = ( (j-1)* Nx ) + (i-1);
            superior = nodo + Nx;
            inferior = nodo - Nx;
            derecha = nodo + 1;
            izquierda = nodo - 1;

            //Si i == 1 || i == Nx || j == 1 || j == Nz
            //Son nodos que limitan con los bordes y por lo tanto 
            //van a tener valores de cero

            //Ti+1,j, si i==Nx limita con un borde y por lo tanto el valor no se incluye en A
            if(i != Nx){
                A(nodo, derecha) = -sz;
            }
            //Ti,j+1, si j==Nz limita con un borde y por lo tanto el valor no se incluye en A
            if(j != Nz){
                A(nodo, superior) = -sx;
            }
            //Ti,j
            A(nodo, nodo) = alfa;
            //Ti-1,j, si i==1 limita con un borde y por lo tanto el valor no se incluye en A
            if(i != 1){
                A(nodo, izquierda) = -sz;
            }
            //Ti,j-1, si j==1 limita con un borde y por lo tanto el valor no se incluye en A
            if(j != 1){
                A(nodo, inferior) = -sx;
            }
        }
    }
}

//Funcion para imprimir solucion
//y graficar con python
void imprimirSolucion(mat &X, int Nx, int Nz, double deltaX, double deltaZ){
    for(double y=0; y< Nz; ++y){
        for(double x=0; x<Nx; ++x){
            //cout<<X(y*Nx+x)<<" ";
            cout<<x*deltaZ<<" "<<y*deltaX<<" "<<X(y*Nx+x)<<endl;
        }
        //cout<<endl;
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
    double deltaX,deltaZ,deltaT,k, temp0;
    int Nx,Nz,T,nodos, x, z;
    cin>>deltaX>>deltaZ>>deltaT>>Nx>>Nz>>T>>k>>x>>z;

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


    /* CODIGO ARRAYFIRE */
    //Se debe especificar la GPU del computador
    int device = 0;
    af::setDevice(device);
    // af::info();

    //Variables para pasar matriz de armadillo a arrayfire
    double *A_mem = (double*)malloc(nodos*nodos*sizeof(double));
    double *B_mem = (double*)malloc(nodos*sizeof(double));
    double *X_mem = (double*)malloc(nodos*sizeof(double));
    //Pasar matrices a punteros de c++
    A_mem = A.memptr();
    B_mem = B.memptr();
    X_mem = X.memptr();
    //Creacion de arrays en ArrayFire
    af::array afA(nodos,nodos,A_mem);
    af::array afB(nodos,f64);
    af::array afX(nodos,f64);

    //Poner temperatura inicial

    //Calculo de temperatura de la malla para cada tiempo
    for(int t = 0; t < T-1; t++){
        //En el tiempo 0, X contiene la malla con las
        //condiciones iniciales
        afB = afX;
        //Solucion del sistema de ecuaciones usando ArrayFire
        //X contiene la temperatura en el tiempo t   
        afX = af::solve(afA,afB);
    }
    return 0;
}

