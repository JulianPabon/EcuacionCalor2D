#include<iostream>
#include<cmath>
#include<armadillo>

using namespace std;
using namespace arma;

//Llenar la matriz A de acuerdo a la formula
void llenarMatrizA(mat &A, int Nx, int Nz,double sx, double sy){
    double tmp = 1 + (2*sy) + (2*sx);
    int node,upper,lower;
    for(int x=1; x<Nz-1; ++x){        // iterate over rows
        for(int y=1; y<Nx-1; ++y){    //iterate over cols
            node = ( (x-1)*(Nx-2) ) + (y-1);
            upper = node + (Nx-2);
            lower = node - (Nx-2);
            A(node, node) =  tmp;
            if(x==1){
                A(node, upper) = -sy;
                if(y==1){                  //bottom-left
                    A(node, node+1) = -sx;
                }
                else if(y==Nx-2){         //bottom-right
                    A(node, node-1) =  -sx;
                }else{                    //bottom-middle
                    A(node, node-1) = -sx;
                    A(node, node+1) = -sx;
                }
            }
            else  if(x==Nz-2){
                A(node, lower) =  -sy;
                if(y==1){                 //top-left
                    A(node, node+1) = -sx;
                }
                else if(y==Nx-2){        //top-right
                    A(node, node-1) = -sx;
                }else{                   //top-middle
                    A(node, node-1) = -sx;
                    A(node, node+1) = -sx;
                }
            }
            else{
                A(node, lower) = -sy;
                A(node, upper) = -sy;
                if(y==1){                //left
                    A(node, node+1) = -sx;
                }
                else if(y==Nx-2){       //right
                    A(node, node-1) = -sx;
                }
                else{    //complete equation central points
                    A(node, node-1) = -sx;
                    A(node, node+1) = -sx;
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
    X((Nx * x + z)) = temp0;

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
