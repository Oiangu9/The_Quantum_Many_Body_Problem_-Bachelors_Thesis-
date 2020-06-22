#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include "LIB_dcomplex.h" // Macro for using dcomplex as std::complex<double> and J as the complex 0.0+i*1.0
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <cmath>
#include <random>
using namespace std::complex_literals;
using namespace Eigen;
using namespace std;
#define PI 3.141592653589793238463
#define INF 1000000.0
int nx=250, ny=150, jmax=10;
double xmin=-10.0, xmax=10.0, ymin=-10.0, ymax=10.0, dx, dy;
cdouble EigenstatesForSectionsInx(double y, double x, int j){ double d=7.0, L, o;if(x<0){ } else{ if(abs(y)>d){return 0.0;} L=2.0*d; o=-d; } return sqrt(2.0/L)*(-j*j*PI*PI/(L*L))*sin(j*PI*(y-o)/L);}
double potential(double x, double y){double sigmax=2.0, sigmay=2.0, mux=-6.0, muy=0.0, kx=0.0, ky=0.0; return pow(1.0/(2*PI*sigmax*sigmax),0.25)*exp(J*kx*x-0.25*pow((x-mux)/sigmax,2))*sin(PI*(y-ymin)/(ymax-ymin))*sqrt(2/(ymax-ymin));}
ArrayXd posx(nx+1), posy(ny+1);
int main(){
ofstream plotFile, potentialFile;
plotFile.open("DATA_adiabaticStatePlot.txt");
//Prepare the grid in x y
dx=(xmax-xmin)/nx;
dy=(ymax-ymin)/ny;
for(int ix=0; ix<=nx; ix++){
posx(ix) = xmin+ix*dx;
}
for(int iy=0; iy<=ny; iy++){
posy(iy) = ymin+iy*dy;
}
for(int ix=0; ix<=nx; ++ix){
for(int iy=0; iy<=ny; ++iy){
plotFile << posy(iy)<<" ";
for(int j=0; j<=jmax; ++j){
plotFile << real(EigenstatesForSectionsInx(posy(iy), posx(ix), j))<<" ";
}
plotFile<<endl;
}
plotFile<<endl<<endl; }
plotFile.close();
potentialFile.open("DATA_potentialToPlot.txt");
for(int i=0; i<=nx; i+=1){
for(int j=0; j<=ny; j+=1){
potentialFile << posx(i) << " " << posy(j) << " " << potential(posx(i), posy(j))<< endl;
}potentialFile<<endl;}
potentialFile.close();
return 0;
}
