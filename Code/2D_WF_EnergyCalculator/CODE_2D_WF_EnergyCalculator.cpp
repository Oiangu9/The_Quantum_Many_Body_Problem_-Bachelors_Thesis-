// 2D WF Energy Caluculator by using the Sandwich method- obv TIME INDEPENDENT POTENTIALS
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
using namespace std::complex_literals;
using namespace Eigen;
using namespace std;
const double PI = 3.141592653589793238463; const double EULER=2.718281828459045;

//defining the particle parameters
double m1=1.0, m2=1.0, hbar=1.0;
//define the SPACE GRID
double x1min=-10.0, x1max=10.0, x2min=-10.0, x2max=10.0;

//We define the initial WF
cdouble psi0(double x, double y){
 double sigmax=2.0, sigmay=2.0, mux=0.0, muy=0.0, kx=0.3, ky=0.0; return pow(1.0/(2*PI*sigmax*sigmax),0.25)*exp(J*kx*x-0.25*pow((x-mux)/sigmax,2))*sqrt(2/(x2max-x2min))*sin(PI*(y-x2min)/(x2max-x2min));;
}
//We define the Potential energy field TIME INDEPENDENT
double V(double x, double y){
 return 0;;
}
int main(){
cdouble Energy=0.0, norma=0.0;
int nx1=850, nx2=850;
double dx1=(x1max-x1min)/nx1, dx2=(x2max-x2min)/nx2;
int gridPoints=(nx1+1)*(nx2+1);

//Build the planed WavePacket (which should be seen as a 2d array
VectorXcd psi_ini(gridPoints); //nx, ny intervalos, then nx+1, ny+1 ptos por cada dim
for(int ix=0; ix<=nx1; ix++){
for (int iy=0; iy<=nx2;iy++){
psi_ini(ix*(nx2+1)+iy)=psi0(x1min+ix*dx1, x2min+iy*dx2);
}
}
//Build the Hamiltonian
SparseMatrix<cdouble> H(gridPoints, gridPoints);
H.reserve(VectorXi::Constant(gridPoints,5));
cdouble a=hbar/(2*m1*dx1*dx1), d=hbar/(2*m2*dx2*dx2);
int j;
for(int ix=0;ix<=nx1;ix++){
for(int iy=0; iy<=nx2;iy++){
j=ix*(nx2+1)+iy;
H.insert(j,j)=(hbar*hbar*(1/(m1*dx1*dx1)+1/(m2*dx2*dx2))+V(x1min+ix*dx1,x2min+iy*dx2))/((cdouble)hbar);
if(iy!=nx2) H.insert(j,j+1)=-d;
if(iy!=0) H.insert(j,j-1)=-d;
if(ix>0) H.insert(j,j-nx2-1)=-a;
if(ix<nx1) H.insert(j,j+nx2+1)=-a;
}
}
H.makeCompressed();
VectorXcd Hpsi(gridPoints), vt(gridPoints);
vt=psi_ini.conjugate();
Hpsi=H*psi_ini;
for(int i=0; i<gridPoints; ++i){Energy = Energy+ vt(i)*Hpsi(i);
 norma=norma+psi_ini(i)*vt(i);}
 Energy=Energy*dx1*dx2;
cout<<endl<<endl<<"Norm="<<norma*dx1*dx2<<endl<<endl;
//cout<<psi_ini<<endl<<endl<<psi_ini.conjugate();
ofstream outsideFile;
outsideFile.open("DATA_ObtainedEnergy.txt");
outsideFile << std::setprecision(17);
outsideFile << Energy<<endl; cout<<endl<<endl <<"Calculated Energy = "<<Energy<<endl<<endl;
outsideFile.close();
return 0;
}
