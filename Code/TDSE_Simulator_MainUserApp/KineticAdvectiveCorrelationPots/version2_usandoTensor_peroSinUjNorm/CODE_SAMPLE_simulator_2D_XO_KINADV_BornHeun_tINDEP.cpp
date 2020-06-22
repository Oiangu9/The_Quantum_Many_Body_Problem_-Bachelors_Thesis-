// 2D SCHRODINGER EQUATION SOLVER - XO ALGORITHM Kinetic and Advective Correlation Potential approximation of G and J:

// The Conditional Single Particle Wave Function (CSPWF) of the dimensions (x,y) will be evolved for each initial conditions using a 1D Cranck Nicolson method for the Pseudo Schrodinger Equations

//We include the necessary libraries

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

//USER INPUT DECLARATION----------------------------------------------------------------------

//We declare the Spatial and Time grids - names are self-explaning
double xmin = -6.0, xmax = 14.0, ymin = -7.0, ymax = 7.0, t0=0.0, dt=0.001, posx, posy, Nx, Ny;
cdouble Uj;
int xDivs=700, yDivs = 700, timeIts=2700, aux;
double dx=(xmax-xmin)/xDivs;
double dy=(ymax-ymin)/yDivs;
int outputDataEvery = 24; // data will be outputed to the external file only every x time iterations, if 1, every time iteration will be outputted

//The constants: hbar and the masses of each dimension are declared
double hbar=1.0, mx=1.0, my=1.0;

//We declare the initial full Wave Function (x,y,t=t0). It will be used in order to obtain the initial CSPWF by evaluating at each dimension the initial position of the other trajectory
cdouble initialFullWF(double x, double y){
    double sigmax=1.0, sigmay=1.0, mux=-2.0, muy=0.0, kx=4.50, ky=0.0; 
    return pow(1.0/(2*PI*sigmax*sigmax),0.25)*exp(J*kx*x-0.25*pow((x-mux)/sigmax,2)) * pow(1.0/(2*PI*sigmay*sigmay),0.25)*exp(J*ky*y-0.25*pow((y-muy)/sigmay,2));
}

//We declare the static external Potential Field
double W(double x,double y){
    double a=2.5, b=3.5, d=1.0;
    if(((x>=a) & (x<=b)) & (y>=d || y<=-d)){ return 1000000.0;}else {return 0.0;}
}

//The highest used adiabatic states' index in the Born Heun expansion's truncated version
int xjmax = 50, yjmax = 5;
double advectiveCor, kineticCor;
cdouble correlPot;

//The "analytical" expression for the adiabatic states for each x as a function of y and j are expressed
double eigenstatesForSectionsInx(double y, double x, int j){ //the so called psi_x^j(y)
    double c=2.5, d=3.5, s=1.0, L, o;
    if(x<c || x>d){ L=ymax-ymin; o=ymin;} 
    else{ if(abs(y)>s){return 0.0;} L=2.0*s; o=-s; } 
    return sqrt(2.0/L)*sin(j*PI*(y-o)/L);
}

//The "analytical" expression for the adiabatic states for each y as a function of x and j are expressed
double eigenstatesForSectionsIny(double x, double y, int j){ //the so called psi_y^j(x)
    double c=2.5, d=3.5, s=1.0, L, o;  
    if((x<c) & (y<-s || y>s)){ L = c-xmin; o=xmin;}  
    else if((x>d) & (y<-s || y>s)){L = xmax-d; o=d;}  
    else{ if(y<-s || y>s){return 0.0;} L=xmax-xmin; o=xmin;}  
    return sqrt(2.0/L)*sin(j*PI*(x-o)/L);
}

//The analytical expression for the above functions' first and second derivatives
double diffyEigenstatesForSectionsInx(double y, double x, int j){ //the so called d psi_x^j(y)/dy
    double c=2.5, d=3.5, s=1.0, L, o;
    if(x<c || x>d){ L=ymax-ymin; o=ymin;} 
    else{ if(abs(y)>s){return 0.0;} L=2.0*s; o=-s; } 
    return sqrt(2.0/L)*(j*PI/L)*cos(j*PI*(y-o)/L);
}

double diffyyEigenstatesForSectionsInx(double y, double x, int j){ //the so called d**2 psi_x^j(y)/dy**2
    double c=2.5, d=3.5, s=1.0, L, o;
    if(x<c || x>d){ L=ymax-ymin; o=ymin;} 
    else{ if(abs(y)>s){return 0.0;} L=2.0*s; o=-s; } 
    return sqrt(2.0/L)*(-j*j*PI*PI/(L*L))*sin(j*PI*(y-o)/L);
}

double diffxEigenstatesForSectionsIny(double y, double x, int j){ //the so called d psi_y^j(x)/dx
    double c=2.5, d=3.5, s=1.0, L, o;  
    if((x<c) & (y<-s || y>s)){L = c-xmin; o=xmin;}  
    else if((x>d) & (y<-s || y>s)){L = xmax-d; o=d;} 
    else{ if(y<-s || y>s){return 0.0;} L=xmax-xmin; o=xmin;} 
    return sqrt(2.0/L)*(j*PI/L)*cos(j*PI*(x-o)/L);
}

double diffxxEigenstatesForSectionsIny(double y, double x, int j){ //the so called d**2 psi_y^j(x)/dx**2
    double c=2.5, d=3.5, s=1.0, L, o;
    if((x<c) & (y<-s || y>s)){L = c-xmin; o=xmin;} 
    else if((x>d) & (y<-s || y>s)){L = xmax-d; o=d;}  
    else{ if(y<-s || y>s){return 0.0;} L=xmax-xmin; o=xmin;} 
    return sqrt(2.0/L)*(-j*j*PI*PI/(L*L))*sin(j*PI*(x-o)/L);
}
//The grid positions in x and y are saved into a vector as they are heavily used in the computation of Uj
ArrayXd xgrid(xDivs+1), ygrid(yDivs+1);
int main(){
for(int k=0; k<=xDivs; ++k){xgrid(k)=xmin+k*dx;}
for(int k=0; k<=yDivs;++k){ygrid(k)=ymin+dy*k;}

//The initial positions of each trajectory that will be evolved using the algorithm are chosen according to the probability distribution given by the modulus squared of the initial wave function
int numTrajs=20; // we choose the number of trajectories that will be evolved
int gridPointsFullWF=(xDivs+1)*(yDivs+1);

double fractional, whole;
ArrayXd probabDensity(gridPointsFullWF);
ArrayXcd initialFullPsi(gridPointsFullWF);

double* initialPosx = (double*) malloc(numTrajs*sizeof(double));
double* initialPosy = (double*) malloc(numTrajs*sizeof(double));

// the initial state of the full wavefunction is generated in order to obtain its modulus squared in each point
for(int i=0; i<=xDivs; ++i){
    for(int j=0; j<=yDivs; ++j){
        initialFullPsi(i*(yDivs+1) + j) = initialFullWF(xmin+i*dx, ymin+j*dy);
    }
}

// the probability associated with each space point is generated
probabDensity =100*abs2(initialFullPsi);

// the random number generator initialised

double* probabClist = probabDensity.data();
std::default_random_engine generator;
std::discrete_distribution<int> distribution (probabClist,probabClist+gridPointsFullWF-1);

// and the initial positions chosen according to the associated probabilities
for(int i=0; i<numTrajs; i++){
    aux=distribution(generator); //returns the winner index of the prob vector -> we must revert the indexing to 2d indexes
    initialPosx[i] = xmin + ((int) aux/(yDivs+1))*dx;
    initialPosy[i] = ymin + (aux%(yDivs+1))*dy;
}


// begin the time iterations for each evolved trajectory - UL, UR must be renamed in every iteration - as if it was a time dependant potential algorithm for a 1D particle
//we declare and prepare the propagator matrices for the Cranck Nicolson (CN) evolution
SparseMatrix<cdouble> U1x(xDivs+1, xDivs+1), U2x(xDivs+1, xDivs+1);
SparseMatrix<cdouble> U1y(yDivs+1, yDivs+1), U2y(yDivs+1, yDivs+1);
U1x.reserve(VectorXi::Constant(xDivs+1,3));
U2x.reserve(VectorXi::Constant(xDivs+1,3));
U1y.reserve(VectorXi::Constant(yDivs+1,3));
U2y.reserve(VectorXi::Constant(yDivs+1,3));

cdouble ax=J*dt*hbar/(4.0*mx*dx*dx);
cdouble ay=J*dt*hbar/(4.0*my*dy*dy);

for(int i=1;i<xDivs;i++){
    U1x.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1x.insert(i,i+1)= -ax;
    U1x.insert(i,i-1)= -ax;
}
U1x.insert(0,0)= 1.0*J;
U1x.insert(0,1)= -ax;
U1x.insert(xDivs,xDivs)= 1.0*J;
U1x.insert(xDivs,xDivs-1)= -ax;
U1x.makeCompressed();
U2x = U1x.conjugate();
U2x.makeCompressed();

for(int i=1;i<yDivs;i++){
    U1y.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1y.insert(i,i+1)= -ay;
    U1y.insert(i,i-1)= -ay;
}
U1y.insert(0,0)= 1.0*J;
U1y.insert(0,1)= -ay;
U1y.insert(yDivs,yDivs)= 1.0*J;
U1y.insert(yDivs,yDivs-1)= -ay;
U1y.makeCompressed();
U2y = U1y.conjugate();
U2y.makeCompressed();


//We initialise the LU solvers

SparseLU<SparseMatrix<cdouble>> LUsolverx;
SparseLU<SparseMatrix<cdouble>> LUsolvery;


//We declare the SPCWFs and the auxiliar vectors for velocity field computation
VectorXcd psiX(xDivs+1), psiY(yDivs+1), U2psix(xDivs+1), U2psiy(yDivs+1);
VectorXcd conjPsix(xDivs+1), conjPsiy(yDivs+1), auxX(xDivs+1), auxY(yDivs+1);
ArrayXd probDensityx(xDivs+1), probDensityy(yDivs+1), velocityFieldx(xDivs+1), velocityFieldy(yDivs+1);


//we define the trajectory matrix
double** traj=new double*[timeIts+1];
for (int i=0; i<=timeIts; ++i){ traj[i]= new double[2];} // the trajectory is saved in an array of timeIts arrays of 2 doubles (xi, yi)

// each of the timeIts arrays contains the value for the trajectory in each of the x,y at that iteration
double vx, vy;

//We open the output streams
ofstream probabDataFile, trajDataFile;
//psiDataFile.open("DATA_rawSimulationData_nD_XO_ZERO_CN_ABC_tDEP.txt");
probabDataFile.open("DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt");
trajDataFile.open("DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt");
//psiDataFile << std::setprecision(17);
probabDataFile << std::setprecision(17);
trajDataFile << std::setprecision(17);

//BEGINNING OF THE ALGORITHM FOR EACH OF THE INITIAL CONDITIONS----------------------------------------------------------------

for(int trajNum=0; trajNum<numTrajs; ++trajNum){ //this is a potential multithreading branching point

//We initialise the SPCWF conditioning the full WF to the intial values of the trajectories of this iteration
for(int i=0; i<=xDivs; ++i){
    psiX(i) = initialFullWF(xmin+i*dx, initialPosy[trajNum]);
}
for(int i=0; i<=yDivs; ++i){
    psiY(i) = initialFullWF(initialPosx[trajNum], ymin+i*dy);
}
traj[0][0]=initialPosx[trajNum];
traj[0][1]=initialPosy[trajNum];

// TIME ITERATIONS BEGIN -----------------------------------------------------------------------------------------------------

for(int it=0; it<timeIts; ++it){
    //NEXT POSITION OF THE TRAJECTORY OBTAINED -----------------------------
    //Using the current SPCWF, the velocity field of each dimension at this time is obtained and the next position of the trajectory calculated

    //first for the x dimension-------------------------------------

    //we first get the probability density function and the inverse psi
    probDensityx =abs2(psiX.array());
    conjPsix = conj(psiX.array());

    // psi_qi^-1* diff(psi_qi,qi) is computed:
    //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
    auxX(0) = conjPsix(0)*(psiX(1)-psiX(0))/(dx*probDensityx(0));
    auxX(1) = conjPsix(1)*(psiX(2)-psiX(0))/(2.0*dx*probDensityx(1));
    //the rest of points are got with a o(dq^4) difference
    for(int i=2; i<=xDivs-2; ++i){
    auxX(i) = conjPsix(i)*(-psiX(i+2) + 8.0*psiX(i+1) - 8.0*psiX(i-1) +psiX(i-2))/(12*dx*probDensityx(i));
    }
    auxX(xDivs-1) = conjPsix(xDivs-1)*(psiX(xDivs)-psiX(xDivs-2))/(2.0*dx*probDensityx(xDivs-1));
    auxX(xDivs) = conjPsix(xDivs)*(psiX(xDivs)-psiX(xDivs-1))/(dx*probDensityx(xDivs));

    // imaginary part is extracted and the velocity field obtained
    velocityFieldx = (hbar/mx)*imag(auxX.array());

    //now the y dimension------------------------------------

    //we first get the probability density function and the inverse psi
    probDensityy =abs2(psiY.array());
    conjPsiy = conj(psiY.array());

    // psi_qi^-1* diff(psi_qi,qi) is computed:
    //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
    auxY(0) = conjPsiy(0)*(psiY(1)-psiY(0))/(dy*probDensityy(0));
    auxY(1) = conjPsiy(1)*(psiY(2)-psiY(0))/(2.0*dy*probDensityy(1));
    //the rest of points are got with a o(dq^4) difference
    for(int i=2; i<=yDivs-2; ++i){
    auxY(i) = conjPsiy(i)*(-psiY(i+2) + 8.0*psiY(i+1) - 8.0*psiY(i-1) +psiY(i-2))/(12*dy*probDensityy(i));
    }
    auxY(yDivs-1) = conjPsiy(yDivs-1)*(psiY(yDivs)-psiY(yDivs-2))/(2.0*dy*probDensityy(yDivs-1));
    auxY(yDivs) = conjPsiy(yDivs)*(psiY(yDivs)-psiY(yDivs-1))/(dy*probDensityy(yDivs));

    // imaginary part is extracted and the velocity field obtained
    velocityFieldy = (hbar/my)*imag(auxY.array());

    //we apply the discretisation of the grid to the traj positions
    fractional = std::modf((traj[it][0]-xmin)/dx, &whole);
    if(whole>=xDivs){whole=xDivs-2;}else if(whole<0){whole=0;}
    vx=(1-fractional)*velocityFieldx(whole)+fractional*velocityFieldx(whole+1);
    traj[it+1][0] = traj[it][0]+vx*dt;
    fractional = std::modf((traj[it][1]-ymin)/dy, &whole);
    if(whole>=yDivs){whole=yDivs-2;}else if(whole<0){whole=0;}
    vy= (1-fractional)*velocityFieldy(whole)+fractional*velocityFieldy(whole+1);
    traj[it+1][1] = traj[it][1]+vy*dt;
    //The Ui propagator matrices of each dimension x,y are updated for this time iteration
    posy = traj[it][1];
    for(int i=0; i<=xDivs; ++i){
        posx = xmin +i*dx;
        correlPot =0.0;
        for(int j=0; j<=xjmax; ++j){ //generate the kinetic and advective correlation potentials for this spatial grid point posx
            kineticCor = -0.5*diffyyEigenstatesForSectionsInx(posy, posx, j);
            advectiveCor = vy*diffyEigenstatesForSectionsInx(posy, posx, j);
            Uj=0.5*(eigenstatesForSectionsInx(ymin,posx,j)*psiY(0) + eigenstatesForSectionsInx(ymax,posx,j)*psiY(yDivs));
            for(int k=1; k<yDivs; ++k){ 
                Uj=Uj+eigenstatesForSectionsInx(ygrid(k),posx,j)*psiY(k);}    
            Uj=Uj*dy/((cdouble) sqrt(Nx*Ny)); 
            correlPot = correlPot + ( Uj )*(kineticCor + J*advectiveCor);
        }
        U1x.coeffRef(i,i) = 1.0+J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) + correlPot )/((cdouble)2.0*hbar);
        U2x.coeffRef(i,i) = 1.0-J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) + correlPot )/((cdouble)2.0*hbar);
    }
    posx = traj[it][0];
    for(int i=0; i<=yDivs; ++i){
        posy = ymin +i*dy;
        correlPot =0.0;
        for(int j=0; j<=yjmax; ++j){ //generate the kinetic and advective correlation potentials for this spatial grid point posx
            kineticCor = -0.5*diffxxEigenstatesForSectionsIny(posx, posy, j);
            advectiveCor = vx*diffxEigenstatesForSectionsIny(posx, posy, j);
            Uj=0.5*(eigenstatesForSectionsIny(xmin,posy,j)*psiX(0) + eigenstatesForSectionsIny(xmax,posy,j)*psiX(xDivs));
            for(int k=1; k<xDivs; ++k){ 
                Uj=Uj+eigenstatesForSectionsIny(xgrid(k),posy,j)*psiX(k);}    
            Uj=Uj*dx/((cdouble) sqrt(Nx*Ny));
            correlPot = correlPot + ( Uj )*(kineticCor + J*advectiveCor);
        }
        U1y.coeffRef(i,i)= 1.0+J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) + correlPot)/((cdouble)2.0*hbar);
        U2y.coeffRef(i,i)= 1.0-J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) + correlPot)/((cdouble)2.0*hbar);
    }
    U1x.makeCompressed();
    U1y.makeCompressed();
    U2x.makeCompressed();
    U2y.makeCompressed();

    //LU decomposition done
    LUsolverx.compute(U1x);
    if(LUsolverx.info()!=Success) {
    cout << "LUx decomposition FAILED!" << endl;
    return 1;
    }
    U2psix= U2x*psiX;
    psiX = LUsolverx.solve(U2psix); //the wavefunction of the next time iteration is generated

    LUsolvery.compute(U1y);
    if(LUsolvery.info()!=Success) {
    cout << "LUy decomposition FAILED!" << endl;
    return 1;
    }
    U2psiy= U2y*psiY;
    psiY = LUsolvery.solve(U2psiy);

    if( it%outputDataEvery == 0){ //then we output the data
        probabDataFile << probDensityx << endl << endl<<endl;
        probabDataFile << probDensityy << endl << endl<<endl; 
        }
    }

    for(int it=0; it<timeIts; ++it){
        if( it%outputDataEvery == 0){ //then we output the data
            trajDataFile << traj[it][0] << " " << traj[it][1] << " ";
            trajDataFile <<" 0"<<endl;
        }
    }
}

probabDataFile.close();
trajDataFile.close();


//We output the shape of the potential in order to be able to plot it
ofstream potentialToPlot;

//we define some output finnes parameters in case it is not necessary to plot the potential to full accuracy (make the algorithm faster)
double potentialPlotFinness=0.013;
int enoughStepx=xDivs*potentialPlotFinness;
int enoughStepy=yDivs*potentialPlotFinness;
double* posArx=new double[xDivs+1];
double* posAry=new double[yDivs+1];

for(int i=0; i<=xDivs; i+=enoughStepx){
    posArx[i]=xmin+i*dx;}

for(int j=0; j<=yDivs; j+=enoughStepy){
    posAry[j]=ymin+j*dy;}

potentialToPlot.open("DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt");

for(int it=0; it<=timeIts; ++it){
    for(int i=0; i<=xDivs; i+=enoughStepx){
        for(int j=0; j<=yDivs; j+=enoughStepy){
            potentialToPlot << posArx[i] << " " << posAry[j] << " " << W(posArx[i], posAry[j])<< endl;
        }
        potentialToPlot<<endl;
    }
    potentialToPlot<<endl;
}
potentialToPlot.close();
return 0;
}
