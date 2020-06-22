#include <iostream>
#include <iomanip>

#include <fstream>
#include <string>
#include <sstream>

#include <complex>
#include "LIB_dcomplex.h"

#include <eigen3/Eigen/Core>
#include <cmath>

#include <random>

using namespace std::complex_literals;
using namespace Eigen;
using namespace std;

int main(int argNum, char **argVec){ 
    /* Expected arguments: path,  trajOption, numTrajs, customTrajs
    - path:: 
        the path of the read file with the raw simulation data
    - trajOption:: 
        "N" if no trajectories to be plotted, 
        "R" if random trajectories to be plotted -according to the initial distribution, 
        "C" if Custom Trajectories will be plotted
    - numTrajs::
        Number of trajectories to plot. It will be 0 if trajOption=="No"
    - customTrajs::
        A sequence of the desired initial trajectories for the "YesC" option as string of numbers separated by spaces
    */
    if (argNum!=5){
        cout << "Error while reading the arguments. Too few arguments introduced? \n";
        return -1;
    }
    char trajOption;
    int numTrajs;
    sscanf(argVec[2], "%c", &trajOption);
    sscanf(argVec[3], "%d", &numTrajs);
    double * trajectories;
    trajectories=(double*) malloc(numTrajs*sizeof(double));
        
    ifstream readFile;
    ofstream writtenFile;
    readFile.open(argVec[1]);
    writtenFile.open("DATA_dataToPlot_1D_CN.txt");
    if(!readFile.good() || !writtenFile.good()){
        cout << "ERROR! Unable to open the read or write files! "<<endl;
        return 1;
    }
    int nx,nt;
    double xmin,xmax,dx, m, hbar, dt, norm;
    string line;
    cdouble arrayEl;
    getline(readFile,line);
    istringstream ss(line);
    ss >> nx >> nt >> dt >> xmin >> xmax >> m >> hbar; //nx is the number of space division, which is one less than the number of spatial points
    
    ArrayXcd currentPsi(nx+1);
    ArrayXd pos(nx+1), probDensity(nx+1);
    MatrixXd results(nx+1, 4);
    dx= (xmax-xmin)/nx;
    for(int xIt=0; xIt<=nx; xIt++){
        pos(xIt) = xmin+xIt*dx;
    }
    writtenFile << std::setprecision(17);

    if(trajOption=='N'){
    
    for(int tIt=1; tIt<=nt; tIt++){
        getline(readFile,line);
        for(int xIt=0; xIt<=nx; xIt++){
            getline(readFile,line);
            istringstream ss(line);
            ss>>arrayEl;
            currentPsi(xIt)=arrayEl;
        }
        
        //obtain the norm of the WF at this time
        probDensity=abs2(currentPsi);
        norm=(probDensity(0)+probDensity(nx))/2.0;
        for(int i=1; i<nx; ++i){norm+=probDensity(i);}
        norm*=dx;
        
        results << pos, probDensity, real(currentPsi), imag(currentPsi);
        writtenFile <<"Norm=" << norm << endl<<results<<endl<< endl << endl;
    }
    } else if (trajOption=='C'){
        istringstream ss(argVec[4]);
        for(int i=0; i<numTrajs; i++){
            ss >> trajectories[i];
        }
        
    } else { //randomly chosen trajectories according to the initial WF
        
        ArrayXd probab(nx+1);
        getline(readFile,line);
        for(int xIt=0; xIt<=nx; xIt++){
            getline(readFile,line);
            istringstream ss(line);
            ss>>arrayEl;
            currentPsi(xIt)=arrayEl;
        }
        probab=100*abs2(currentPsi);
        double* probabClist = probab.data();
        std::default_random_engine generator;
        std::discrete_distribution<int> distribution (probabClist,probabClist+nx);
        
        for(int i=0; i<numTrajs; i++){
            trajectories[i] = pos(distribution(generator));
        }
        
        readFile.seekg(0);
        
    }
    /*
    //computar trayectorias con la accion S
    if(trajOption=='R' || trajOption=='C'){
        double whole, fractional, baseDif, minDif;
        int kbase;
        ArrayXd velocityField(nx+1), phase(nx+1);

        ofstream trajFile;
        trajFile.open("DATA_trajectoriesToPlot_1D_CN.txt");
        if(!trajFile.good()){
            cout << "ERROR! Unable to open the read or write files! "<<endl;
            return 1;
        }
        
        for(int tIt=1; tIt<=nt; tIt++){
            getline(readFile,line);
            for(int xIt=0; xIt<=nx; xIt++){
                getline(readFile,line);
                istringstream ss(line);
                ss>>arrayEl;
                currentPsi(xIt)=arrayEl;
            }
            //sacar la fase compleja de currentPsi, sacar el velocity field y computar el siguiente tiempo de las trayectorias
            
            phase=arg(currentPsi);//hbar
            //lets make it a continous function a fuersa bruta
            kbase=0;
            for(int i=1; i<=nx;++i){
                baseDif=phase(i)-phase(i-1)+kbase*2*M_PI;
                minDif=min({abs(baseDif),abs(baseDif+2*M_PI),abs(baseDif-2*M_PI)});
                if(minDif==abs(baseDif)){
                }else if(minDif==abs(baseDif+2*M_PI)){ kbase++;
                }else{ kbase--;}
                phase(i)=phase(i)+kbase*2*M_PI;
            }
            velocityField(0)=hbar*(phase(1)-phase(0))/(m*dx);
            velocityField(1)=hbar*(phase(2)-phase(0))/(2*m*dx);
            for(int i=2; i<(nx-2); i++){
                velocityField(i)=hbar*(-phase(i+2)+8*phase(i+1)-8*phase(i-1)+phase(i-2))/(12*m*dx);
            }
            velocityField(nx-1)=hbar*(phase(nx)-phase(nx-2))/(2*m*dx);
            velocityField(nx)=hbar*(phase(nx)-phase(nx-1))/(m*dx);
            */
    //Compute the velocity field using that:
    /*
      Jk= hbar* Im{psi^-1 * diff(psi, qk)}/mk
      vk = Jk / psi**2
    */
            
    if(trajOption=='R' || trajOption=='C'){
        double whole, fractional;
        ArrayXd velocityField(nx+1);
        ArrayXcd auxComplexVector(nx+1), conjPsi(nx+1);

        ofstream trajFile;
        trajFile.open("DATA_trajectoriesToPlot_1D_CN.txt");
        if(!trajFile.good()){
            cout << "ERROR! Unable to open the read or write files! "<<endl;
            return 1;
        }
        
        for(int tIt=1; tIt<=nt; tIt++){
            getline(readFile,line);
            for(int xIt=0; xIt<=nx; xIt++){
                getline(readFile,line);
                istringstream ss(line);
                ss>>arrayEl;
                currentPsi(xIt)=arrayEl;
            }
            
            //Get the probability density function and the inverse psi
            
            probDensity = abs2(currentPsi);
            conjPsi= conj(currentPsi);
            
            // psi^-1* diff(psi,x) is computed
            
            auxComplexVector(0)=conjPsi(0)*(currentPsi(1)-currentPsi(0))/(dx*probDensity(0));
            auxComplexVector(1)=conjPsi(1)*(currentPsi(2)-currentPsi(0))/(2.0*dx*probDensity(1));
            for(int i=2; i<(nx-2); i++){
                auxComplexVector(i)=conjPsi(i)*(-currentPsi(i+2)+8.0*currentPsi(i+1)-8.0*currentPsi(i-1)+currentPsi(i-2))/(12.0*dx*probDensity(i));
            }
            auxComplexVector(nx-1)=conjPsi(nx-1)*(currentPsi(nx)-currentPsi(nx-2))/(2.0*dx*probDensity(nx-1));
            auxComplexVector(nx)=conjPsi(nx)*(currentPsi(nx)-currentPsi(nx-1))/(dx*probDensity(nx));
            
            // imaginary part is extracted and Jk obtained
            
            velocityField = (hbar/m)*imag(auxComplexVector);
            
            //velocity field obtained
               
            //obtain the norm
             //obtain the norm of the WF at this time
            norm=probDensity(0)/2.0+probDensity(nx)/2.0;
            for(int i=1; i<nx; ++i){norm+=probDensity(i);}
            norm*=dx;

            //print out the data of this iteration            
            results << pos, probDensity, real(currentPsi), imag(currentPsi);
            
            //and compute the trajectories
            
            for(int i=0; i<numTrajs;i++){
                //we apply the discretisation of the grid to the traj positions
                //trajectories[i] = trajectories[i]+velocityField( round((trajectories[i]-xmin)/dx) )*dt;
                //trajectories[i] = xmin + dx*round(( trajectories[i]+velocityField( round((trajectories[i]-xmin)/dx) )*dt -xmin)/dx);
                fractional = std::modf((trajectories[i]-xmin)/dx, &whole);
                if(whole>=nx){whole=nx-2;}else if(whole<0){whole=0;}
                trajectories[i] = trajectories[i]+( (1-fractional)*velocityField(whole)+fractional*velocityField(whole+1) )*dt;
                trajFile << trajectories[i] << " 0" << endl;
            }
            trajFile<<endl<<endl;
            writtenFile <<"Norm=" << norm << endl<<results<<endl<< endl << endl;
        }
       trajFile.close(); 
        
    }
    
    
    writtenFile.close();
    readFile.close();  
    return 0;
}
