#include <iostream>
#include <iomanip>

#include <fstream>
#include <string>
#include <sstream>

#include <complex>
#include "LIB_dcomplex.h"

#include <eigen3/Eigen/Core>

#include <random>

using namespace std::complex_literals;
using namespace Eigen;
using namespace std;

int main(int argNum, char **argVec){
    if (argNum!=5){
        cout << "Error while reading the arguments. Too few arguments introduced? \n";
        return -1;
    }
    char trajOption;
    int numTrajs;
    sscanf(argVec[2], "%c", &trajOption);
    sscanf(argVec[3], "%d", &numTrajs);
    double * trajectoriesx1, * trajectoriesx2;
    int aux;
    trajectoriesx1=(double*) malloc(numTrajs*sizeof(double));
    trajectoriesx2=(double*) malloc(numTrajs*sizeof(double));
    
    ifstream readFile;
    ofstream writtenFile;
    readFile.open(argVec[1]);
    writtenFile.open("DATA_dataToPlot_2D_CN.txt");
    if(!readFile.good() || !writtenFile.good()){
        cout << "ERROR! Unable to open the read or write files! "<<endl;
        return 1;
    }
    int nx1,nx2,nt,gridPoints;
    double x1min,x1max,x2min,x2max,dx1,dx2,m1,m2,hbar,dt, norm, sumCorners, sumBorders, sumInterior;
    string line;
    cdouble arrayEl;
    getline(readFile,line);
    istringstream ss(line);
    ss >> nx1 >> nx2 >> nt >> dt >>x1min >> x1max >> x2min >> x2max>>m1>>m2>>hbar; //nx is the number of space division, which is one less than the number of spatial points
    gridPoints=(nx1+1)*(nx2+1);
    ArrayXcd currentPsi(gridPoints);
    ArrayXd posx1(gridPoints), posx2(gridPoints), absolute(gridPoints);
    MatrixXd results(nx2+1, 3);
    dx1=(x1max-x1min)/nx1;
    dx2=(x2max-x2min)/nx2;
    for(int ix=0; ix<=nx1; ix++){
        for(int iy=0; iy<=nx2; iy++){
            posx1(ix*(nx2+1)+iy) = x1min+ix*dx1;
            posx2(ix*(nx2+1)+iy) = x2min+iy*dx2;
        }
    }
    writtenFile << std::setprecision(17);
    if(trajOption=='N'){
    for(int tIt=1; tIt<=nt; tIt++){
        getline(readFile,line);
        for(int xIt=0; xIt<gridPoints; xIt++){
            getline(readFile,line);
            istringstream ss(line);
            ss>>arrayEl;
            currentPsi(xIt)=arrayEl;
        }
        absolute=abs2(currentPsi);
        
        // calculate the norm in this time iteration using the 2d trapezium rule
        
        sumCorners = 0;
        sumBorders = 0;
        sumInterior = 0;
        
        sumCorners = absolute(0) + absolute(nx2) + absolute(gridPoints-1) + absolute(nx1*(nx2+1));
        for(int ix=1; ix<nx1; ++ix){
            for(int iy=1; iy<nx2; ++iy){
                sumInterior += absolute(ix*(nx2+1)+iy);
            }
        }
        for(int iy=1; iy<nx2; ++iy){
            sumBorders += absolute(iy) + absolute(nx1*(nx2+1) + iy); // ix=0 row and ix=nx1 row
        }
        for(int ix=1; ix<nx1; ++ix){
            sumBorders += absolute(ix*(nx2+1)) + absolute(ix*(nx2+1)+nx2); //iy=0 row and iy=nx2 row
        }
        
        norm=0.25*dx1*dx2*(sumCorners + 2*sumBorders +4* sumInterior);
        writtenFile <<"Norm="<< norm << endl;
        
        //output the position together with the probability
        for(int ix=0; ix<=nx1; ix++){
            results << posx1.block(ix*nx2+ix,0, nx2+1, 1), posx2.block(ix*nx2+ix,0, nx2+1, 1), absolute.block(ix*nx2+ix,0, nx2+1, 1);
            writtenFile <<results<<endl<< endl;
        }
        writtenFile <<endl;
    }
    }else if(trajOption=='C'){
        // vField: R2->R2; vField(x1,x2)=(vFieldx1(x1,x2),vFieldx2(x1,x2)) with vFieldxi:R2->R
        
        istringstream ss(argVec[4]);
        for(int i=0; i<numTrajs; i++){
            ss >> trajectoriesx1[i];
            ss >> trajectoriesx2[i];
        }
                
        
    }else{
        // vField: R2->R2; vField(x1,x2)=(vFieldx1(x1,x2),vFieldx2(x1,x2)) with vFieldxi:R2->R
    
        for(int xIt=0; xIt<gridPoints; xIt++){
            getline(readFile,line);
            istringstream ss(line);
            ss>>arrayEl;
            currentPsi(xIt)=arrayEl;
        }
        absolute=100*abs2(currentPsi);
        double* probabClist = absolute.data();
        std::default_random_engine generator;
        std::discrete_distribution<int> distribution (probabClist,probabClist+gridPoints-1);
        
        for(int i=0; i<numTrajs; i++){
            aux=distribution(generator);
            trajectoriesx1[i] = posx1(aux);
            trajectoriesx2[i] = posx2(aux);
        }
        
        readFile.seekg(0);
    }
    if(trajOption=='C' || trajOption=='R'){
        double wholex1, wholex2, fractionalx1, fractionalx2;
        ArrayXd vFieldx1(gridPoints), vFieldx2(gridPoints), probDensity(gridPoints); 
        ArrayXcd auxComplexVectorx1(gridPoints), auxComplexVectorx2(gridPoints), conjPsi(gridPoints);
        
        ofstream trajFile;
        trajFile.open("DATA_trajectoriesToPlot_2D_CN.txt");
        if(!trajFile.good()){
            cout << "ERROR! Unable to open the read or write files! "<<endl;
            return 1;
        }
        
        for(int tIt=1; tIt<=nt; tIt++){
        getline(readFile,line);
        for(int xIt=0; xIt<gridPoints; xIt++){
            getline(readFile,line);
            istringstream ss(line);
            ss>>arrayEl;
            currentPsi(xIt)=arrayEl;
        }        
        //compute velocity field and trajectories
        
        //Get the probability density function and the inverse psi
            
            probDensity = abs2(currentPsi);
            conjPsi= conj(currentPsi);
            
            // psi^-1* diff(psi,x)/ probDensity is computed
            
        for(int i=0; i<=nx1; i++){
            for(int j=0; j<=nx2;j++){
                if(i==0){
                    auxComplexVectorx1(j)=conjPsi(j)*(currentPsi(nx2+1+j)-currentPsi(j))/(probDensity(j)*dx1);
                } else if(i==nx1){
                    auxComplexVectorx1(nx1*(nx2+1)+j)=conjPsi(nx1*(nx2+1)+j)*(currentPsi(nx1*(nx2+1)+j)-currentPsi((nx1-1)*(nx2+1)+j))/(probDensity(nx1*(nx2+1)+j)*dx1);
                } else if(i==1 || i==(nx1-1)){
                    auxComplexVectorx1(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(currentPsi((i+1)*(nx2+1)+j)-currentPsi((i-1)*(nx2+1)+j))/(2*probDensity(i*(nx2+1)+j)*dx1);
                } else{
                    auxComplexVectorx1(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(-currentPsi((i+2)*(nx2+1)+j)+8.0*currentPsi((i+1)*(nx2+1)+j)-8.0*currentPsi((i-1)*(nx2+1)+j)+currentPsi((i-2)*(nx2+1)+j))/(12.0*probDensity(i*(nx2+1)+j)*dx1);
                }
                if(j==0){
                    auxComplexVectorx2(i*(nx2+1))=conjPsi(i*(nx2+1))*(currentPsi(i*(nx2+1)+1)-currentPsi(i*(nx2+1)))/(probDensity(i*(nx2+1))*dx2);
                } else if(j==nx2){
                    auxComplexVectorx2(i*(nx2+1)+nx2)=conjPsi(i*(nx2+1)+nx2)*(currentPsi(i*(nx2+1)+nx2)-currentPsi(i*(nx2+1)+nx2-1))/(probDensity(i*(nx2+1)+nx2)*dx2);
                } else if(j==1 || j==(nx2-1)){
                    auxComplexVectorx2(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(currentPsi(i*(nx2+1)+j+1)-currentPsi(i*(nx2+1)+j-1))/(2*probDensity(i*(nx2+1)+j)*dx2);
                }else{
                    auxComplexVectorx2(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(-currentPsi(i*(nx2+1)+j+2)+8.0*currentPsi(i*(nx2+1)+j+1)-8.0*currentPsi(i*(nx2+1)+j-1)+currentPsi(i*(nx2+1)+j-2))/(12.0*probDensity(i*(nx2+1)+j)*dx2);
                } //THERE IS SOMETHING WRONG WITH THE X1 AXIS!
            }
        }
                
            // imaginary part is extracted and Jk obtained
            
            vFieldx1 = (hbar/m1)*imag(auxComplexVectorx1);
            vFieldx2 = (hbar/m2)*imag(auxComplexVectorx2);

        
        for(int i=0; i<numTrajs;i++){
                //we apply the discretisation of the grid to the traj positions
               
                fractionalx1 = std::modf((trajectoriesx1[i]-x1min)/dx1, &wholex1);
                fractionalx2 = std::modf((trajectoriesx2[i]-x2min)/dx2, &wholex2);
               
                if(wholex1>=nx1-1){wholex1=nx1-2;}else if(wholex1<0){wholex1=0;}
                if(wholex2>=nx2-1){wholex2=nx2-2;}else if(wholex2<0){wholex2=0;}

                trajectoriesx1[i] = trajectoriesx1[i]+( (1-fractionalx1)*vFieldx1( wholex1*(nx2+1) + wholex2)+fractionalx1*vFieldx1( (wholex1+1)*(nx2+1) + wholex2 ))*dt;
                trajectoriesx2[i] = trajectoriesx2[i]+( (1-fractionalx2)*vFieldx2( wholex1*(nx2+1) + wholex2)+fractionalx2*vFieldx2( wholex1*(nx2+1) + wholex2+1 ))*dt;
              
            /*
                wholex1=round((trajectoriesx1[i]-x1min)/dx1);
                wholex2=round((trajectoriesx2[i]-x2min)/dx2);
                if(wholex1>=nx1-1){wholex1=nx1-2;}else if(wholex1<0){wholex1=0;}
                if(wholex2>=nx2-1){wholex2=nx2-2;}else if(wholex2<0){wholex2=0;}
                trajectoriesx1[i] = trajectoriesx1[i] + vFieldx1( wholex1*(nx2+1) + wholex2 )*dt;
                trajectoriesx2[i] = trajectoriesx2[i] + vFieldx2( wholex1*(nx2+1) + wholex2 )*dt;
            */    
                trajFile << trajectoriesx1[i] << " " << trajectoriesx2[i] << " 0"<< endl;
            
        }
        trajFile<<endl<<endl;
        
        absolute=abs2(currentPsi);
        
        //The norm of the surface is calculated with the 2D trapezium method
        sumCorners = 0;
        sumBorders = 0;
        sumInterior = 0;
        
        sumCorners = absolute(0) + absolute(nx2) + absolute(gridPoints-1) + absolute(nx1*(nx2+1));
        for(int ix=1; ix<nx1; ++ix){
            for(int iy=1; iy<nx2; ++iy){
                sumInterior += absolute(ix*(nx2+1)+iy);
            }
        }
        for(int iy=1; iy<nx2; ++iy){
            sumBorders += absolute(iy) + absolute(nx1*(nx2+1) + iy); // ix=0 row and ix=nx1 row
        }
        for(int ix=1; ix<nx1; ++ix){
            sumBorders += absolute(ix*(nx2+1)) + absolute(ix*(nx2+1)+nx2); //iy=0 row and iy=nx2 row
        }
        
        norm=0.25*dx1*dx2*(sumCorners + 2*sumBorders +4* sumInterior);
        writtenFile << "Norm=" << norm << endl;
        
        //the probaility density is outputted
        for(int ix=0; ix<=nx1; ix++){
            results << posx1.block(ix*nx2+ix,0, nx2+1, 1), posx2.block(ix*nx2+ix,0, nx2+1, 1), absolute.block(ix*nx2+ix,0, nx2+1, 1);
            writtenFile <<results<<endl<< endl;
        }
        
        writtenFile <<endl;
        
        }
        trajFile.close();
        
    }
    
    writtenFile.close();
    readFile.close();
    return 0;
}
