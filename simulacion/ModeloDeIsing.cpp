#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "Random64.h" 
using namespace std;

const int L=10;
const int L2=L*L;
const double kB=1;

// Class Spin System

class SpinSystem{

private:

  int s[L][L],E,M;
public:

  void Inicie(int ll, Crandom & ran64);
  void UnPasoDeMetropolis(double Beta,Crandom & ran64);
  double GetE(void){return (double) E;};
  double GetMR(void){return (double) M;}
  double GetM(void){return fabs(M);};
  int Getij(int ii, int jj){return s[ii][jj];}
};
void SpinSystem::Inicie(int ll,Crandom & ran64){
  int ii,jj;
  if(ll==1){//T=0 con spin up
  for( ii=0;ii<L;ii++)
    for( jj=0;jj<L;jj++)
	s[ii][jj]=1;
  E=2*L2; M=L2;
  }
  else{//Condiciones iniciales T = inf
    for( ii=0;ii<L;ii++){
      for( jj=0;jj<L;jj++){
	s[ii][jj]=1;
      }
    }

    E=2*L2; M=L2;

    int n,i,j,dE;
    
    for(ii=0;ii<2*L2;ii++){//Cambiar aleatoriamente 2*L2 spines
      n=(int)(L2*ran64.r());
      i = n/L; j=n%L;
      dE=2*s[i][j]*(s[i][(j-1+L)%L]+s[i][(j+1)%L]+s[(i-1+L)%L][j]+s[(i+1)%L][j]);
      s[i][j]*=-1; E+=dE; M+=2*s[i][j];//Se registra el cambio de las variables del sistema
     }
	
       
     }

  }
  

void SpinSystem::UnPasoDeMetropolis(double Beta,Crandom & ran64){
	
  int n,i,j,dE;
		
  //Escoger un espin al azar;
  n=(int)(L2*ran64.r());
  i=n/L;
  j=n%L;
  //Calcular el dE que se produciria si yo cambio el spin;
  dE=2*s[i][j]*(s[i][(j-1+L)%L]+s[i][(j+1)%L]+s[(i-1+L)%L][j]+s[(i+1)%L][j]);
  
  if(dE<=0)
    {s[i][j]*=-1; E+=dE; M+=2*s[i][j];} //lo volteo;
  else if(ran64.r()<exp(-Beta*dE))
    {s[i][j]*=-1; E+=dE; M+=2*s[i][j];} //lo volteo;

}

string make_filename( const string& basename, double KT, const string& ext )
  {
  ostringstream result;
  result << basename << KT << ext;
  return result.str();
  }

const int teq = (int)(200*pow(L/8.0,2.125));
const int tcorr = (int)(50*pow(L/8.0,2.125));
const int Nmuestras = 10000;

/*
const int teq = 6000;
const int tcorr = 5;
const int Nmuestras = 4;
*/
int main(void){
  
  SpinSystem Ising1;//Tinicial=0
  SpinSystem Ising2;//Tinicial=inf

  Crandom ran641(1);
  Crandom ran642(0);

  int t,mcs,cmuestras;
  
  double E1,E2,M1,M11,M2,M22;
  
  
  double kT, Beta;
  
  ofstream datafile;
  ofstream Arreglo0;
  ofstream Arregloinf;
  
  long int contador = 0;
  int CantDatos = teq + Nmuestras*tcorr +1;
  
  
  for(kT=0.2; kT<=0.2;kT+=0.1){
    
    contador = 0;
    
    datafile.open( make_filename( "IsingL10_",kT, ".dat" ).c_str() );
                  
    
    Beta = 1.0/kT;
    
    
    
    // Inicio y equilibrio el sistema 
    Ising1.Inicie(1,ran641); 
    Ising2.Inicie(0,ran642); 

    M1 = Ising1.GetM();
    M11 = Ising1.GetMR();
    E1 = Ising1.GetE();
    M2 = Ising2.GetM();
    M22 = Ising2.GetMR();
    E2 = Ising2.GetE();

    //Condiciones iniciales

    datafile<<CantDatos<<" \n";
    
    datafile<<contador<<" "<<E1<<" "<<M1<<" "<<M11<<" "<<E2<<" "<<M2<<" "<<M22<<" "<<E1/L2<<" "<<M1/L2<<" "<<M11/L2<<" "<<E2/L2<<" "<<M2/L2<<" "<<M22/L2<<" \n";
	
    for(t=0;t<teq;t++){//Se llega al equilibrio
      for(mcs=0;mcs<L2;mcs++){//Paso de Metropolis por lugar
	Ising1.UnPasoDeMetropolis(Beta,ran641);
	Ising2.UnPasoDeMetropolis(Beta,ran642);
      }
      contador +=1;
      M1 = Ising1.GetM();
      M11 = Ising1.GetMR();
      E1 = Ising1.GetE();
      M2 = Ising2.GetM();
      M22 = Ising2.GetMR();
      E2 = Ising2.GetE();
	 
      datafile<<contador<<" "<<E1<<" "<<M1<<" "<<M11<<" "<<E2<<" "<<M2<<" "<<M22<<" "<<E1/L2<<" "<<M1/L2<<" "<<M11/L2<<" "<<E2/L2<<" "<<M2/L2<<" "<<M22/L2<<" \n";

    }


    //tomo datos y evoluciono
    for(cmuestras=0;cmuestras<Nmuestras;cmuestras++){
      
      //evoluciono hasta las siguiente muestra, equivalente a tcorr pasos de Monte carlo
      for(t=0;t<tcorr;t++){
	for(mcs=0;mcs<L2;mcs++){//Un paso por lugar
	  Ising1.UnPasoDeMetropolis(Beta,ran641);
	  Ising2.UnPasoDeMetropolis(Beta,ran642);
	}
	contador+=1;
	M1 = Ising1.GetM();
	M11 = Ising1.GetMR();
	E1 = Ising1.GetE();
	M2 = Ising2.GetM();
	M22 = Ising2.GetMR();
	E2 = Ising2.GetE();
	
	
	datafile<<contador<<" "<<E1<<" "<<M1<<" "<<M11<<" "<<E2<<" "<<M2<<" "<<M22<<" "<<E1/L2<<" "<<M1/L2<<" "<<M11/L2<<" "<<E2/L2<<" "<<M2/L2<<" "<<M22/L2<<" \n";
	
      }
    }
    
    datafile.close( );
    Arreglo0.open( make_filename( "Arreglo0_IsingL10_",kT, ".dat" ).c_str() );
    for(int ii =0;ii<L;ii++){
      for(int jj=0;jj<L;jj++){
	if(jj==L-1){
	  Arreglo0<<Ising1.Getij(ii,jj)<<"\n";
	}
	else{
	  Arreglo0<<Ising1.Getij(ii,jj)<<" ";}
      }
    }
    Arreglo0.close();

    Arregloinf.open( make_filename( "Arregloinf_IsingL10_",kT, ".dat" ).c_str() );
    for(int ii =0;ii<L;ii++){
      for(int jj=0;jj<L;jj++){
	if(jj==L-1){
	  Arregloinf<<Ising2.Getij(ii,jj)<<"\n";
	}
	else{
	  Arregloinf<<Ising2.Getij(ii,jj)<<" ";}
      }
    }
    Arregloinf.close();

  }
		
	

  return 0;


}
