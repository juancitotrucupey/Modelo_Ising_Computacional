#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "linear.h"
using namespace std;

const int E1 = 7;
const int AbsM1 = 8;
const int M1 = 9;
const int E2 = 10;
const int AbsM2 = 11;
const int M2 = 12;


string make_filename( const string& basename, double KT, const string& ext )
  {
  ostringstream result;
  result << basename << KT << ext;
  return result.str();
  }

void Linear_Fit(const double * Mag, int & contador, double & Coefficient, double & Slope){

  double x[contador];      double logy[contador];
  
  for(int ii=0 ;ii<contador;ii++){
    x[ii]=ii; logy[ii]=log(Mag[ii]);
  }
  
      
  Maths::Regression::Linear A(contador,x,logy);

  Slope=A.getSlope(); Coefficient=A.getCoefficient();

 if(Coefficient<0.998){
   contador = contador - max(1,(contador/10));
   if(contador < 2){contador = 2;}
 }
  
}



int main (void) {

  //Obtener datos del archivo
  int cant=800322;
  
  int ii,jj;
  

  int temp = (int)((5.0-0.2)/0.1 + 1);

  //Tabla para guardar los datos
  
  double* Datos = new double[cant*13];
  
  
  int Equ = 15000;

  int Corr;
  
  Corr=(int)((cant-Equ)/1000);

  double Correl[6][Corr];

  int Columna[6];
  Columna[0]=E1;
  Columna[1]=AbsM1;
  Columna[2]=M1;
  Columna[3]=E2;
  Columna[4]=AbsM2;
  Columna[5]=M2;
  
  double kT;
  
  double SumCorr[6]={0.0};
  double SumPro1[6]={0.0};
  double SumPro2[6]={0.0};

  double Slope[6],Coefficient[6];

  double Correlacion[6][temp*2];
  
  int contador[6];

  //Para leer y escribir en los datos
  
  ifstream myfile;
  ofstream Correlation;  
  string line;

  //Delimitadores del ciclo, temp0>=0.2, tempfi<=5.0
  double temp0 = 0.7;
  double tempfi = 0.7;
  
  for(kT=temp0;kT<=tempfi;kT+=0.1){
    
    ii=0; jj=0;
    
    temp=(int)((kT-0.2)/0.1);
    
    myfile.open(make_filename( "IsingL10_",kT, ".dat" ).c_str());
    
    // Leer los datos de los archivos .dat
        
    getline(myfile, line); //Quitar el primer dato, que corresponde a la cantidad de pasos realizados en la simulación  
    
    while (getline(myfile, line))
      {     
        // `istringstream` behaves like a normal input stream
        // but can be initialized from a string
        istringstream iss(line);
        // The input operator `>>` returns the stream
        // And streams can be used as a boolean value
        // A stream is "true" as long as everything is okay
	jj=0;
        while (iss >> Datos[ii*13+jj])
	  {
	                jj+=1;
	  }

        // Flush the standard output stream and print a newline
        ii+=1;
      }
    
    myfile.close();
    
    
    //Calcular la función de correlación magnetica
    
    Correlation.open( make_filename("Ising_Correlacion_",kT,".dat").c_str());

    for(int kk=0;kk<6;kk++){
      
    SumCorr[kk]=0.0;
    SumPro1[kk]=0.0;
    SumPro2[kk]=0.0;
    
    contador[kk]=0;
    
    for(ii=0;ii<Corr;ii++){
      
      for(jj=0;jj+Equ<cant-ii;jj++){
	SumCorr[kk]+=(Datos[(jj+Equ)*13+Columna[kk]]*Datos[(Equ+jj+ii)*13+Columna[kk]]);
	SumPro1[kk]+=Datos[(jj+Equ)*13+Columna[kk]];
	SumPro2[kk]+=Datos[(Equ+jj+ii)*13+Columna[kk]];
      }
      
      SumCorr[kk]/=(cant-(ii+Equ));
      SumPro1[kk]/=(cant-(ii+Equ));
      SumPro2[kk]/=(cant-(ii+Equ));
      
      Correl[kk][ii]=SumCorr[kk]-SumPro1[kk]*SumPro2[kk];
      
      if((Correl[kk][ii]>0)&&(contador[kk]>=ii)){contador[kk]+=1;}
      
    }
  }
    for(ii=0;ii<Corr;ii++){
      for(int kk=0;kk<6;kk++){
	if(kk==0){
	  Correlation<<ii<<" "<<Correl[kk][ii]<<" ";
	}
	else if(kk<5){Correlation<<Correl[kk][ii]<<" ";}
	else{Correlation<<Correl[kk][ii]<<" \n";}
      }
    }

    Correlation.close();
    
    
    //Realizar el ajuste exponencial a la función de correlación magnetica
    for(int kk=0;kk<6;kk++){
      
    if(contador[kk] >= 2){

      Coefficient[kk] = 0.0;
      
      while(Coefficient[kk]<0.998){
	
      Linear_Fit(Correl[kk],contador[kk],Coefficient[kk],Slope[kk]);

      }
      Correlacion[kk][2*temp]=kT; Correlacion[kk][2*temp+1]=-1/Slope[kk];
    }

    else{Correlacion[kk][2*temp]=kT; Correlacion[kk][2*temp+1]=1; }
    }
  }
  
//Se imprime el archivo con las correlaciones de los datos obtenidos
  Correlation.open("Ising_Correlacion_.dat");

  for(kT=temp0;kT<=tempfi;kT+=0.1){

    temp=(int)((kT-0.2)/0.1);
    
      for(int kk=0;kk<6;kk++){
	if(kk==0){ Correlation<<kT<<" "<<Correlacion[kk][temp*2+1]<<" ";}
	else if(kk<5){Correlation<<Correlacion[kk][temp*2+1]<<" ";}
	else{Correlation<<Correlacion[kk][temp*2+1]<<" \n";}

      }
    }

   
  Correlation.close();
  
  return 0;
}
