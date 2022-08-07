#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "linear.h"
#include "Random64.h"
using namespace std;

//Columnas de los datos
const int cant=800322;
const int Equ = 10000;//con 1500 se obtiene resultados extraños para la correlación T=0.7xs
const int Corr = 100;//Tau=25
const int NumSamples = 1000; //Cantidad de datos para aplicar el método de Jacknife en el cálculo del error 
 
const int E1 = 7;
const int AbsM1 = 8;
const int M1 = 9;
const int E2 = 10;
const int AbsM2 = 11;
const int M2 = 12;
const double kB=1.0;



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

void MeanEnergy(const double *Mag, double & En1, double &dE1,double &En2, double &dE2){
  double Mean=0.0, Mean2=0.0;
  int contador=0;
  for(int ii=Equ;ii<cant;ii+=2*Corr){
    Mean+=Mag[13*ii+E1];
    Mean2+=Mag[13*ii+E1]*Mag[13*ii+E1];
    contador+=1;
  }
  Mean/=contador;
  Mean2/=contador;
  En1=Mean;
  dE1=sqrt((Mean2-Mean*Mean)/(contador-1));

  Mean=0.0, Mean2=0.0;
  contador=0;
  for(int ii=Equ;ii<cant;ii+=2*Corr){
    Mean+=Mag[13*ii+E2];
    Mean2+=Mag[13*ii+E2]*Mag[13*ii+E2];
    contador+=1;
  }
  Mean/=contador;
  Mean2/=contador;
  En2=Mean;
  dE2=sqrt((Mean2-Mean*Mean)/(contador-1));  

}

void MeanMagnetization(const double *Mag, double &Ma1, double &dM1,double &Ma2, double &dM2){

  double Mean=0.0, Mean2=0.0;
  int contador=0;
  for(int ii=Equ;ii<cant;ii+=2*Corr){
    Mean+=Mag[13*ii+AbsM1];
    Mean2+=Mag[13*ii+AbsM1]*Mag[13*ii+AbsM1];
    contador+=1;
  }
  Mean/=contador;
  Mean2/=contador;
  Ma1=Mean;
  dM1=sqrt((Mean2-Mean*Mean)/(contador-1));

  Mean=0.0, Mean2=0.0;
  contador=0;
  for(int ii=Equ;ii<cant;ii+=2*Corr){
    Mean+=Mag[13*ii+AbsM2];
    Mean2+=Mag[13*ii+AbsM2]*Mag[13*ii+AbsM2];
    contador+=1;
  }
  Mean/=contador;
  Mean2/=contador;
  Ma2=Mean;
  dM2=sqrt((Mean2-Mean*Mean)/(contador-1));

}


//Se Calcular las variables termodinamicas usado el método de Jacknife 
void CalorificCapacity(const double *Mag, double &Cc1, double &dCc1, double &Cc2, double &dCc2, double Beta){

  int Size = (int)((cant-Equ)/NumSamples);

  int contador=0;
  
  double MeanE=0.0,MeanE2=0.0,delta=0.0;

  for(int ii=Equ;ii<cant;ii+=Size){
  
    MeanE+=Mag[13*ii+E1];
    MeanE2+=Mag[13*ii+E1]*Mag[13*ii+E1];
    
   contador+=1;
  }

  MeanE/=contador;
  MeanE2/=contador;

  Cc1=Beta*Beta*kB*(MeanE2-MeanE*MeanE);

  contador=0;
  MeanE=0.0,MeanE2=0.0;

  for(int ii=Equ;ii<cant;ii+=Size){
    
    MeanE+=Mag[13*ii+E2];
    MeanE2+=Mag[13*ii+E2]*Mag[13*ii+E2];
    
    contador+=1;
  }
  MeanE/=contador;
  MeanE2/=contador;

  Cc2=Beta*Beta*kB*(MeanE2-MeanE*MeanE);

  
  
  double Datos[contador][2];
  double Muestra[contador];

  for(int ii=0;ii<contador;ii++){Datos[ii][0]=Mag[13*(Equ+(ii*Size))+E1]; Datos[ii][1]=Mag[13*(Equ+(ii*Size))+E2]; }

  //Datos Tinicial=0
  for(int ii=0;ii<contador;ii++){
    MeanE=0.0;MeanE2=0.0;
    for(int jj=0;jj<contador;jj++){
      if(jj!=ii){MeanE+=Datos[jj][0]; MeanE2 += Datos[jj][0]*Datos[jj][0];} 
    }
    MeanE/=(contador-1); MeanE2/=(contador -1);
    Muestra[ii]=Beta*Beta*kB*(MeanE2-MeanE*MeanE);
  }

  delta=0.0;
  
  for(int ii=0;ii<contador;ii++){delta+=(Muestra[ii]-Cc1)*(Muestra[ii]-Cc1);}
  delta=sqrt(delta);
  dCc1=delta;

  //Datos Tinicial=inf
  for(int ii=0;ii<contador;ii++){
    MeanE=0.0;MeanE2=0.0;
    for(int jj=0;jj<contador;jj++){
      if(jj!=ii){MeanE+=Datos[jj][1]; MeanE2 += Datos[jj][1]*Datos[jj][1];} 
    }
    MeanE/=(contador-1); MeanE2/=(contador -1);
    Muestra[ii]=Beta*Beta*kB*(MeanE2-MeanE*MeanE);
  }

  delta=0.0;
  
  for(int ii=0;ii<contador;ii++){delta+=(Muestra[ii]-Cc2)*(Muestra[ii]-Cc2);}
  delta=sqrt(delta);

  dCc2=delta;
}

void MagneticSucep(const double *Mag, double &Xs1, double &dXs1, double &Xs2, double &dXs2,double Beta){

  
  
   int Size = (int)((cant-Equ)/NumSamples);

  int contador=0;
  
  double MeanM=0.0,MeanM2=0.0,delta=0.0;

  for(int ii=Equ;ii<cant;ii+=Size){
  
    MeanM+=Mag[13*ii+AbsM1];
    MeanM2+=Mag[13*ii+AbsM1]*Mag[13*ii+AbsM1];
    
   contador+=1;
  }

  MeanM/=contador;
  MeanM2/=contador;

  Xs1=Beta*(MeanM2-MeanM*MeanM);

  contador=0;
  MeanM=0.0,MeanM2=0.0;

  for(int ii=Equ;ii<cant;ii+=Size){
    
    MeanM+=Mag[13*ii+AbsM2];
    MeanM2+=Mag[13*ii+AbsM2]*Mag[13*ii+AbsM2];
    
    contador+=1;
  }
  MeanM/=contador;
  MeanM2/=contador;

  Xs2=Beta*(MeanM2-MeanM*MeanM);

  
  
  double Datos[contador][2];
  double Muestra[contador];

  for(int ii=0;ii<contador;ii++){Datos[ii][0]=Mag[13*(Equ+(ii*Size))+AbsM1]; Datos[ii][1]=Mag[13*(Equ+(ii*Size))+AbsM2]; }

  //Datos Tinicial=0
  for(int ii=0;ii<contador;ii++){
    MeanM=0.0;MeanM2=0.0;
    for(int jj=0;jj<contador;jj++){
      if(jj!=ii){MeanM+=Datos[jj][0]; MeanM2 += Datos[jj][0]*Datos[jj][0];} 
    }
    MeanM/=(contador-1); MeanM2/=(contador -1);
    Muestra[ii]=Beta*(MeanM2-MeanM*MeanM);
  }

  delta=0.0;
  
  for(int ii=0;ii<contador;ii++){delta+=(Muestra[ii]-Xs1)*(Muestra[ii]-Xs1);}
  delta=sqrt(delta);
  dXs1=delta;

  //Datos Tinicial=inf
  for(int ii=0;ii<contador;ii++){
    MeanM=0.0;MeanM2=0.0;
    for(int jj=0;jj<contador;jj++){
      if(jj!=ii){MeanM+=Datos[jj][1]; MeanM2 += Datos[jj][1]*Datos[jj][1];} 
    }
    MeanM/=(contador-1); MeanM2/=(contador -1);
    Muestra[ii]=Beta*(MeanM2-MeanM*MeanM);
  }

  delta=0.0;
  
  for(int ii=0;ii<contador;ii++){delta+=(Muestra[ii]-Xs2)*(Muestra[ii]-Xs2);}
  delta=sqrt(delta);

  dXs2=delta;


}

void UB(const double *Mag, double &Ub1, double &dUb1, double &Ub2, double &dUb2){

  int Size = (int)((cant-Equ)/NumSamples);

  int contador=0;
  
  double MeanM2=0.0,MeanM4=0.0,delta=0.0;

  for(int ii=Equ;ii<cant;ii+=Size){
  
    MeanM2+=Mag[13*ii+AbsM1]*Mag[13*ii+AbsM1];
    MeanM4+=Mag[13*ii+AbsM1]*Mag[13*ii+AbsM1]*Mag[13*ii+AbsM1]*Mag[13*ii+AbsM1];
    
   contador+=1;
  }

  MeanM2/=contador;
  MeanM4/=contador;

  Ub1=1.0-(1.0/3.0)*MeanM4/(MeanM2*MeanM2);
  
  contador=0;
  MeanM2=0.0,MeanM4=0.0;

  for(int ii=Equ;ii<cant;ii+=Size){
    
    MeanM2+=Mag[13*ii+AbsM2]*Mag[13*ii+AbsM2];
    MeanM4+=Mag[13*ii+AbsM2]*Mag[13*ii+AbsM2]*Mag[13*ii+AbsM2]*Mag[13*ii+AbsM2];
    
    contador+=1;
  }
  MeanM2/=contador;
  MeanM4/=contador;

  Ub2=1.0-(1.0/3.0)*MeanM4/(MeanM2*MeanM2);

  
  
  double Datos[contador][2];
  double Muestra[contador];

  for(int ii=0;ii<contador;ii++){Datos[ii][0]=Mag[13*(Equ+ii*Size)+AbsM1]; Datos[ii][1]=Mag[13*(Equ+ii*Size)+AbsM2]; }

  //Datos Tinicial=0
  for(int ii=0;ii<contador;ii++){
    MeanM4=0.0;MeanM2=0.0;
    for(int jj=0;jj<contador;jj++){
      if(jj!=ii){MeanM2+=(Datos[jj][0]*Datos[jj][0]); MeanM4 += (Datos[jj][0]*Datos[jj][0]*Datos[jj][0]*Datos[jj][0]);} 
    }
    MeanM2/=(contador-1); MeanM4/=(contador -1);
    Muestra[ii]=1.0-(1.0/3.0)*MeanM4/(MeanM2*MeanM2);
  }

  delta=0.0;
  
  for(int ii=0;ii<contador;ii++){delta+=(Muestra[ii]-Ub1)*(Muestra[ii]-Ub1);}
  delta=sqrt(delta);
  dUb1=delta;

  //Datos Tinicial=inf
  for(int ii=0;ii<contador;ii++){
    MeanM4=0.0;MeanM2=0.0;
    for(int jj=0;jj<contador;jj++){
      if(jj!=ii){MeanM2+=(Datos[jj][1]*Datos[jj][1]); MeanM4 += (Datos[jj][1]*Datos[jj][1]*Datos[jj][1]*Datos[jj][1]);} 
    }
    MeanM2/=(contador-1); MeanM4/=(contador -1);
    Muestra[ii]=1.0-(1.0/3.0)*MeanM4/(MeanM2*MeanM2);
  }

  delta=0.0;
  
  for(int ii=0;ii<contador;ii++){delta+=(Muestra[ii]-Ub2)*(Muestra[ii]-Ub2);}
  delta=sqrt(delta);

  dUb2=delta;


}



int main (void) {

  //Obtener datos del archivo
   
  int ii,jj;
  
  int temp = (int)((5.0-0.2)/0.1 + 1);

  //Tabla para guardar los datos
  
  double* Datos = new double[cant*13];

  double kT;

  double E[temp][4];
  double M[temp][4];
  double Cv[temp][4];
  double Xs[temp][4];
  double Ub[temp][4];

  double Dato1, delta1, Dato2, delta2;
   //Para leer y escribir en los datos
  
  ifstream myfile;
  ofstream Caracterizacion;  
  string line;

  //Delimitadores del ciclo, temp0>=0.2, tempfi<=5.0
  double temp0 = 0.2;
  double tempfi = 5.0;
  
  for(kT=temp0;kT<=tempfi;kT+=0.1){
    
    ii=0; jj=0;
    cout<<kT<<endl;
    
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

    MeanEnergy(Datos,Dato1,delta1,Dato2,delta2);

    E[temp][0]=Dato1; E[temp][1]=delta1; E[temp][2]=Dato2; E[temp][3] = delta2;

    MeanMagnetization(Datos,Dato1,delta1,Dato2,delta2);
    
    M[temp][0]=Dato1; M[temp][1]=delta1; M[temp][2]=Dato2; M[temp][3] = delta2;
    
    CalorificCapacity(Datos,Dato1,delta1,Dato2,delta2,1.0/kT);

    Cv[temp][0]=Dato1; Cv[temp][1]=delta1; Cv[temp][2]=Dato2; Cv[temp][3] = delta2;

    MagneticSucep(Datos,Dato1,delta1,Dato2,delta2,1.0/kT);

    Xs[temp][0]=Dato1; Xs[temp][1]=delta1; Xs[temp][2]=Dato2; Xs[temp][3] = delta2;

    UB(Datos,Dato1,delta1,Dato2,delta2);

    Ub[temp][0]=Dato1; Ub[temp][1]=delta1; Ub[temp][2]=Dato2; Ub[temp][3] = delta2;

  }

  //Energia

  Caracterizacion.open("Energia.dat");
  for(kT=temp0;kT<=tempfi;kT+=0.1){
        
    temp=(int)((kT-0.2)/0.1);

    Caracterizacion<<kT<<" "<<E[temp][0]<<" "<<E[temp][1]<<" "<<E[temp][2]<<" "<<E[temp][3]<<" \n";
  }
  Caracterizacion.close();

  
      //Magnetizacion
  Caracterizacion.open("Magnetizacion.dat");
  for(kT=temp0;kT<=tempfi;kT+=0.1){
        
    temp=(int)((kT-0.2)/0.1);

    Caracterizacion<<kT<<" "<<M[temp][0]<<" "<<M[temp][1]<<" "<<M[temp][2]<<" "<<M[temp][3]<<" \n";
  }
  Caracterizacion.close();

  //Capacidad Calorifica

  Caracterizacion.open("CapacidadCalorifica.dat");
  for(kT=temp0;kT<=tempfi;kT+=0.1){
        
    temp=(int)((kT-0.2)/0.1);

    Caracterizacion<<kT<<" "<<Cv[temp][0]<<" "<<Cv[temp][1]<<" "<<Cv[temp][2]<<" "<<Cv[temp][3]<<" \n";
  }
  Caracterizacion.close();
  
      
      //Suceptibilidad

   Caracterizacion.open("SuceptibilidadMagnetica.dat");
  for(kT=temp0;kT<=tempfi;kT+=0.1){
        
    temp=(int)((kT-0.2)/0.1);

    Caracterizacion<<kT<<" "<<Xs[temp][0]<<" "<<Xs[temp][1]<<" "<<Xs[temp][2]<<" "<<Xs[temp][3]<<" \n";
  }
  Caracterizacion.close();
  

      //Ub

   Caracterizacion.open("Ub.dat");
  for(kT=temp0;kT<=tempfi;kT+=0.1){
        
    temp=(int)((kT-0.2)/0.1);

    Caracterizacion<<kT<<" "<<Ub[temp][0]<<" "<<Ub[temp][1]<<" "<<Ub[temp][2]<<" "<<Ub[temp][3]<<" \n";
  }
  Caracterizacion.close();
  
      
  
return 0;
}
