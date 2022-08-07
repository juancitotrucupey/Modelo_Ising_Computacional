#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;



string make_filename( const string& basename, double KT, const string& ext )
  {
  ostringstream result;
  result << basename << KT << ext;
  return result.str();
  }

int main(void){

ofstream datafile;
 ofstream uno;
 for(double ii=0.0;ii<3;ii+=1.0){
  datafile.open( make_filename( "body",ii, ".txt" ).c_str() );
  datafile<<ii<<"uno \n";
  datafile<<ii<<"  uno "<<ii;
  datafile.close( );
  uno.open( make_filename( "dato",ii+2, ".dat" ).c_str() );
  uno<<ii<<"bla \n";
  uno<<ii<<"  bla "<<ii;
  uno.close( );
 }
  return 0;
}
