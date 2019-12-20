#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

using namespace std;

int main(int argc, char** argv){
  //3 arguments :
  //database_source, new_database et lettre(quelle lettre au début du nom)

  string db = argv[1];
  string new_db = argv[2];
  string first_letters = argv[3];
  int lenstr = strlen(first_letters);

  string line;
  string outputString = "";
  ifstream file(db);
  bool detected = false;
  if (file.is_open()){
    while(getline(file, line)){
      if(line[0] == '>'){
        if(compare(line.substr(1+line.find_last_of('|'), lenstr),first_letters)){
            outputString+=line;
            outputString+='\n';
            detected = true;
        }
        else {
          detected = false;
        }
      }
      else{
        if(detected){
          outputString+=line;
          outputString+='\n';
        }
      }
    }
    file.close();
  }
  ofstream outFile(new_db);
  if(outFile.is_open()){
    outFile << outputString;
    outFile.close();
  }
  return 0;


}
