#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

using namespace std;

int main(int argc, char** argv){
  //3 arguments : database_source, new_database et quelle lettre
  //Il faut mettre uniprot.fasta
  string db = argv[1];
  string new_db = argv[2];
  string first_letter = argv[3];

  string line;
  string outputString = "";
  ifstream file(db);
  if (file.is_open()){

    while(getline(file, line)){
      if(line[0] == '>'){

        if(line.at(1+line.find_last_of('|')) == first_letter.at(0)){
            outputString+=line;
            cout << line << endl;
            outputString+='\n';
            while(getline(file,line)){
              if (line.at(0) == '>'){
                break;
              }
              outputString+=line;
              outputString+='\n';
            }
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
