/*
Find 1.1
This project has been made by Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan
SMEESTERS.
The purpose of this code is to align a query sequence with all the proteins
from the database. Use the Smith-Watermann Algorithm to get scores and to
get the alignements.

Usage : ./launch [OPTIONS]
Usage
-q      name of query file (required)
-d      name of data file (required)
-o      name of output file (result.txt)
-n      number of results showed (10)
-m      scoring matrix used for Smith-Waterman(blosum62)
-gpo    gap penality opening (11)
-gpe    gap penality expansion (1)
*/

#include "iomanagement.h"

int main( int argc, char **argv ){
  // It verifies if there is less than the maximum parameters allowed
  if( argc > 15 ){ //Problem because we allow maximum 7 parameters (+ 7 flags)
      cerr << "Program need 7 parameters maximum" << endl;
      return 0;
  }

  //Read arguments

  string queryFileName; //required so, no default value
  string dataFileName; //required so, no default value
  string outputFile = "result.txt"; //default value;
  int numberOfResults = 10; //default value;
  string smMatrix = "blosum62"; //default value;
  int gapPenalityOpening = 11; // default value;
  int gapPenalityExpansion = 1; // default value;

  bool queryGiven = false; //We absolutely need the query
  bool databaseGiven = false; //We also absolutely need the database

  for (int i = 1; i < argc-1; i++){
    if ((string)argv[i] == "-q"){
      queryFileName = argv[i+1];
      queryGiven = true;
    }
    else if ((string)argv[i] == "-d"){
      dataFileName = argv[i+1];
      databaseGiven = true;
    }
    else if ((string)argv[i] == "-o"){
      outputFile = argv[i+1];
    }
    else if ((string)argv[i] == "-n"){
      numberOfResults = atoi(argv[i+1]);
    }
    else if((string)argv[i] == "-m"){
      smMatrix = argv[i+1];
    }
    else if((string)argv[i] == "-gpo"){
      gapPenalityOpening = atoi(argv[i+1]);
    }
    else if((string)argv[i] == "-gpe"){
      gapPenalityExpansion = atoi(argv[i+1]);
    }
  }

  if(!queryGiven || !databaseGiven){
    cerr << "No query file given or no databaseGiven" << endl;
    return 0;
  }

  //Let's start the timer
   chrono::time_point<chrono::system_clock> begin = chrono::system_clock::now();

  //It reads the name and the content of query sequence
  string name, content = "";
  int state;

  tie(name,content) = readFasta(queryFileName);
  if (name == "" && content == ""){
    cerr<< "Unable to read the query file (.fasta)" << endl;
    return 0;
  }
  // Create an object PIN which reads the file *.pin
  PIN *filePIN = new PIN();
  state = filePIN->read(dataFileName);
  if(state == EXIT_FAILURE){
    cerr << "the blast data file (.pin) in parameter is empty or inaccessible" << endl;
    return 0;
  }
  // create an object PSQ which read the file *.psq

  PSQ *filePSQ = new PSQ();
  // We charge the content of the file in the memory
  if(filePSQ->charge(filePIN, dataFileName)){
    cout << "Error in PSQ charging" << endl;
    return 0;
  }

  //We can now start the Smith-Watermann Algorithm
  vector<vector<int>> results;
  results = dbAlignment(dataFileName, content, filePSQ, smMatrix, gapPenalityOpening,
    gapPenalityExpansion,numberOfResults);
  if(results.size() == 0){
    cout << "No results from the algorithm" << endl;
    return 0;
  }
  //Create the PHR object, which is needed to find the titles
  PHR* filePHR = new PHR();
  if(filePHR->charge(filePIN, dataFileName)){
    cout << "Error in PHR charging" << endl;
    return 0;
  }
  //Now, we write the results in the outputFile
  writeOutput(results, outputFile, queryFileName, dataFileName, name, content, begin,filePIN,filePSQ, filePHR);

  //We can now delete the different objects made on the heap
  filePSQ->end();
  filePHR->end();
  delete filePIN;
  delete filePSQ;
  delete filePHR;
  return 0;
}
