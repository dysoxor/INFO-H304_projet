/*
This project has been made by Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan
SMEESTERS.
The purpose of this code is to find correspondance between a query sequence and
a protein from datafile which should be given in parameter such as the query.

execution in terminal: ./main -q query.fasta [-d database.fasta]
*/

#include "fasta.h"
#include "binary.h"




int main( int argc, char **argv ){
  // It verify if the query file and the data file is given in parameter
  if( argc > 5 ){ //Problem because we allow maximum 2 parameter (+2flags)
      cerr << "need 2 parameter" << endl;
      return EXIT_FAILURE;
  }

  //Read arguments

  string queryFileName;
  string dataFileName = "uniprot_sprot.fasta"; //default value
  bool queryGiven = false; //We absolutely need the query
  for (int i = 1; i < argc-1; i++){
    if ((string)argv[i] == "-d"){
      dataFileName = argv[i+1];
    }
    else if ((string)argv[i] == "-q"){
      queryFileName = argv[i+1];
      queryGiven = true;
    }
  }

  if(!queryGiven){
    cerr << "No query file given" << endl;
    return EXIT_FAILURE;
  }


  //it read the name and the content of query sequence
  string name, content = "";
  int state;

  List *listProtein = readFasta(queryFileName);
  if ( listProtein->getNumOfProtein() == 0){
    cerr << "the query file in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  //following commented lines print the query name and sequence
  //cout << "The name is : " << listProtein->getHead()->getName() << endl;
  //cout << "The sequence is " << listProtein->getHead()->getSequence() << endl;

  // create an object PIN wich read the file *.pin
  PIN *filePIN = new PIN();
  state = filePIN->read(dataFileName);
  if(state == EXIT_FAILURE){
    cerr << "the blast data file (.pin) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  // create an object PSQ which read the file *.psq
  PSQ *filePSQ = new PSQ();
  // get the index of the corresponding sequence in datafile
  int index = filePSQ->read(filePIN, listProtein->getHead()->getSeqence(), dataFileName);

  if (index != -1){
    // if index exists it read the info about the query sequence from *.phr
    PHR *filePHR = new PHR();
    state = filePHR->read(filePIN,index, dataFileName);
    delete filePHR;
  }
  if(index == EXIT_FAILURE){
    cerr << "the blast data file (.psq) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  if(state == EXIT_FAILURE){
    cerr << "the blast data file (.phr) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }


  delete filePIN;
  delete filePSQ;



  return 0;
}
