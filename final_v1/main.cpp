/*
This project has been made by Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan
SMEESTERS.
The purpose of this code is to find correspondance between a query sequence and
a protein from datafile which should be given in parameter such as the query.

execution in terminal: ./main query.fasta database.fasta
*/

#include "fasta.h"
#include "smith.h"



int main( int argc, char **argv ){
  // It verify if the query file and the data file is given in parameter
  if( argc < 3 ){
      cerr << "need 2 parameter" << endl;
      return EXIT_FAILURE;
  }
  //it read the name and the content of query sequence
  string name, content = "";
  string queryFileName = argv[1];
  string dataFileName = argv[2];
  string sequence1 = "";
  string sequence2 = "";
  if( argc == 5){
    List *seq1 = readFasta(argv[3]);
    if ( seq1->getNumOfProtein() == 0){
      sequence1 = argv[3];
    }
    else{
      sequence1 = seq1->getHead()->getSequence();
    }

    List *seq2 = readFasta(argv[4]);
    if(seq2->getNumOfProtein() == 0){
      sequence2 = argv[4];
    }
    else{
      sequence2 = seq2->getHead()->getSequence();
    }
    delete seq2;
    delete seq1;
  }
  int state;

  List *listProtein = readFasta(queryFileName);
  if ( listProtein->getNumOfProtein() == 0){
    cerr << "the query file in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  else{
    //following commented lines print the query name and sequence
    //cout << "The name is : " << listProtein->getHead()->getName() << endl;
    //cout << "The sequence is " << listProtein->getHead()->getSeqence() << endl;
  }

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
  int index = filePSQ->read(filePIN, listProtein->getHead()->getSequence(), dataFileName);

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

  dbAlignment(dataFileName, queryFileName, filePIN, filePSQ);
  delete filePIN;
  delete filePSQ;
  delete listProtein;

  //int score = matching(sequence1, sequence2);

  return 0;
}
