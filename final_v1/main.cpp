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
  /*if( argc == 5){
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

  }*/
  int state;

  //List *listProtein = readFasta(queryFileName);
  content = readFasta2(queryFileName);
  //if ( listProtein->getNumOfProtein() == 0){
  if (content == ""){
    cerr << "the query file in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  else{
    //following commented lines print the query name and sequence
    //cout << "The name is : " << listProtein->getHead()->getName() << endl;
    //cout << "The sequence is " << listProtein->getHead()->getSeqence() << endl;
  }
  if (argc == 4){
    cout << "On fait la sequence avec elle meme" << endl;
      dbAlignmentTest(content, content);


      return EXIT_SUCCESS;

  }
  // create an object PIN wich read the file *.pin
  PIN *filePIN = new PIN();
  cout << "Reading PIN ..."<< endl;
  state = filePIN->read(dataFileName);
  if(state == EXIT_FAILURE){
    cerr << "the blast data file (.pin) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  cout<<"PIN done"<<endl;
  // create an object PSQ which read the file *.psq

  cout << "Reading PSQ ..."<< endl;
  PSQ *filePSQ = new PSQ();
  // get the index of the corresponding sequence in datafile
  int index = filePSQ->read(filePIN, content, dataFileName);
  cout << "PSQ done"<<endl;
  if (index != -1){
    // if index exists it read the info about the query sequence from *.phr
    cout <<"Reading PHR ..."<< endl;
    PHR *filePHR = new PHR();
    state = filePHR->read(filePIN,index, dataFileName);
    cout << "PHR done"<<endl;
  }
  if(index == EXIT_FAILURE){
    cerr << "the blast data file (.psq) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  /*if(state == EXIT_FAILURE){
    cerr << "the blast data file (.phr) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }*/
  cout << "Algorithm in process ..."<< endl;
  int indexOfBestSequence = dbAlignment(dataFileName, content, filePIN, filePSQ);
  cout <<"Reading PHR ..."<< endl;
  PHR *filePHR = new PHR();
  state = filePHR->read(filePIN,indexOfBestSequence, dataFileName);
  cout << "PHR done"<<endl;
  cout << "Done"<<endl;

  delete filePIN;
  delete filePSQ;
  //delete listProtein;

  //int score = matching(sequence1, sequence2);

  return 0;
}
