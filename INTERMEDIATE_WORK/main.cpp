/*
This project has been made by Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan
SMEESTERS.
The purpose of this code is to find correspondance between a query sequence and
a protein from datafile which should be given in parameter such as the query.

Usage : ./main [OPTIONS]
Usage
-q      name of query file (required)
-d      name of data file (uniprot_sprot.fasta)
-o      name of output file (result.txt)
          "false" value means no output file
*/

#include "binary.h"
#include <ctime>

pair<string, string> readFasta(string file){
  ifstream input(file);
  string line = "";
  string name = "";
  string content = "";
  //If it is not able to open the file it returns error and in the main it exits
  // on failure
  if(!input){
    return pair<string,string> ("","");
  }
  //while it not in end of file it read the line
  while(getline( input, line )){
      //in fasta the symbol '>' precedes the name and nexts lines are the sequences
      // of it
      if( line[0] == '>' ){
        name += line.substr(1);
      } else if (!line.empty()){
        content+=line;
      }
  }
  input.close();
  return pair<string,string> (name,content);
}

void writeOutputInfo(string outputFile, string queryFileName, string dataBaseFileName, string queryName, string querySequence, double timeElapsed){
  ofstream output(outputFile);
  if (!output){
    cerr<<"Output file is not readable or accessible"<<endl;
  }
  //We write the date of execution
  time_t actualTime = time(nullptr);
  output << "Date : " << asctime(localtime(&actualTime)) << endl;
  output << "Database : "<<dataBaseFileName<<endl;
  output << "Query file : " << queryFileName << endl;
  output << "Elapsed time : " << timeElapsed <<"s" << endl;
  output << "Query length : " << querySequence.size() << endl;
  output << "Query full name : " << queryName << endl;
  output << "Query sequence : " << querySequence << endl;
  output.close();
}



int main( int argc, char **argv ){
  // It verify if the query file and the data file is given in parameter
  if( argc > 7 ){ //Problem because we allow maximum 3 parameters (+ 3 flags)
      cerr << "Program need 3 parameters maximum" << endl;
      return EXIT_FAILURE;
  }

  //Read arguments

  string queryFileName;
  string dataFileName = "uniprot_sprot.fasta"; //default value
  string outputFile = "result.txt";
  bool outputResult = true;
  bool queryGiven = false; //We absolutely need the query
  for (int i = 1; i < argc-1; i++){
    if ((string)argv[i] == "-q"){
      queryFileName = argv[i+1];
      queryGiven = true;
    }
    else if ((string)argv[i] == "-d"){
      dataFileName = argv[i+1];
    }
    else if ((string)argv[i] == "-o"){
      if ((string)argv[i+1] == "false"){
        outputResult = false;
      }
      else {
        outputFile = argv[i+1];
      }
    }
  }

  if(!queryGiven){
    cerr << "No query file given" << endl;
    return EXIT_FAILURE;
  }

  //Let's start the timer
  clock_t begin = clock();

  //it read the name and the content of query sequence
  string name, content = "";
  int state;

  tie(name,content) = readFasta(queryFileName);
  if (name == "" && content == ""){
    cerr<< "Unable to read the query file (.fasta)" << endl;
    return EXIT_FAILURE;
  }
  //following commented lines print the query name and sequence
  //cout << "The name is : " << listProtein->getHead()->getName() << endl;
  //cout << "The sequence is " << listProtein->getHead()->getSequence() << endl;

  // create an object PIN wich read the file *.pin
  PIN *filePIN = new PIN();
  state = filePIN->read(dataFileName);
  if(state == EXIT_FAILURE){
    cerr << "Blast data file (.pin) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  // create an object PSQ which read the file *.psq
  PSQ *filePSQ = new PSQ();
  // get the index of the corresponding sequence in datafile
  int index = filePSQ->read(filePIN, content, dataFileName);

  if (index != -1){
    // if index exists it read the info about the query sequence from *.phr
    PHR *filePHR = new PHR();
    state = filePHR->read(filePIN,index, dataFileName);
    delete filePHR;
  }
  if(index == EXIT_FAILURE){
    cerr << "Blast data file (.psq) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  if(state == EXIT_FAILURE){
    cerr << "Blast data file (.phr) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }

  delete filePIN;
  delete filePSQ;

  //Calculate the elapsed time
  clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;
  cout << "Time elapsed : " << time << "s" << endl;

  //Write the results in the outPutFile
  if (outputResult){
    writeOutputInfo(outputFile, queryFileName, dataFileName, name, content, time);
  }
  return EXIT_SUCCESS;
}
