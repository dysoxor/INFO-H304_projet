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
-n      number of results showed (10);
*/

#include "smith.h"

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
  //while it is not in end of file it read the line
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



string scoreString(int index, int score, string db, PIN* filePIN, PHR* filePHR, int maxLine){
  string res = "";
  string name = filePHR->read(filePIN, index, db);
  //We only take the first char of the name;
  if(name.size()>maxLine){
    name = name.substr(0,maxLine-1)+"...";
  }
  else {
    while(name.size()<= maxLine+3){
      name+=" ";
    }
  }
  res += (">"+name+ "     " +to_string(score));
  return res;
}

string alignementString(vector<int> result ,string query,string db, PIN* filePIN, PHR* filePHR, PSQ* filePSQ, int maxLine){
  string res = "";
  int index = result[0];
  int score = result[1];
  int startX = result[result.size()-1];
  int startY = result[result.size()-2];
  string name = filePHR->read(filePIN, index, db);
  res+=">";
  while (name.size()>maxLine){
    res+=(name.substr(0, maxLine-1)+"\n");
    name = name.substr(maxLine);
  }
  res +=(name+"\n");
  string dbSeq = filePSQ->getSequence(index);
  res += ("Longueur = " + to_string(dbSeq.size()));
  res +=("Bitscore = " + to_string(score)+"\n");
  string alignement;
  int indexX;
  int indexY;
  int x = startX;
  int y = startY;
  string line1 = "";//query
  string line2 = "";//alignement
  string line3 = "";//dbSeq
  int value;
  int sizeAlignement = result.size()-4;
  int posScore = 0;
  int gapScore = 0;
  int idScore = 0;
  for (int i = result.size()-3; i >= 2; i--){
      value = result[i];
      if (value == 1){//up : gap query
        gapScore++;
        line1+="-";
        line2+=" ";
        line3+=dbSeq[y];
        y++;
      }
      else if(value == 2){//left : gap dbSeq
        gapScore++;
        line1+=query[x];
        x++;
        line2+=" ";
        line3+= "-";
      }
      else if(value ==3){//diag : match negatif
        line1+=query[x];
        x++;
        line2+=" ";
        line3+=dbSeq[y];
        y++;
      }
      else if (value == 4){//diag : match parfait
        posScore++;
        idScore++;
        line1+=query[x];
        line2+=query[x];
        x++;
        line3+=dbSeq[y];
        y++;
      }
      else if (value == 5){//diag : match positif
        posScore++;
        line1+=query[x];
        x++;
        line2+="+";
        line3+=dbSeq[y];
        y++;
      }
  }
  double posRatio = (double)posScore/sizeAlignement;
  double gapRatio = (double) gapScore/sizeAlignement;
  double idRatio = (double)idScore/sizeAlignement;
  if ((int)idRatio > 0)
    res += ("Identities : "+to_string(idScore)+"/"+to_string(sizeAlignement)+"("+to_string((int)(idRatio*100))+"%), ");
  if ((int)posRatio > 0)
    res += ("Positives : "+to_string(posScore)+"/"+to_string(sizeAlignement)+"("+to_string((int)(posRatio*100))+"%), ");
  if ((int)gapRatio > 0)
    res += ("Gaps : "+to_string(gapScore)+"/"+to_string(sizeAlignement)+"("+to_string((int)(gapRatio*100))+"%), ");
  res+="\n";

  int endX = x;
  int endY = y;
  x = startX;
  y = startY;
  int maxLenX = to_string(endX).size();
  int maxLenY = to_string(endY).size();
  int maxLenNumber = max(maxLenX, maxLenY);
  string allLine = "";
  while (line1.size() > maxLine && line2.size() > maxLine && line3.size() > maxLine){

    //Query
    allLine += "Query :";
    for (int i = 0; i <= maxLenNumber-to_string(x).size()+2; i++){
      allLine+=" ";
    }
    allLine+=to_string(x);
    allLine+= (" " +line1.substr(0,maxLine-1) + " ");
    line1 = line1.substr(maxLine);
    x+=maxLine;
    allLine+= (to_string(x-1) + "\n");

    //alignement

    for (int i = 0; i < 10+ maxLenNumber; i++){
      allLine+=" ";
    }
    allLine+= (" "+line2.substr(0,maxLine-1)+"\n");
    line2 = line2.substr(maxLine);

    //Subject

    allLine+= "Subject :";
    for (int i = 0; i <= maxLenNumber-to_string(y).size(); i++){
      allLine+=" ";
    }
    allLine+=to_string(y);
    allLine+= (" " +line3.substr(0,maxLine-1) + " ");
    line3 = line3.substr(maxLine);
    y+=(maxLine);
    allLine+= (to_string(y-1) + "\n");


    allLine+="\n";
  }
  int restLineLenght = max(line1.size(), max(line2.size(), line3.size()));
  //query

  allLine += "Query : ";
  for (int i = 0; i <= maxLenNumber-to_string(x).size()+2; i++){
    allLine+=" ";
  }
  allLine+=to_string(x);
  allLine+= (" " +line1.substr(0,restLineLenght-1) + " ");
  x+=restLineLenght;
  allLine+= (to_string(x) + "\n");

  //alignement

  for (int i = 0; i < 10 + maxLenNumber; i++){
    allLine+=" ";
  }
  allLine+= (" "+line2.substr(0, restLineLenght-1)+"\n");

  //subject

  allLine+= "Subject : ";
  for (int i = 0; i <= maxLenNumber-to_string(y).size(); i++){
    allLine+=" ";
  }
  allLine+=to_string(y);
  allLine+= (" " +line3.substr(0,restLineLenght-1) + " ");
  y+=restLineLenght;
  allLine+= (to_string(y) + "\n");


  allLine+= "\n";
  res+=allLine;
  return res;

}
void writeOutput(vector<vector<int>> results, string outputFile, string queryFileName, string dataBaseFileName, string queryName, string querySequence, clock_t begin, PIN* filePIN, PSQ* filePSQ){
  string res = "";
  string temp = "";
  PHR* filePHR = new PHR();

  int maxLine = 60;
  res+="Name";
  string scoreTitle = "Bitscore";
  while(res.size()<=maxLine+12-scoreTitle.size()){
    res+=" ";
  }
  res+=scoreTitle;
  res+= "\n";


  for (int i = 0; i < results.size(); i++){
    temp = scoreString(results[i][0], results[i][1], dataBaseFileName, filePIN, filePHR, maxLine);
    res+=temp;
    res+="\n";
  }
  res+="\n-------------\n\n";
  for (int i = 0; i < results.size(); i++){
    temp = alignementString(results[i], querySequence, dataBaseFileName, filePIN, filePHR, filePSQ, maxLine);
    res+=temp;
    res+= "\n";
  }
  delete filePHR;
  ofstream output(outputFile);
  if (!output){
    cerr<<"Output file is not readable or accessible"<<endl;
  }
  //We write the date of execution
  time_t actualTime = time(nullptr);
  output << "FIND 1.0.0" << endl;
  output << "Authors : Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan SMEESTERS" << endl;
  output << endl;
  output << "Date : " << asctime(localtime(&actualTime)) << endl;
  output << "Database : "<< dataBaseFileName <<endl;
  output << "Number of sequences in database : " << filePIN->getNumSeq()<<endl;
  output << "Query file : " << queryFileName << endl;
  clock_t end = clock();
  double timeElapsed= double(end - begin)/CLOCKS_PER_SEC;
  output << "Elapsed time : " << timeElapsed <<"s" << endl;
  output << "Query length : " << querySequence.size() << endl;
  output << "Query full name : " << queryName << endl;
  output << "Query sequence : " << querySequence << endl;
  output<< endl << res<< endl;
  output.close();
}
int main( int argc, char **argv ){
  // It verifies if the query file and the data file are given in parameters
  if( argc > 9 ){ //Problem because we allow maximum 4 parameters (+ 4 flags)
      cerr << "Program need 4 parameters maximum" << endl;
      return EXIT_FAILURE;
  }

  //Read arguments

  string queryFileName;
  string dataFileName = "uniprot_sprot.fasta"; //default value
  string outputFile = "result.txt";
  int numberOfResults = 10;
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
    else if ((string)argv[i] == "-n"){
      numberOfResults = atoi(argv[i+1]);
    }
  }

  if(!queryGiven){
    cerr << "No query file given" << endl;
    return EXIT_FAILURE;
  }

  //Let's start the timer
  clock_t begin = clock();

  //it reads the name and the content of query sequence
  string name, content = "";
  int state;

  tie(name,content) = readFasta(queryFileName);
  if (name == "" && content == ""){
    cerr<< "Unable to read the query file (.fasta)" << endl;
    return EXIT_FAILURE;
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
  /*if (index != -1){
    // if index exists it read the info about the query sequence from *.phr
    cout <<"Reading PHR ..."<< endl;
    PHR *filePHR = new PHR();
    state= filePHR->read(filePIN,index, dataFileName);
    cout << "PHR done"<<endl;
  }*/
  if(index == EXIT_FAILURE){
    cerr << "the blast data file (.psq) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }
  /*if(state == EXIT_FAILURE){
    cerr << "the blast data file (.phr) in parameter is empty or inaccessible" << endl;
    return EXIT_FAILURE;
  }*/

  //index de P00533 = 116939
  cout << "Algorithm in process ..."<< endl;
  vector<vector<int>> results;
  results = dbAlignment(dataFileName, content, filePIN, filePSQ, numberOfResults, 116000,117000);
  cout << "Algorithm done" << endl;
  writeOutput(results, outputFile, queryFileName, dataFileName, name, content, begin,filePIN,filePSQ);




  /*cout <<"Reading PHR ..."<< endl;
  PHR *filePHR = new PHR();
  state = filePHR->read(filePIN,indexOfBestSequence, dataFileName);
  cout << "PHR done"<<endl;
  cout << "Done"<<endl;*/

  delete filePIN;
  delete filePSQ;
  //Calculate the elapsed time
  /*clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;
  cout << "Time elapsed : " << time << "s" << endl;*/

  return 0;
}
