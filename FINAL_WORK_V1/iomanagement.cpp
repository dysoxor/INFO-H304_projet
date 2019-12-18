#include "iomanagement.h"

/*
This file is used for the input and output interactions
It writes the output file once the algorithm is done
*/

char intToChar[28] ={
  '-','A','B','C','D','E','F','G','H','I',
  'K','L','M','N','P','Q','R','S','T','V',
  'W','X','Y','Z','U','*','O','J'};


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
  //string name = filePHR->read(filePIN, index, db);
  string name = filePHR->getTitle(index);
  //We only take the first #maxLine characters of the name;
  if(name.size()>maxLine){
    name = name.substr(0,maxLine)+"...";
  }
  else {
    while(name.size()<= maxLine+3){
      name+=" ";
    }
  }
  res += (">"+name+ "     " +to_string(score));
  return res;
}

string alignementString(vector<int> result ,string query,string db, PIN* filePIN, PHR* filePHR, char dataBase[] , int maxLine){
  string res = "";
  //We extract useful information from the result vector
  //Index and score are in the 2 first position
  int index = result[0];
  int score = result[1];
  //We get the offSet in the sequence of query and subject
  int startX = result[result.size()-1];
  int startY = result[result.size()-2];
  //string name = filePHR->read(filePIN, index, db);
  string name = filePHR->getTitle(index);
  res+=">";
  //We don't want a too long line so we seperate it into lines of
  //maximum #maxLine characters
  while (name.size()>maxLine){
    res+=(name.substr(0, maxLine)+"\n");
    name = name.substr(maxLine);
  }

  res+=name;
  res+='\n';
  int seqOffset = filePIN->getSqOffset(index);//position in .psq file of the found sequence
  int size = filePIN->getSqOffset(index+1)-seqOffset;//size of the sequence's header

  //Conversion of the sequence of int from the database into a string readable for us
  string dbSeq = "";
  for (int i = 0; i <size; i++){
    dbSeq+=intToChar[dataBase[i+seqOffset]];
  }

  res += ("Length = " + to_string(dbSeq.size())+"\n");
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
      else if(value ==3){//diag : negative match
        line1+=query[x];
        x++;
        line2+=" ";
        line3+=dbSeq[y];
        y++;
      }
      else if (value == 4){//diag : perfect match
        posScore++;
        idScore++;
        line1+=query[x];
        line2+=query[x];
        x++;
        line3+=dbSeq[y];
        y++;
      }
      else if (value == 5){//diag : positive match
        posScore++;
        line1+=query[x];
        x++;
        line2+="+";
        line3+=dbSeq[y];
        y++;
      }
  }
  //We also show the ratio of positive, gap, and identity matches
  double posRatio = (double)posScore/sizeAlignement;
  double gapRatio = (double) gapScore/sizeAlignement;
  double idRatio = (double)idScore/sizeAlignement;
  if ((int)(idRatio*100) > 0)
    res += ("Identities : "+to_string(idScore)+"/"+to_string(sizeAlignement)+"("+to_string((int)(idRatio*100))+"%), ");
  if ((int)(posRatio*100) > 0)
    res += ("Positives : "+to_string(posScore)+"/"+to_string(sizeAlignement)+"("+to_string((int)(posRatio*100))+"%), ");
  if ((int)(gapRatio*100) > 0)
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
  int gapX;
  int gapY;
  string tempString;
  //While all the lines are too long, we continue to write lines of #maxLine characters
  while (line1.size() > maxLine && line2.size() > maxLine && line3.size() > maxLine){

    //Query
    allLine += "Query :";
    for (int i = 0; i <= maxLenNumber-to_string(x).size()+2; i++){
      allLine+=" ";
    }
    allLine+=to_string(x);
    tempString = line1.substr(0,maxLine);
    allLine+= (" " +tempString + " ");
    gapX = count(tempString.begin(), tempString.end(), '-');
    line1 = line1.substr(maxLine);
    x+=(maxLine-gapX);
    allLine+= (to_string(x-1) + "\n");

    //Alignement

    for (int i = 0; i < 10+ maxLenNumber; i++){
      allLine+=" ";
    }
    allLine+= (" "+line2.substr(0,maxLine)+"\n");
    line2 = line2.substr(maxLine);

    //Subject

    allLine+= "Subject :";
    for (int i = 0; i <= maxLenNumber-to_string(y).size(); i++){
      allLine+=" ";
    }
    allLine+=to_string(y);
    tempString = line3.substr(0,maxLine);
    allLine+= (" " +tempString + " ");
    gapY = count(tempString.begin(), tempString.end(), '-');

    line3 = line3.substr(maxLine);
    y+=(maxLine-gapY);
    allLine+= (to_string(y-1) + "\n");


    allLine+="\n";
  }
  //Now, we have to do it one more for the rest of the line
  int restLineLength = max(line1.size(), max(line2.size(), line3.size()));

  //Query

  allLine += "Query :";
  for (int i = 0; i <= maxLenNumber-to_string(x).size()+2; i++){
    allLine+=" ";
  }
  allLine+=to_string(x);
  tempString = line1.substr(0,restLineLength);
  allLine+= (" " +tempString + " ");
  gapX = count(tempString.begin(), tempString.end(), '-');
  x+=(restLineLength-gapX);
  allLine+= (to_string(x) + "\n");

  //Alignement

  for (int i = 0; i < 10 + maxLenNumber; i++){
    allLine+=" ";
  }
  allLine+= (" "+line2.substr(0, restLineLength)+"\n");

  //Subject

  allLine+= "Subject :";
  for (int i = 0; i <= maxLenNumber-to_string(y).size(); i++){
    allLine+=" ";
  }
  allLine+=to_string(y);
  tempString = line3.substr(0,restLineLength);
  allLine+= (" " +tempString + " ");
  gapY = count(tempString.begin(), tempString.end(), '-');

  y+=(restLineLength-gapY);
  allLine+= (to_string(y) + "\n");


  allLine+= "\n";
  res+=allLine;
  return res;

}
void writeOutput(vector<vector<int>> results, string outputFile, string queryFileName, string dataBaseFileName, string queryName, string querySequence, chrono::time_point<chrono::system_clock> begin, PIN* filePIN, PSQ* filePSQ, PHR* filePHR){
  string res = "";
  string temp = "";

  //Maximum characters per line
  int maxLine = 60;
  res+="Name";
  string scoreTitle = "Bitscore";
  while(res.size()<=maxLine+15-scoreTitle.size()){
    res+=" ";
  }
  res+=scoreTitle;
  res+= "\n";

  //We write the score
  for (int i = 0; i < results.size(); i++){
    temp = scoreString(results[i][0], results[i][1], dataBaseFileName, filePIN, filePHR, maxLine);
    res+=temp;
    res+="\n";
  }

  res+="\n-------------\n\n";

  //We write the alignements
  for (int i = 0; i < results.size(); i++){
    temp = alignementString(results[i], querySequence, dataBaseFileName, filePIN, filePHR, filePSQ->getDatabase(), maxLine);
    res+=temp;
    res+= "\n";
  }

  //We open the outputFile
  ofstream output(outputFile);
  if (!output){
    cerr<<"Output file is not readable or accessible"<<endl;
  }

  //We write the date of execution
  time_t actualTime = time(nullptr);
  output << "FIND 1.1.0" << endl;
  output << "Authors : Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan SMEESTERS" << endl;
  output << endl;
  output << "Date : " << asctime(localtime(&actualTime)) << endl;
  output << "Database : "<< dataBaseFileName <<endl;
  output << "Number of sequences in database : " << filePIN->getNumSeq()<<endl;
  output << "Query file : " << queryFileName << endl;
  chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
  auto timeElapsed = chrono::duration_cast<chrono::microseconds>(end - begin).count();
  output << "Elapsed time : " << (double)(timeElapsed/1000)/1000 <<"s" << endl;
  output << "Number of cores : " << thread::hardware_concurrency()<< endl;
  output << "Query length : " << querySequence.size() << endl;

  int tempSize;
  temp ="Query full name : ";
  //We don't want a too long line
  if(queryName.size() > maxLine-temp.size()){
    tempSize = maxLine-temp.size();
    temp+=queryName.substr(0,tempSize);
    temp+="\n";
    queryName = queryName.substr(tempSize);
  }
  while (queryName.size() > maxLine){
    temp+=queryName.substr(0,maxLine);
    temp+="\n";
    queryName = queryName.substr(maxLine);
  }
  output << temp << queryName << endl;

  temp ="Query sequence : ";
  if(querySequence.size() > maxLine-temp.size()){
    tempSize = maxLine-temp.size();
    temp+=querySequence.substr(0,tempSize);
    temp+="\n";
    querySequence = querySequence.substr(tempSize);
  }
  while (querySequence.size() > maxLine){
    temp+=querySequence.substr(0,maxLine);
    temp+="\n";
    querySequence = querySequence.substr(maxLine);
  }
  output << temp << querySequence << endl;
  output<< endl << res<< endl;
  output.close();
}
