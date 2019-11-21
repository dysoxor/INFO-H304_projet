#include "smith.h"

int findMax(int tableau[], int size){
	int res = tableau[0];
	for (int i = 1; i < size; i++){
		if (tableau[i]>res){
			res = tableau[i];
		}
	}
	return res;
}

vector<vector<int>> blosumMatrix;
map<char,int> charToInt;
vector<int> indexList;
vector<int> scoreList;
vector<bool> queryAlignement;
vector<bool> dbAlignement;

void setupBlosumMatrix(string pathToBlosumMatrix){
	ifstream file(pathToBlosumMatrix);
	string line;
	bool first_line = true;
	int value;
	int n_line = 0;
	int n_column = 0;
	if (file.is_open()){
		while(getline(file,line)){ //We check every line
			if (line.at(0) != '#'){ //We don't do anything we the comments of the file
				n_column=0;
				if (first_line){
					//If we are at the first line, we can read the character of the amino acid
					first_line=false;
					for (int i = 0; i < line.length();i++){
						if (line.at(i) != ' '){
							// If we get a character different than a void char, we put it in a map
							charToInt.insert(pair<char,int>(line.at(i),n_column));
							n_column++;
						}
					}
					blosumMatrix.assign(charToInt.size(), vector<int> (charToInt.size(),0));
					//We initialize the matrix with full 0
				} else {
					line.erase(0,1);//remove the letter of the line
					for (int i = 0; i< line.length();i++){
						if (line.at(i) != ' ' && line.at(i) != '-'){//We skip the '-' character
							value = (int)line.at(i) -48; //ASCII digits starts at 48
							if (i != 0 && line.at(i-1) == '-'){
								//If there is a '-' before the int, we have a negative value
								value = -value;
							}
							blosumMatrix[n_line][n_column] = value; //Set the value
							n_column++;
						}
					}
					n_line++;
				}

			}
		}
		file.close();
	} else {
		cout << "Problem while opening the blosum file" << endl;
	}
};

int matching(vector<int> seq1, string prot2){
  //clock_t begin = clock();

	/*//creation of the blosum matrix
	BlosumMatrix* blosum = new BlosumMatrix("blosum62");*/



	//Initialization of the ScoringMatrix
	ScoringMatrix* matrix = new ScoringMatrix();


  clock_t begin = clock();
  //setupBlosumMatrix("blosum62");
  //cout <<"here" << endl;
  //define the 2 constants : gap opening and gap expansion
  const int gap_op = 11;
  const int gap_ex = 1;

  //Define the 2 sequence and place them into 2 variables prot1 and prot2

  int len1 = seq1.size();
  int len2 = prot2.size();
  /*vector<int> seq1;
  for (int i = 0; i < len1; i++){
    //seq1.push_back(blosum->charToIntConversion(prot1.at(i)));
    seq1.push_back(charToInt[prot1.at(i)]);
  }*/
  vector<int> seq2;
  for (int i = 0; i < len2; i++){
    //seq2.push_back(blosum->charToIntConversion(prot2.at(i)));
    seq2.push_back(charToInt[prot2.at(i)]);
  }
  /*string prot1 = argv[1];
  string prot2 = argv[2];*/

  //Define variables used in order to fill the ScoringMatrix
  int up = 0;
  int left = 0;
  int diag = 0;
  int maxValue = 0;
	vector<vector<int>> rootAlignement;
  vector<vector<int>> matrice;
  vector<int> tempv;
  tempv.assign(len2+1,0);
  matrice.assign(len1+1,tempv);
  int maxLine[len2+1];
  int maxColumn[len1+1];
  int posXmaxLine[len2+1];
  int posYmaxColumn[len1+1];
  for (int i = 0; i < len1+1; i++){
    posYmaxColumn[i] = 0;
    maxColumn[i] = 0;
    matrice[i][0] = 0;
  }
  for (int j = 0; j < len2+1; j++){
    posXmaxLine[j] = 0;
    maxLine[j] = 0;
    matrice[0][j] = 0;
  }
  int temp[4];
  temp[3] = 0;
  int tempVal;
  for (int i = 1; i < len1+1; i++){
    for (int j = 1; j < len2+1; j++){
      temp[0] = maxColumn[i] - gap_op - gap_ex*(j - posYmaxColumn[i]);
      temp[1] = maxLine[j] - gap_op - gap_ex*(i - posXmaxLine[j]);
      temp[2] = matrice[i-1][j-1] + blosumMatrix[seq1[i-1]][seq2[j-1]];

      //temp[2] = matrice[i-1][j-1] + blosum->get(prot1.at(i-1), prot2.at(j-1));
      tempVal = findMax(temp,4);
      matrice[i][j] = tempVal;
      if (tempVal >= maxColumn[i]-gap_ex*(j-posYmaxColumn[i])){
        maxColumn[i] = tempVal;
        posYmaxColumn[i] = j;
      }
      if (tempVal>=maxLine[j]-gap_ex*(i-posXmaxLine[j])){
        maxLine[j] = tempVal;
        posXmaxLine[j] = i;
      }
      if (tempVal > maxValue){
        maxValue = matrice[i][j];
      }
    }
  }

  clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;
  double lambda = 0.267;
  double k = 0.041;
  double bitscore = double(maxValue);
  bitscore = (lambda*bitscore - log(k))/log(2);
  //cout << "Score donne : "<< maxValue << "("<<bitscore<<" bitscore) en " << time << " secondes"<< endl;
	//free(allPos);
	//Pas oublier de delete !!!!

	//matrix->print();
	//string res = findPath(maxPos, prot1, prot2, blosum);
	//cout << "The alignement of " << argv[1] << " and " << argv[2] << " gives " << res << " and gets a score of "<< maxPos->getValue() << endl;
	return bitscore;
}

void dbAlignment(string db, string query, PIN* filePIN, PSQ* filePSQ){
  int dbSize = filePIN->getNumSeq();
  clock_t begin = clock();
	vector<int> vquery;
	for (int i = 0; i < query.size(); i++)
		vquery.push_back(charToInt[query.at(i)]);
  setupBlosumMatrix("blosum62");
	int score;
  for(int i=0; i < dbSize; i++){
		//score = matching(vquery, filePSQ->getSequence(i));
		indexList.push_back(i);
		scoreList.push_back(matching(vquery, filePSQ->getSequence(i)));

    /*if(i%1000 == 0)
      cout << "[score: " << i << "] " << score << endl;*/

  }
  clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;
  cout << "The time of matching is : " << time << endl;

}
