#include "smith.h"



vector<vector<int>> blosumMatrix;
map<char,int> charToInt;
vector<int> indexList;
vector<int> scoreList;
vector<vector<int>> alignementList;


int findMax(int tableau[], int size){
	int res = 0;
	for (int i = 1; i < size; i++){
		if (tableau[i]>tableau[res]){
			res = i;
		}
	}
	return res;
}

int findMax(vector<int> tableau, int size){
	int res = 0;
	for (int i = 1; i < size; i++){
		if (tableau[i]>tableau[res]){
			res = i;
		}
	}
	return res;
}

void merge(vector<int> &scorev, vector<int> &indexv, int left, int mid, int right)
{
	/** Merge two sorted vectors **/
	int nl = mid-left;
	int nr = right - mid +1;
	vector<int> Lscore;
	vector<int> Lindex;
	vector<int> Rscore;
	vector<int> Rindex;
	for (int i = 0; i <= nl-1; i++){
		Lscore.push_back(scorev[left+i]);
		Lindex.push_back(indexv[left+i]);
	}
	for (int i = 0; i <= nr-1; i++){
		Rscore.push_back(scorev[mid+i]);
		Rindex.push_back(indexv[mid+i]);
	}

	const int infini = 9999999;
	Lscore.push_back(infini);
	Rscore.push_back(infini);
	int il = 0;
	int ir = 0;
	for (int i=left; i <= right; i++){
		if (Lscore[il] <= Rscore[ir]){
			scorev[i] = Lscore[il];
			indexv[i] = Lindex[il];
			il++;
		} else {
			scorev[i] = Rscore[ir];
			indexv[i] = Rindex[ir];
			ir++;
		}
	}
}

void insertion_sortmerge(vector<int> & scorev, vector<int> &indexv,int left, int right){
	int tmps;
	int tmpi;
	for(int i=left; i<right; i++)
	{
		int j=i;
		while ( j>left && scorev[j-1]>scorev[j] )
		{
			tmps = scorev[j];
			tmpi = indexv[j];
			scorev[j] = scorev[j-1];
			indexv[j] = indexv[j-1];
			scorev[j-1] = tmps;
			indexv[j-1] = tmpi;
			j--;
		}
	}
}

void merge_sort(vector<int> &scorev, vector<int> &indexv, int left, int right)
{
	/** Sort a vector by calling merge() recursively **/
	const int seuil = 256;
	if(right-left+1 < seuil) { insertion_sortmerge(scorev, indexv,left, right); }
	if (left < right){
		int mid = ceil((float)(left+right)/2);
		merge_sort(scorev, indexv, mid, right);
		merge_sort(scorev, indexv, left, mid-1);
		merge(scorev, indexv,left,mid,right);

	}
}

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
void traceback(vector<vector<int>> rootAlignement, int maxX, int maxY, int sizeX, int sizeY){
	int x = maxX;
	int y = maxY;
	vector<int> alignement;//on part de la fin
	int value = rootAlignement[x][y];
	while (value != 0){ //temp[0] = 0
		alignement.push_back(value);
		switch(value){
			case 1 	: y--;break;//gap dans query
			case 2 	: x--;break;//gap dans dbSeq
			case 3 	: x--; y--;break;//match negatif
			case 4	: x--; y--;break;//match parfait
			case 5	: x--; y--; break; //match positif
		}
		value = rootAlignement[x][y];
	}
	//offset par rapport au debut
	//d abord la db puis la query
	alignement.push_back(y); //offset de la db
	alignement.push_back(x); //offset de la query

	alignementList.push_back(alignement);
}
int matching(vector<int> seq1, string prot2, int len1){
  //clock_t begin = clock();
	//cout << seq1.size() << " "<< prot2.size() << endl;
	//int len2 = seq2.size();

	/*//creation of the blosum matrix
	BlosumMatrix* blosum = new BlosumMatrix("blosum62");*/



	//Initialization of the ScoringMatrix
	//ScoringMatrix* matrix = new ScoringMatrix();


  clock_t begin = clock();
  //setupBlosumMatrix("blosum62");
  //cout <<"here" << endl;
  //define the 2 constants : gap opening and gap expansion
  const int gap_op = 11;
  const int gap_ex = 1;

  //Define the 2 sequence and place them into 2 variables prot1 and prot2

  /*int len1 = prot1.size();
  int len2 = prot2.size();
  /*vector<int> seq1;
  for (int i = 0; i < len1; i++){
    //seq1.push_back(blosum->charToIntConversion(prot1.at(i)));
    seq1.push_back(charToInt[prot1.at(i)]);
  }*/
  vector<int> seq2;
	int len2 = prot2.size();
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
	int maxX = 0;
	int maxY = 0;
	vector<vector<int>> rootAlignement;
  vector<vector<int>> matrice;
  vector<int> tempv;
  tempv.assign(len2+1,0);
  matrice.assign(len1+1,tempv);
	vector<int> tempv2;
	tempv2.assign(len2+1,0);

	rootAlignement.assign(len1+1,tempv2);
  int maxLine[len2+1];
  int maxColumn[len1+1];
  int posXmaxLine[len2+1];
  int posYmaxColumn[len1+1];
  for (int i = 0; i < len1+1; i++){
    posYmaxColumn[i] = 0;
    maxColumn[i] = 0;
    matrice[i][0] = 0;
		rootAlignement[i][0] = 0;
  }
  for (int j = 0; j < len2+1; j++){
    posXmaxLine[j] = 0;
    maxLine[j] = 0;
    matrice[0][j] = 0;
		rootAlignement[0][j] = 0;
  }
  int temp[4];
  temp[0] = 0;
  int tempVal;
	int tempIndex;
	int blosumGet;
	int aa1;
	int aa2;
  for (int i = 1; i < len1+1; i++){
    for (int j = 1; j < len2+1; j++){
      temp[1] = maxColumn[i] - gap_op - gap_ex*(j - posYmaxColumn[i]);
      temp[2] = maxLine[j] - gap_op - gap_ex*(i - posXmaxLine[j]);
			aa1=seq1[i-1];
			aa2=seq2[j-1];
			blosumGet = blosumMatrix[aa1][aa2];
      temp[3] = matrice[i-1][j-1] + blosumGet;
			tempIndex = findMax(temp,4);

			tempVal = temp[tempIndex];
			matrice[i][j] = tempVal;
			/*if(tempVal>10)
				cout<<tempVal << " "<<tempIndex << endl;*/
			if (tempIndex == 3){ //si on a match negatif, on laisse 3
				if(aa1==aa2){//match parfait
					tempIndex=4;
				}
				else if(blosumGet>0){ //match positif mais aa1 != aa2
					tempIndex=5;
				}
			}
			rootAlignement[i][j] = tempIndex;
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
				maxX = i;
				maxY = j;
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
	if(bitscore > 2000){
		cout <<"plus de 2000"<< endl;
	}
	traceback(rootAlignement, maxX,maxY, len1+1, len2+1);
	return bitscore;
}

void dbAlignment(string db, string query, PIN* filePIN, PSQ* filePSQ){
  int dbSize = filePIN->getNumSeq();
  clock_t begin = clock();
	int len1 = query.size();
	vector<int> vquery;
	for (int i = 0; i < len1; i++){
		vquery.push_back(charToInt[query.at(i)]);
	}

  setupBlosumMatrix("blosum62");
	int score;
	int pcent = 0;
	clock_t inter;
	double interTime;
	double estimatedTime;
	cout << "0% ..."<< endl;
  for(int i=0; i <dbSize ; i++){
		//score = matching(vquery, filePSQ->getSequence(i));
		indexList.push_back(i);
		scoreList.push_back(matching(vquery, filePSQ->getSequence(i), len1));
		if (i%(dbSize/100) == 0 && i != 0){
			inter = clock();
			pcent++;
		 	interTime = double(inter - begin)/CLOCKS_PER_SEC;
			estimatedTime = interTime*(100-pcent);

			cout << pcent << "% ... (estimated time remaining : "<<estimatedTime<< " s)" << endl;

		}
  }
	merge_sort(scoreList, indexList, 0, scoreList.size()-1);

		cout << "best score" << scoreList[scoreList.size()-1] << endl;
		cout << "best score" << scoreList[findMax(scoreList, scoreList.size())] << endl;
	ofstream output("res.txt");
	for (int i = 0; i < alignementList[indexList[indexList.size()-1]].size(); i++){
		output << alignementList[indexList[indexList.size()-1]][i];
	}
	output<< endl;
	output.close();
  clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;
  cout << "The time of matching is : " << time << endl;

}
