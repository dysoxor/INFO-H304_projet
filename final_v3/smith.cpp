#include "smith.h"



vector<vector<int>> blosumMatrix;
map<char,int> charToInt;
vector<int> indexList;
vector<int> scoreList;
vector<vector<int>> alignementList;
vector<vector<int>> matrix;
vector<vector<int>> rootAlignement;
int gap_op;
int gap_ex;


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
		char conversionTable[28] ={
		  '-','A','B','C','D','E','F','G','H','I',
		  'K','L','M','N','P','Q','R','S','T','V',
		  'W','X','Y','Z','U','*','O','J'};
		vector<vector<int>> otherBlosumMatrix;
		vector<int> tempv;
		tempv.assign(28,0);
		otherBlosumMatrix.assign(28,tempv);
		for (int i = 1; i < 28; i++){
			for (int j = 1; j <28; j++){
				otherBlosumMatrix[i][j] = blosumMatrix[charToInt[conversionTable[i]]][charToInt[conversionTable[j]]];
			}
		}
		blosumMatrix = otherBlosumMatrix;


	} else {
		cout << "Problem while opening the blosum file" << endl;
	}
};
void traceback(int maxX, int maxY, int sizeX, int sizeY){
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
//int matching(vector<int> seq1, vector<int> seq2, int len1){
int matching(int seq1[], int index, char db[], int len1, int len2){

	//clock_t begin = clock();
	matrix.clear();
	rootAlignement.clear();
	//clock_t bal1 = clock();

  //Define variables used in order to fill the ScoringMatrix
  int up = 0;
  int left = 0;
  int diag = 0;
  int maxValue = 0;
	int maxX = 0;
	int maxY = 0;
	//vector<vector<int>> rootAlignement;
  //vector<vector<int>> matrix;
  vector<int> tempv;
  tempv.assign(len2+1,0);
  matrix.assign(len1+1,tempv);
	vector<int> tempv2;
	tempv2.assign(len2+1,0);

	rootAlignement.assign(len1+1,tempv2);

	//clock_t bal2 = clock();

  int maxLine[len2+1];
  int maxColumn[len1+1];
  int posXmaxLine[len2+1];
  int posYmaxColumn[len1+1];
  for (int i = 0; i < len1+1; i++){
    posYmaxColumn[i] = 0;
    maxColumn[i] = 0;
    /*matrix[i][0] = 0;
		rootAlignement[i][0] = 0;*/
  }
  for (int j = 0; j < len2+1; j++){
    posXmaxLine[j] = 0;
    maxLine[j] = 0;
    /*matrix[0][j] = 0;
		rootAlignement[0][j] = 0;*/
  }
	//clock_t bal3 = clock();

  int temp[4];
  temp[0] = 0;
  int tempVal;
	int tempIndex;
	int blosumGet;
	int aa1;
	int aa2;
	//clock_t beg = clock();
	//clock_t bal4 = clock();

  for (int i = 1; i < len1+1; i++){
    for (int j = 1; j < len2+1; j++){
			//clock_t fortime1 = clock();
      temp[1] = maxColumn[i] - gap_op - gap_ex*(j - posYmaxColumn[i]);
			//clock_t fortime2 = clock();

      temp[2] = maxLine[j] - gap_op - gap_ex*(i - posXmaxLine[j]);
			//clock_t fortime3 = clock();

			aa1=seq1[i-1];
			//clock_t fortime4 = clock();

			aa2=db[index+j-1];//seq2[j-1];
			//clock_t fortime5 = clock();

			blosumGet = blosumMatrix[aa1][aa2];
			//clock_t fortime6 = clock();

      temp[3] = matrix[i-1][j-1] + blosumGet;
			//clock_t fortime7 = clock();


			tempIndex = findMax(temp,4);
			//clock_t fortime8 = clock();


			tempVal = temp[tempIndex];
			//clock_t fortime9 = clock();

			matrix[i][j] = tempVal;
			//clock_t fortime10 = clock();

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
			//clock_t fortime11 = clock();

			rootAlignement[i][j] = tempIndex;

			//clock_t fortime12 = clock();

      if (tempVal >= maxColumn[i]-gap_ex*(j-posYmaxColumn[i])){
        maxColumn[i] = tempVal;
        posYmaxColumn[i] = j;
      }
      if (tempVal>=maxLine[j]-gap_ex*(i-posXmaxLine[j])){
        maxLine[j] = tempVal;
        posXmaxLine[j] = i;
      }
      if (tempVal > maxValue){
        maxValue = matrix[i][j];
				maxX = i;
				maxY = j;
      }


    }
  }
  double lambda = 0.267;
  double logk = -3.34;
  double bitscore = double(maxValue);
  bitscore = (lambda*bitscore - logk)/log(2);
	traceback(maxX,maxY, len1+1, len2+1);
	return bitscore;
}

vector<vector<int>> dbAlignment(string db, string query, PSQ* filePSQ, string smMatrix, int gpo, int gpe, int nbResults){
	gap_op = gpo;
	gap_ex = gpe;
	PIN* filePIN = filePSQ->getPIN();
	int dbSize = filePIN->getNumSeq();
	/*if (endIndex == -1){
		endIndex = dbSize;
	}*/
  clock_t begin = clock();
	setupBlosumMatrix(smMatrix);
	int len1 = query.size();
	map<int,char> conversionTable {
		{'A',1},{'B',2},{'C',3},{'D',4},
		{'E',5},{'F',6},{'G',7},{'H',8},{'I',9},
		{'K',10},{'L',11},{'M',12},{'N',13},{'P',14},
		{'Q',15},{'R',16},{'S',17},{'T',18},{'V',19},
		{'W',20},{'X',21},{'Y',22},{'Z',23},{'U',24},
		{'*',25},{'O',26},{'J',27}
	};
	//vector<int> vquery;
	int* vquery = new int[len1];
	for (int i = 0; i < len1; i++){
		vquery[i] = conversionTable[query.at(i)];//.push_back(conversionTable[query.at(i)]);
	}

	int score;
	int pcent = 0;

	//int nbSqTraveled = endIndex-beginIndex;
	int tempScore;
	int seqOffset;
	int size;
	cout << "0% ..."<< endl;
	clock_t inter;
	double interTime;
	double estimatedTime;
  for(int i=0; i <dbSize ; i++){
		//score = matching(vquery, filePSQ->getSequence(i));
		//clock_t inter1 = clock();
		indexList.push_back(i);
		seqOffset = filePIN->getSqOffset(i);//position in .psq file of the found sequence
	  size = filePIN->getSqOffset(i+1)-seqOffset;//size of the sequence's header
		tempScore = matching(vquery, seqOffset, filePSQ->getDatabase(), len1, size);
		scoreList.push_back(tempScore);
		//inter = clock();
		//cout <<"inter " << (double)(inter-inter1)/CLOCKS_PER_SEC << endl;

		if (i%(dbSize/20) == 0 && i != 0){
			inter = clock();
			//cout <<"inter " << (double)(inter-beg)/CLOCKS_PER_SEC << endl;
			pcent+=5;
		 	interTime = double(inter - begin)/CLOCKS_PER_SEC;
			estimatedTime = interTime*(100-pcent)/pcent;
			cout << pcent << "% ... (estimated time remaining : "<<(int)estimatedTime/60<< "m"<<(int)estimatedTime%60<<"s)" << endl;
		}
  }
	merge_sort(scoreList, indexList, 0, scoreList.size()-1);
	/*ofstream output("res.txt");
	for (int i = 0; i < alignementList[indexList[indexList.size()-1]-beginIndex].size(); i++){
		output << alignementList[indexList[indexList.size()-1]-beginIndex][i];
	}
	output<< endl;
	output.close();*/
  /*clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;
  cout << "The time of matching is : " << time << endl;*/
	vector<vector<int>> results;
	vector<int> tempRes;
	for (int i = indexList.size()-1; i > indexList.size()-1 -nbResults; i--){
			tempRes.clear();
			tempRes.push_back(indexList[i]);
			tempRes.push_back(scoreList[i]);
			tempRes.insert(tempRes.end(), alignementList[indexList[i]].begin(), alignementList[indexList[i]].end());
			results.push_back(tempRes);
	}
	//filePSQ->end();
	return results;

}

/*void dbAlignmentTest(string query, string dbSeq){
  setupBlosumMatrix("blosum62");
	int len1 = query.size();
	vector<int> vquery;
	for (int i = 0; i < len1; i++){
		vquery.push_back(charToInt[query.at(i)]);
	}
	vector<int> vdbSeq;
	for (int i = 0; i < dbSeq.size(); i++){
		vdbSeq.push_back(charToInt[dbSeq.at(i)]);
	}

	int score = matching(vquery, dbSeq,len1);
	cout << "Score " << score << endl;
}*/


/*
*---------------------------------- SIMD ----------------------------------*
*/
/*
int dbAlignmentSIMD(string query, PIN* filePIN, PSQ* filePSQ){

	int dbSize = filePIN->getNumSeq();

  clock_t begin = clock();

	setupBlosumMatrix("blosum62");

	int len1 = query.size();
	vector<int> vquery;
	for (int i = 0; i < len1; i++){
		vquery.push_back(charToInt[query.at(i)]);
	}

	int score;

	int pcent = 0;
	clock_t inter;
	double interTime;
	double estimatedTime;

	vector<vector<int>> sequences = filePSQ->getAllSequences();
	filePSQ->clearSequences();
	int tempScore;
	int maxSeq = filePIN->getmaxSeq();

	__attribute__((aligned (8))) int8_t analysedSequences[16][maxSeq+1];
	int len[16];
	int analysedSeqIndex[16];
	int freePosition = 0;

	__attribute__((aligned (8))) int8_t residue [16];


	__attribute__((aligned (8))) int8_t dGap [len1];
	__attribute__((aligned (8))) int8_t lGap [len1];

	__attribute__((aligned (8))) int8_t scoresRow [len1+1];

	int maxScoreX;
	int maxScoreY;
	int maxScoreVal;


	for(int i = 1; i < dbSize; i++){
		if(freePosition != -1){
			analysedSequences[freePosition] = sequences[i];
			len[freePosition] = sequence[i].size();
			analysedSeqIndex[freePosition] = 0;
			freePosition = -1;
			for(int k = 0; k<16; k++){
				residue[k] = analysedSequences[analysedSeqIndex[k]];
				if(residue[k] == 0)
					freePosition = k;
			}
		}

		while(freePosition == -1){
			matching_SIMD(vquery, len1, residue, scoresRow, maxScoreX, maxScoreY, maxScoreVal);
			for(int k = 0; k<16; k++){
				analysedSeqIndex[k]++;
				residue[k] = analysedSequences[analysedSeqIndex[k]];
				if(residue[k] == 0)
					freePosition = k;
			}
		}
	}


}

union {
    __m128i m128;
    int8_t i8[16];
} left;
const int gap_op = 11;
const int gap_ex = 1;
void matching_SIMD(vector<int> vquery, int len1, int8_t* residue, int8_t* scoresRow,
	 									int maxScoreX,int maxScoreY, int maxScoreVal){
	__attribute__((aligned (8))) int8_t diag [len1+1];
	for(int i = 0; i < len1 ; i++){
		diag[i+1] = blosumMatrix[vquery[i]][residue];
	}
	diag[0] = 0;

	__m128i  = _mm_set1_epi8(gap_ex);
	for(int i = 1; i < len1+1; i+=16){
		left.m128 = _mm_load_si128(scoresRow[i]);
		vsum = _mm_add_epi32(vsum, v);
	}

	residuePtr.m128 = _mm_load_si128((__m128i*) residues);
	_mm_store_si128((__m128i*)sum, _mm_add_epi8(residuePtr.m128,toAddPtr.m128));
  sumPtr.m128 = _mm_load_si128((__m128i*) sum);

}
*/
/*
*--------------------------------------------------------------------------*
*/
