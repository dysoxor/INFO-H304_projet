#include "smith.h"



vector<vector<int>> blosumMatrix;
map<char,int> charToInt;
vector<int> indexList;
vector<int> scoreList;
vector<vector<int>> alignementList;
int gap_op;
int gap_ex;



char conversionTable[28] ={
  '-','A','B','C','D','E','F','G','H','I',
  'K','L','M','N','P','Q','R','S','T','V',
  'W','X','Y','Z','U','*','O','J'};

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
void traceback(int maxX, int maxY, vector<vector<int>> rootAlignement){
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

int matching1(int seq1[], int index, char db[], int len1, int len2){

	int *ligne1 = new int[len2+1];
	int *ligne2 = new int[len2+1];
	int *tempLine;

  //Define variables used in order to fill the ScoringMatrix

  int maxValue = 0;
	int maxX = 0;
	int maxY = 0;
  int maxLine[len2+1];
  int maxColumn[len1+1];
  int posXmaxLine[len2+1];
  int posYmaxColumn[len1+1];
  for (int i = 0; i < len1+1; i++){
    posYmaxColumn[i] = 0;
    maxColumn[i] = 0;
  }
  for (int j = 0; j < len2+1; j++){
    posXmaxLine[j] = 0;
    maxLine[j] = 0;
		ligne1[j] = 0;
		ligne2[j] = 0;
  }
  int temp[4];
  temp[0] = 0;
	int tempIndex;

  for (int i = 1; i < len1+1; i++){
		tempLine = ligne1;
		ligne1 = ligne2;
		ligne2 = tempLine;
    for (int j = 1; j < len2+1; j++){
      temp[1] = maxColumn[i] - gap_ex*(j - posYmaxColumn[i]);
      temp[2] = maxLine[j] - gap_ex*(i - posXmaxLine[j]);
			//blosumGet = blosumMatrix[seq1[i-1]][db[index+j-1]];
			temp[3] = ligne1[j-1] + blosumMatrix[seq1[i-1]][db[index+j-1]];
			tempIndex = findMax(temp,4);
			ligne2[j] = temp[tempIndex];
      if (temp[tempIndex] >= temp[1]){
        maxColumn[i] = temp[tempIndex] - gap_op;
        posYmaxColumn[i] = j;
      }
      if (temp[tempIndex]>= temp[2]){
        maxLine[j] = temp[tempIndex] - gap_op;
        posXmaxLine[j] = i;
      }
      if (temp[tempIndex] > maxValue){
				maxValue = temp[tempIndex];
				maxX = i;
				maxY = j;
      }
    }
  }
  double lambda = 0.267;
  double logk = -3.34;
  double bitscore = double(maxValue);
  bitscore = (lambda*bitscore - logk)/log(2);
	delete ligne1;
	delete ligne2;
	return bitscore;
}

void matching2(int seq1[], int index, char db[], int len1, int len2){

	int *ligne1 = new int[len2+1];
	int *ligne2 = new int[len2+1];
	int *tempLine;

  //Define variables used in order to fill the ScoringMatrix

  int maxValue = 0;
	int maxX = 0;
	int maxY = 0;

	vector<vector<int>> rootAlignement;
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
  }
  for (int j = 0; j < len2+1; j++){
    posXmaxLine[j] = 0;
    maxLine[j] = 0;
		ligne1[j] = 0;
		ligne2[j] = 0;
  }
  int temp[4];
  temp[0] = 0;
	int tempIndex;
	int blosumGet;
	int aa1;
	int aa2;

  for (int i = 1; i < len1+1; i++){
		tempLine = ligne1;
		ligne1 = ligne2;
		ligne2 = tempLine;
    for (int j = 1; j < len2+1; j++){
      temp[1] = maxColumn[i] - gap_ex*(j - posYmaxColumn[i]);
      temp[2] = maxLine[j] - gap_ex*(i - posXmaxLine[j]);
			aa1 = seq1[i-1];
			aa2 = db[index+j-1];
			blosumGet = blosumMatrix[aa1][aa2];
			temp[3] = ligne1[j-1] + blosumGet;
			tempIndex = findMax(temp,4);
			ligne2[j] = temp[tempIndex];

			if (tempIndex == 3){ //if negative match, we keep 3
				if(aa1==aa2){//perfect match
					tempIndex=4;
				}
				else if(blosumGet>0){ //positif match but aa1 != aa2
					tempIndex=5;
				}
			}
			rootAlignement[i][j] = tempIndex;

      if (temp[tempIndex] >= temp[1]){
        maxColumn[i] = temp[tempIndex] - gap_op;
        posYmaxColumn[i] = j;
      }
      if (temp[tempIndex]>= temp[2]){
        maxLine[j] = temp[tempIndex] - gap_op;
        posXmaxLine[j] = i;
      }
      if (temp[tempIndex] > maxValue){
				maxValue = temp[tempIndex];
				maxX = i;
				maxY = j;
      }
    }
  }
  double lambda = 0.267;
  double logk = -3.34;
  double bitscore = double(maxValue);
  bitscore = (lambda*bitscore - logk)/log(2);
	traceback(maxX,maxY,rootAlignement);
	delete ligne1;
	delete ligne2;
}



vector<vector<int>> dbAlignment(string db, string Squery, PSQ* filePSQ, string smMatrix, int gpo, int gpe, int nbResults){
	gap_op = gpo;
	gap_ex = gpe;
	PIN* filePIN = filePSQ->getPIN();
	int dbSize = filePIN->getNumSeq();
	/*if (endIndex == -1){
		endIndex = dbSize;
	}*/
  clock_t begin = clock();
	setupBlosumMatrix(smMatrix);
	int len1 = Squery.size();
	map<int,char> charToInt{
		{'A',1},{'B',2},{'C',3},{'D',4},
		{'E',5},{'F',6},{'G',7},{'H',8},{'I',9},
		{'K',10},{'L',11},{'M',12},{'N',13},{'P',14},
		{'Q',15},{'R',16},{'S',17},{'T',18},{'V',19},
		{'W',20},{'X',21},{'Y',22},{'Z',23},{'U',24},
		{'*',25},{'O',26},{'J',27}
	};
	//vector<int> vquery;
	int* query = new int[len1];
	for (int i = 0; i < len1; i++){
		query[i] = charToInt[Squery.at(i)];//.push_back(conversionTable[query.at(i)]);
	}

	int pcent = 0;

	//int nbSqTraveled = endIndex-beginIndex;
	int tempScore;
	int seqOffset;
	int size;
	cout << "0% ..."<< endl;
	clock_t inter;
	double interTime;
	double estimatedTime;

  //------------------------------------------------------ SIMD PRELIM ------------------------------------------
  __attribute__((aligned (16))) int16_t score [16];
  __attribute__((aligned (16))) int16_t diagGap[16];
  __attribute__((aligned (16))) int16_t opGap[16];
  __attribute__((aligned (16))) int16_t leftScore[len1+1][16];
  __attribute__((aligned (16))) int16_t exGap[16];
  __attribute__((aligned (16))) int16_t maxScore[16];
  __attribute__((aligned (16))) int16_t zero[16];

  __attribute__((aligned (16))) int16_t maxPosLine[len1][16];
  __attribute__((aligned (16))) int16_t maxPosCol[16];
  __attribute__((aligned (16))) int16_t maxColVal[16];
  __attribute__((aligned (16))) int16_t maxLineVal[len1][16];
  __attribute__((aligned (16))) int16_t ColPos[16];
  __attribute__((aligned (16))) int16_t LinePos[16];
  __attribute__((aligned (16))) int16_t colScore[16];
  __attribute__((aligned (16))) int16_t lineScore[16];
  __attribute__((aligned (16))) int16_t unit[16];


  for(int i = 0; i < 16; i++){
    zero[i] = 0;
    unit[i] = 1;
  }


  vector<int> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  fill_n(len, 16, 0);
  int freePosition = 0;
  int residue[16];

  bool done = false;

  bool firstLine = false;

  /*int Root[16][len1+1][filePIN->getmaxSeq()+1];
  for(int i = 0; i < 16; i++){
    for(int j = 0; j < len1+1; j++)
      Root[i][j][0] = 0;
    for(int k = 0; k < filePIN->getmaxSeq() + 1; k++)
      Root[i][0][k] = 0;
  }
  int maxX[16];
  int maxY[16];*/


  //--------------------------------------- End PRELIM -------------------------------------------------
  char* seqContainer = filePSQ->getDatabase();



  double lambda = 0.267;
  double logk = -3.34;
  double bitscore;
  int index[16];

  for(int i=1; i <dbSize ; i++){
		//tempScore = matching1(query, seqOffset, filePSQ->getDatabase(), len1, size);


    //--------------------------------------------- SIMD ------------------------------------------
    if(freePosition != -1){
      index[freePosition] = i;
  		seqOffset = filePIN->getSqOffset(i);//position in .psq file of the found sequence

      analysedSequences[freePosition] = seqOffset;
      len[freePosition] = filePIN->getSqOffset(i+1)-seqOffset;
      analysedSeqIndex[freePosition] = 0;

      freePosition = -1;

      for(int t = 0; t < 16; t++){

        if(len[t] > analysedSeqIndex[t]){
          residue[t] = seqContainer[analysedSequences[t]+analysedSeqIndex[t]];
          done = false;
        }
        else{
          residue[t] = -1;
          if(freePosition == -1 && i != dbSize - 1)
            freePosition = t;
        }

      }

    }
    while(!done && (freePosition == -1 || i == dbSize-1)){
      done = true;
      for(int l = 0; l < len1; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            if(analysedSeqIndex[j] == 0){
              for(int i = 0; i < len1+1; i++){
                leftScore[i][j] = 0;
                if(i != len1){
                  maxPosLine[i][j] = 0;
                  maxLineVal[i][j] = 0;
                }
              }
              opGap[j] = 11;
              exGap[j] = 1;
              maxScore[j] = 0;
              /*maxX[j] = 0;
              maxY[j] = 0;*/
            }
            maxPosCol[j] = 0;
            maxColVal[j] = 0;
          }
          ColPos[j] = l;
          LinePos[j] = analysedSeqIndex[j];

          colScore[j] = 0;
          lineScore[j] = 0;
          if(residue[j] != -1)
            diagGap[j] = blosumMatrix[query[l]][residue[j]];
          else{
            diagGap[j] = 0;
          }
          /*if(query[l] == residue[j])
            match[j] = 4;
          else if(diagGap[j] > 0)
            match[j] = 5;
          else
            match[j] = 3;*/


        }

        //for(int w = 0; w < 2; w++){
        __attribute__((aligned (16))) int16_t mlv[16];
        __attribute__((aligned (16))) int16_t mpl[16];
        __attribute__((aligned (16))) int16_t ls[16];

        for(int c = 0; c < 16; c++){
          mlv[c] = maxLineVal[l][c];
          mpl[c] = maxPosLine[l][c];
          ls[c] = leftScore[l][c];
        }



        __m256i Z = _mm256_load_si256((__m256i*) &zero);
        __m256i RESC = _mm256_load_si256((__m256i*) &colScore);
        __m256i RESL = _mm256_load_si256((__m256i*) &lineScore);
        __m256i S = _mm256_load_si256((__m256i*) &score);

        __m256i OP = _mm256_load_si256((__m256i*) &opGap);
        __m256i COL = _mm256_load_si256((__m256i*) &maxColVal);
        __m256i PMCOL = _mm256_load_si256((__m256i*) &maxPosCol);
        __m256i PCOL = _mm256_load_si256((__m256i*) &ColPos);

        __m256i EX = _mm256_load_si256((__m256i*) &exGap);
        __m256i LIN = _mm256_load_si256((__m256i*) &mlv);
        __m256i PMLIN = _mm256_load_si256((__m256i*) &mpl);
        __m256i PLIN = _mm256_load_si256((__m256i*) &LinePos);

        __m256i D = _mm256_load_si256((__m256i*) &ls);
        __m256i P = _mm256_load_si256((__m256i*) &diagGap);

        __m256i M = _mm256_load_si256((__m256i*) &maxScore);
        __m256i UNIT = _mm256_load_si256((__m256i*) &unit);



        __m256i SUBC2 = _mm256_sub_epi16(PCOL,PMCOL);
        __m256i MULTC = _mm256_mullo_epi16(EX,SUBC2);
        RESC = _mm256_sub_epi16(COL,MULTC);

        __m256i SUBL2 = _mm256_sub_epi16(PLIN,PMLIN);
        __m256i MULTL = _mm256_mullo_epi16(EX,SUBL2);
        RESL = _mm256_sub_epi16(LIN,MULTL);

        S = _mm256_add_epi16(D,P);

        S = _mm256_max_epi16(S,RESC);
        S = _mm256_max_epi16(S,RESL);
        S = _mm256_max_epi16(S,Z);

        __m256i SUB1 = _mm256_sub_epi16(S, OP);

        __m256i CMPLTC = _mm256_add_epi16(UNIT,_mm256_cmpgt_epi16(RESC, S));
        __m256i MULTC2 = _mm256_mullo_epi16(PCOL,CMPLTC);
        __m256i NEWMPC = _mm256_max_epi16(PMCOL, MULTC2);
        __m256i NEWMC = _mm256_max_epi16(_mm256_mullo_epi16(_mm256_sub_epi16(UNIT,CMPLTC), COL),_mm256_mullo_epi16(CMPLTC, SUB1));

        __m256i CMPLTL = _mm256_add_epi16(UNIT,_mm256_cmpgt_epi16(RESL, S));
        __m256i MULTL2 = _mm256_mullo_epi16(PLIN,CMPLTL);
        __m256i NEWMPL = _mm256_max_epi16(PMLIN, MULTL2);
        __m256i NEWML = _mm256_max_epi16(_mm256_mullo_epi16(_mm256_sub_epi16(UNIT,CMPLTL), LIN),_mm256_mullo_epi16(CMPLTL, SUB1));

        _mm256_store_si256((__m256i*)&score, S);
        _mm256_store_si256((__m256i*)&ls, S);
        _mm256_store_si256((__m256i*)&maxPosCol, NEWMPC);
        _mm256_store_si256((__m256i*)&maxColVal, NEWMC);
        _mm256_store_si256((__m256i*)&mpl, NEWMPL);
        _mm256_store_si256((__m256i*)&mlv, NEWML);
        _mm256_store_si256((__m256i*)&maxScore, M = _mm256_max_epi16(S,M));
        _mm256_store_si256((__m256i*)&colScore, RESC);
        _mm256_store_si256((__m256i*)&lineScore, RESL);

        for(int c = 0; c < 16; c++){
          maxLineVal[l][c] = mlv[c];
          maxPosLine[l][c] = mpl[c];
          leftScore[l][c] = ls[c];
        }
        //}
        /*for(int x = 0; x < 16; x++){
          if(score[x] == colScore[x])
            match[x] = 1;
          else if (score[x] == lineScore[x])
            match[x] = 2;
          if(score[x] == 0)
            match[x] = 0;
          if(score[x] == maxScore[x]){
            maxX[x] = analysedSeqIndex[x];
            maxY[x] = l;
          }
          Root[x][l+1][analysedSeqIndex[x]+1] = match[x];
        }*/

        if(l == len1-1){
          for(int f = len1-1; f >= 0; f--){
            if(f == 0)
              firstLine = true;
            offset(leftScore[f+1], leftScore[f], zero, firstLine);
            firstLine = false;
          }
        }
      }

      for(int k = 0 ; k < 16; k++){
        if(residue[k] != -1)
          analysedSeqIndex[k]++;
        if(len[k] > analysedSeqIndex[k]){
          residue[k] = seqContainer[analysedSequences[k]+analysedSeqIndex[k]];
          done = false;
        }
        else{
          residue[k] = -1;
          if(freePosition == -1 && i != dbSize - 1){
            freePosition = k;
            bitscore = double(maxScore[k]);
            bitscore = (lambda*bitscore - logk)/log(2);
            scoreList.push_back(bitscore);
            indexList.push_back(index[k]);
            //traceback(maxX[k], maxY[k], Root[k]);
          }
        }
      }

    }
    //------------------------------------------- END SIMD ----------------------------------------


		if (i%(dbSize/20) == 0 && i != 0){
			inter = clock();
			//cout <<"inter " << (double)(inter-beg)/CLOCKS_PER_SEC << endl;
			pcent+=5;
		 	interTime = double(inter - begin)/CLOCKS_PER_SEC;
			estimatedTime = interTime*(100-pcent)/pcent;
			cout << pcent << "% ... (estimated time remaining : "<<(int)estimatedTime/60
      << "m"<<(int)estimatedTime%60<<"s)" << endl;
		}
  }



	merge_sort(scoreList, indexList, 0, scoreList.size()-1);
	vector<vector<int>> results;
	vector<int> tempRes;
	int compt = 0;
	for (int i = indexList.size()-1; i > indexList.size()-1 -nbResults; i--){
			tempRes.clear();
			tempRes.push_back(indexList[i]);
			tempRes.push_back(scoreList[i]);
			seqOffset = filePIN->getSqOffset(indexList[i]);
			size = filePIN->getSqOffset(indexList[i]+1)-filePIN->getSqOffset(indexList[i]);
			matching2(query,seqOffset,filePSQ->getDatabase(),len1,size);
			tempRes.insert(tempRes.end(), alignementList[compt].begin(), alignementList[compt].end());
			results.push_back(tempRes);
			compt++;
	}
	delete query;
	return results;

}

void offset(int16_t* leftScore, int16_t* diagonalScore, int16_t* zero, bool firstLine){
  for(int i = 0; i < 2; i++){

    __m128i L = _mm_load_si128((__m128i*) &leftScore[i*8]);
    __m128i D = _mm_load_si128((__m128i*) &diagonalScore[i*8]);
    _mm_store_si128((__m128i*)&leftScore[i*8], D);
    if(firstLine){
      __m128i Z = _mm_load_si128((__m128i*) &zero[i*8]);
      _mm_store_si128((__m128i*)&diagonalScore[i*8], Z);
    }
  }
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
