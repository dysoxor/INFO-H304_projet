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
		/*map<int,char> conversionTable {
			{'A',1},{'B',2},{'C',3},{'D',4},
	    {'E',5},{'F',6},{'G',7},{'H',8},{'I',9},
	    {'K',10},{'L',11},{'M',12},{'N',13},{'P',14},
	    {'Q',15},{'R',16},{'S',17},{'T',18},{'V',19},
	    {'W',20},{'X',21},{'Y',22},{'Z',23},{'U',24},
	    {'*',25},{'O',26},{'J',27}
	  };*/
		/*map<int,char> conversionTable {
	    {0,'-'},{1,'A'},{2,'B'},{3,'C'},{4,'D'},
	    {5,'E'},{6,'F'},{7,'G'},{8,'H'},{9,'I'},
	    {10,'K'},{11,'L'},{12,'M'},{13,'N'},{14,'P'},
	    {15,'Q'},{16,'R'},{17,'S'},{18,'T'},{19,'V'},
	    {20,'W'},{21,'X'},{22,'Y'},{23,'Z'},{24,'U'},
	    {25,'*'},{26,'O'},{27,'J'}
	  };*/
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
void traceback(int maxX, int maxY, vector<vector<int>> r){
	int x = maxX;
	int y = maxY;
	vector<int> alignement;//on part de la fin
	int value = r[x][y];
	while (value != 0){ //temp[0] = 0
		alignement.push_back(value);
		switch(value){
			case 1 	: y--;break;//gap dans query
			case 2 	: x--;break;//gap dans dbSeq
			case 3 	: x--; y--;break;//match negatif
			case 4	: x--; y--;break;//match parfait
			case 5	: x--; y--; break; //match positif
		}
		value = r[x][y];
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
			/*clock_t fortimeend = clock();
			double time_length = (double)(fortimeend-fortime1);
			/*double time1 = (double)(fortime2-fortime1);
			double time2 = (double)(fortime3-fortime2);
			double time3 = (double)(fortime4-fortime3);
			double time4 = (double)(fortime5-fortime4);
			double time5 = (double)(fortime6-fortime5);
			double time6 = (double)(fortime7-fortime6);
			double time7 = (double)(fortime8-fortime7);
			double time8 = (double)(fortime9-fortime8);
			double time9 = (double)(fortime10-fortime9);
			double time10 = (double)(fortime11-fortime10);
			double time11 = (double)(fortime12-fortime11);
			double time12 = (double)(fortimeend-fortime12);
			cout << "1 " << time1 << endl;
			cout << "2 " << time2 << endl;
			cout << "3 " << time3 << endl;
			cout << "4 " << time4 << endl;
			cout << "5 " << time5 << endl;
			cout << "6 " << time6 << endl;
			cout << "7 " << time7 << endl;
			cout << "8 " << time8 << endl;
			cout << "9 " << time9 << endl;
			cout << "10 " << time10 << endl;
			cout << "11 " << time11 << endl;
			cout << "12 " << time12 << endl;
			cout << "time_length " << time_length << endl;*/
			/*cout << "clocks  " << time_length << endl;
			cout << "clock estimated if 1 clock per for for " << len1*len2 << endl;
			cout << "clock estimated " << len1*len2*time_length << endl;
			while(true){

			}*/


    }
  }
	//clock_t bal5 = clock();

  /*clock_t end = clock();
  double time = double(end - begin)/CLOCKS_PER_SEC;*/
  double lambda = 0.267;
  double logk = -3.34;
  double bitscore = double(maxValue);
  bitscore = (lambda*bitscore - logk)/log(2);
	/*if(bitscore > 2000){

		cout << bitscore << " " << index << endl;
	}
	//cout << "Score donne : "<< maxValue << "("<<bitscore<<" bitscore) en " << time << " secondes"<< endl;
	//free(allPos);
	//Pas oublier de delete !!!!

	//matrix->print();
	//string res = findPath(maxPos, prot1, prot2, blosum);
	//cout << "The alignement of " << argv[1] << " and " << argv[2] << " gives " << res << " and gets a score of "<< maxPos->getValue() << endl;
	//clock_t bal6 = clock();*/


	//traceback(maxX,maxY, len1+1, len2+1);
	/*clock_t end = clock();

	int length = (int)(end-begin);
	cout << "0 " << (int)(bal1-begin) << endl;
	cout << "1 " << (int)(bal2-bal1) << endl;
	cout << "2 " << (int)(bal3-bal2) << endl;
	cout << "3 " << (int)(bal4-bal3) << endl;
	cout << "4 " << (int)(bal5-bal4) << endl;
	cout << "5 " << (int)(bal6-bal5) << endl;
	cout << "6 " << (int)(end-bal6) << endl;
	cout << "all time" << (double)(end-begin) << endl;
	while (true);*/
	return bitscore;
}

vector<vector<int>> dbAlignment(string db, string query, PSQ* filePSQ,
	string smMatrix, int gpo, int gpe, int nbResults, int beginIndex,
	int endIndex){
	gap_op = gpo;
	gap_ex = gpe;
	PIN* filePIN = filePSQ->getPIN();
	int dbSize = filePIN->getNumSeq();
	if (endIndex == -1){
		endIndex = dbSize;
	}
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


	int nbSqTraveled = endIndex-beginIndex;
	int tempScore;
	int seqOffset;
	int size;
	cout << "0% ..."<< endl;
	clock_t inter;
	double interTime;
	double estimatedTime;
	__attribute__((aligned (16))) int16_t score [16];
  __attribute__((aligned (16))) int16_t diagGap[16];
  __attribute__((aligned (16))) int16_t opGap[16];
  __attribute__((aligned (16))) int16_t HScore [16];
  __attribute__((aligned (16))) int16_t leftScore[21][16];
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
  __attribute__((aligned (16))) int16_t match[16];

	int maxX[16];
	int maxY[16];

	vector<int> initial(len1,0);


  for(int i = 0; i < 16; i++){
    zero[i] = 0;
    unit[i] = 1;
    match[i] = 0;
  }


  vector<int> analysedSequences(16);
  int analysedSeqIndex[16];
  int len[16];
  fill_n(len, 16, 0);
  int freePosition = 0;
  char residue[16];

  bool done = false;
  bool firstQuery = false;

  vector<vector<vector<int>>> tracebacks(16, vector<vector<int>>(len1, vector<int>(filePIN->getmaxSeq(),0)));
  bool WantSeq = true;
  bool firstLine = false;

  int maxS[16];

	char* seqContainer = filePSQ->getDatabase();


  for(int i = beginIndex; i < endIndex; i++){
		indexList.push_back(i);
		seqOffset = filePIN->getSqOffset(i);//position in .psq file of the found sequence
	  size = filePIN->getSqOffset(i+1)-seqOffset;//size of the sequence's header
    if(freePosition != -1){
			if(i>=16)
				traceback(maxX[freePosition], maxY[freePosition], tracebacks[freePosition]);
      analysedSequences[freePosition] = i;
      len[freePosition] = size;
      analysedSeqIndex[freePosition] = 0;
      freePosition = -1;

      for(int t = 0; t < 16; t++){

        if(len[t] > 0){
          residue[t] = seqContainer[analysedSequences[t]+analysedSeqIndex[t]];
        }
        else{
          residue[t] = '#';
          if(freePosition == -1)
            freePosition = t;
        }
      }

    }
    while(!done && (freePosition == -1 || i == endIndex-1)){
      done = true;
      if(analysedSeqIndex[0] > len[0]-1)
        WantSeq = false;
      for(int l = 0; l < len1; l++){
        for(int j = 0; j < 16; j++){
          if(l == 0){
            if(analysedSeqIndex[j] == 0){
              for(int i = 0; i < 21; i++){
                leftScore[i][j] = 0;
                if(i != 20){
                  maxPosLine[i][j] = 0;
                  maxLineVal[i][j] = 0;
                }
              }
              opGap[j] = 0;
              exGap[j] = 1;
              maxScore[j] = 0;
            }
            maxPosCol[j] = 0;
            maxColVal[j] = 0;
          }
          ColPos[j] = l;
          LinePos[j] = analysedSeqIndex[j];

          colScore[j] = 0;
          lineScore[j] = 0;
          diagGap[j] = blosumMatrix[query[l]][residue[j]];
          if(query[l] == residue[j])
            match[j] = 4;
          else if(diagGap[j] > 0)
            match[j] = 5;
          else
            match[j] = 3;


        }
        if(WantSeq && analysedSeqIndex[0] == 0){
          cout << "-------------------- AVANT[" << l << "]-------"<< query[l] <<"---" << residue[0] <<"-------------" << endl;
          cout << "score: " << score[0] << " diagap: " << diagGap[0] << " opGap: " << opGap[0]
          << " maxColVal: " << maxColVal[0] << " maxPosCol: " << maxPosCol[0] << " maxLineVal: "
          << maxLineVal[l][0] << "\nmaxPosLine: " << maxPosLine[l][0] << " exGap: " << exGap[0]
          << " leftScore: " << leftScore[l][0] << " maxScore: " << maxScore[0] << " zero: " << zero[0]
          << " ColPos: " << ColPos[0] << "\nLinePos: " << LinePos[0] << " colScore: " << colScore[0]
          << " lineScore: " << lineScore[0] << " unit: " << unit[0] << endl;
        }
        matching_SIMD(score, diagGap, opGap, maxColVal, maxPosCol, maxLineVal[l],
          maxPosLine[l], exGap, leftScore[l], maxScore, zero, ColPos, LinePos, colScore, lineScore, unit);
        for(int x = 0; x < 16; x++){
          if(score[x] == colScore[x])
            match[x] = 1;
          else if (score[x] == lineScore[x])
            match[x] = 2;
					tracebacks[x][l][analysedSeqIndex[x]] = match[x];
        }
        /*vector<int> root;
        fo(int d = 0; d < 16; d++){
          if(score[d] != diagGap[d]){
            if(score[d] == colScore[d])
              match[d] = 1;
            else
              match[d] = 2;
          }
          root
        }*/
        if(WantSeq && analysedSeqIndex[0] == 0){
          cout << "-------------------- APRES[" << l << "]-------"<< query[l] <<"---" << residue[0] <<"-------------" << endl;
          cout << "score: " << score[0] << " diagap: " << diagGap[0] << " opGap: " << opGap[0]
          << " maxColVal: " << maxColVal[0] << " maxPosCol: " << maxPosCol[0] << " maxLineVal: "
          << maxLineVal[l][0] << "\nmaxPosLine: " << maxPosLine[l][0] << " exGap: " << exGap[0]
          << " leftScore: " << leftScore[l][0] << " maxScore: " << maxScore[0] << " zero: " << zero[0]
          << " ColPos: " << ColPos[0] << "\nLinePos: " << LinePos[0] << " colScore: " << colScore[0]
          << " lineScore: " << lineScore[0] << " unit: " << unit[0] << " match: " << match[0] << endl << endl;
        }
        for(int h = 0; h < 16; h++){
					if(maxS[h] > maxScore[h]){
						maxS[h] = maxScore[h];
						maxX[h] = l;
						maxY[h] = analysedSeqIndex[h];
					}
				}
        if(l == len1-1){
          if(WantSeq && analysedSeqIndex[0] == 0)
            cout << "Décalage de LEFT" << endl;
          for(int f = 19; f >= 0; f--){
            if(f == 0)
              firstLine = true;
            offset(leftScore[f+1], leftScore[f], zero, firstLine);
            firstLine = false;
          }
          if(WantSeq && analysedSeqIndex[0] == 0){
            for(int f = 0; f < 20; f++){
              cout << "L : " << leftScore[f+1][0] << " D : " << leftScore[f][0] << endl;
            }
          }
        }
        /*if(firstResidu && firstQuery){
          //printf("[i:%d] maxScore ", i);
          for(int j = 0; j < 16; j++) {
              //printf("[j:%d]%d ", j, static_cast<int16_t>(saveScore[j]));
              printf("maxScore[%d] = %d score[%d] = %d HScore[%d] = %d upGap[%d] = %d leftScore[%d,%d] = %d leftScore[%d,%d] = %d\n", j, static_cast<int16_t>(maxScore[j]), j, static_cast<int16_t>(score[j]), j, static_cast<int16_t>(HScore[j]), j, static_cast<int16_t>(upGap[j]), l+1, j, static_cast<int16_t>(leftScore[l+1][j]), l, j, static_cast<int16_t>(leftScore[l][j]));
          }
          printf("\n");
        }*/
        //firstQuery = false;
      }


      /*if(firstResidu){
        cout << "residue : ";
        for(int m = 0; m< 16; m++){
          cout << residue[m] << " ";
        }
        cout << endl;
      }*/
      //firstResidu = false;
      for(int k = 0 ; k < 16; k++){
        analysedSeqIndex[k]++;
        if(len[k] > analysedSeqIndex[k]){
          residue[k] = seqContainer[analysedSequences[k]+analysedSeqIndex[k]];
          done = false;
        }
        else{
          residue[k] = '#';
          if(freePosition == -1 && i != endIndex - 1)
            freePosition = k;
        }
      }

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
			tempRes.insert(tempRes.end(), alignementList[indexList[i]-beginIndex].begin(), alignementList[indexList[i]-beginIndex].end());
			results.push_back(tempRes);
	}
	//filePSQ->end();
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

void matching_SIMD(int16_t* score, int16_t* diagGap, int16_t* opGap,
  int16_t* maxColVal, int16_t* maxPosCol, int16_t* maxLineVal,
  int16_t* maxPosLine, int16_t* exGap, int16_t* diagonalScore,
  int16_t* maxScore, int16_t* zero, int16_t* PosCol, int16_t* PosLine,
  int16_t* colScore, int16_t* lineScore, int16_t* unit){

  /*--------------------------------------------------------------------------
  paddsb P[q], H // H = H + P[q]
  pmaxub F, H // H = max(H, F)
  pmaxub E, H // H = max(H, E)
  pmaxub H, S // S = max(S, H)
  psubsb R, F // F = F – R
  psubsb R, E // E = E – R
  movdqa H, N // N = H
  psubsb Q, H // H = H – Q
  pmaxub H, E // E = max(H, E)
  pmaxub H, F // F = max(H, F)
  ----------------------------------------------------------------------------*/

  for(int i = 0; i < 2; i++){


    __m128i Z = _mm_load_si128((__m128i*) &zero[i*8]);
    __m128i RESC = _mm_load_si128((__m128i*) &colScore[i*8]);
    __m128i RESL = _mm_load_si128((__m128i*) &lineScore[i*8]);
    __m128i S = _mm_load_si128((__m128i*) &score[i*8]);

    __m128i OP = _mm_load_si128((__m128i*) &opGap[i*8]);
    __m128i COL = _mm_load_si128((__m128i*) &maxColVal[i*8]);
    __m128i PMCOL = _mm_load_si128((__m128i*) &maxPosCol[i*8]);
    __m128i PCOL = _mm_load_si128((__m128i*) &PosCol[i*8]);

    __m128i EX = _mm_load_si128((__m128i*) &exGap[i*8]);
    __m128i LIN = _mm_load_si128((__m128i*) &maxLineVal[i*8]);
    __m128i PMLIN = _mm_load_si128((__m128i*) &maxPosLine[i*8]);
    __m128i PLIN = _mm_load_si128((__m128i*) &PosLine[i*8]);

    __m128i D = _mm_load_si128((__m128i*) &diagonalScore[i*8]);
    __m128i P = _mm_load_si128((__m128i*) &diagGap[i*8]);

    __m128i M = _mm_load_si128((__m128i*) &maxScore[i*8]);
    __m128i UNIT = _mm_load_si128((__m128i*) &unit[i*8]);


    __m128i SUBC1 = _mm_sub_epi16(COL, OP);
    __m128i SUBC2 = _mm_sub_epi16(PCOL,PMCOL);
    __m128i MULTC = _mm_mullo_epi16(EX,SUBC2);
    RESC = _mm_sub_epi16(SUBC1,MULTC);
    __m128i SUBL1 = _mm_sub_epi16(LIN, OP);
    __m128i SUBL2 = _mm_sub_epi16(PLIN,PMLIN);
    __m128i MULTL = _mm_mullo_epi16(EX,SUBL2);
    RESL = _mm_sub_epi16(SUBL1,MULTL);

    S = _mm_add_epi16(D,P);

    S = _mm_max_epi16(S,RESC);
    S = _mm_max_epi16(S,RESL);
    S = _mm_max_epi16(S,Z);


    __m128i CMPLTC = _mm_add_epi16(UNIT,_mm_cmplt_epi16(S, COL));
    __m128i MULTC2 = _mm_mullo_epi16(PCOL,CMPLTC);
    __m128i NEWMPC = _mm_max_epi16(PMCOL, MULTC2);

    __m128i CMPLTL = _mm_add_epi16(UNIT,_mm_cmplt_epi16(S, LIN));
    __m128i MULTL2 = _mm_mullo_epi16(PLIN,CMPLTL);
    __m128i NEWMPL = _mm_max_epi16(PMLIN, MULTL2);

    _mm_store_si128((__m128i*)&score[i*8], S);
    _mm_store_si128((__m128i*)&diagonalScore[i*8], S);
    _mm_store_si128((__m128i*)&maxPosCol[i*8], NEWMPC);
    _mm_store_si128((__m128i*)&maxColVal[i*8], _mm_max_epi16(COL, S));
    _mm_store_si128((__m128i*)&maxPosLine[i*8], NEWMPL);
    _mm_store_si128((__m128i*)&maxLineVal[i*8], _mm_max_epi16(LIN, S));
    _mm_store_si128((__m128i*)&maxScore[i*8], M = _mm_max_epi16(S,M));
    _mm_store_si128((__m128i*)&colScore[i*8], RESC);
    _mm_store_si128((__m128i*)&lineScore[i*8], RESL);
  }
}
