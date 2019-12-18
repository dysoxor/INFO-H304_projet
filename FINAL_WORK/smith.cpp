#include "smith.h"

/*
This file contains all the calculations for the algorithm
This program use SIMD and multithreading if the computer allows it.
*/

//We declare global variables which are constantly used

vector<vector<int>> blosumMatrix;
map<char,int> charToInt;
vector<int> indexList;
vector<int> scoreList;
vector<vector<int>> alignementList;
int gap_op;
int gap_ex;
PIN* globalPIN;
PSQ* globalPSQ;
mutex locker;
char* seqContainer;


char conversionTable[28] ={
  '-','A','B','C','D','E','F','G','H','I',
  'K','L','M','N','P','Q','R','S','T','V',
  'W','X','Y','Z','U','*','O','J'};

struct thread_data{
  int* query;
  int len1;
  int begin;
  int end;
};

void job(struct thread_data data){
  /*----------------------------------------------------------------------------------------------------------
    score calculation using SIMD for 16 sequences in parallel
  -----------------------------------------------------------------------------------------------------------*/

    int begin = data.begin;
    int end = data.end;
    int len1 = data.len1;
    int *query = data.query;
    int dbSize = end;
    //--------------------------------INITIALIZING-------------------------------------------------------------

    __attribute__((aligned (16))) int16_t score[16]; //current score for 16 different sequences
    __attribute__((aligned (16))) int16_t diagGap[16]; //score on diagonal from the current position
    __attribute__((aligned (16))) int16_t opGap[16]; //opening gap
    __attribute__((aligned (16))) int16_t leftScore[len1+1][16]; // score on left from the current position
    __attribute__((aligned (16))) int16_t exGap[16]; //extend gap
    __attribute__((aligned (16))) int16_t maxScore[16]; //memorize the maximum score
    __attribute__((aligned (16))) int16_t zero[16]; //need a self-made matrix of 0 for the SIMD part

    __attribute__((aligned (16))) int16_t maxPosLine[len1][16]; //to memorize the position of the maximum score for each lines
    __attribute__((aligned (16))) int16_t maxPosCol[16]; //to memorize the position of the maximum score of the current column
    __attribute__((aligned (16))) int16_t maxColVal[16]; //to memorize the maximum of the column
    __attribute__((aligned (16))) int16_t maxLineVal[len1][16];//to memorize the maximum of the line
    __attribute__((aligned (16))) int16_t ColPos[16];//need current position on the column to calculate score
    __attribute__((aligned (16))) int16_t LinePos[16];//position on the line
    __attribute__((aligned (16))) int16_t colScore[16];//calculate the score from up
    __attribute__((aligned (16))) int16_t lineScore[16];//calculate the score from left
    __attribute__((aligned (16))) int16_t unit[16];//matrix of 1

    //initializing constant variables
    for(int i = 0; i < 16; i++){
      zero[i] = 0;
      unit[i] = 1;
      opGap[i] = 11;
      exGap[i] = 1;
    }


    vector<int> analysedSequences(16);//position of 16 analysed sequences in database
    int analysedSeqIndex[16];//offset of these analysed sequences
    int len[16];//size
    fill_n(len, 16, 0);
    int freePosition = 0;//if a sequence finish the calculations it is replaced by another sequence in analysedSequences
    int residue[16];//each nucleon of the sequence
  	int seqOffset;

    bool done = false;//if nothing more to analyse it becomes true

    bool firstLine = false;//managing the left and diagonal scores because they are stored in the same matrix 'leftScore'

    //normalized scores
    double lambda = 0.267;
    double logk = -3.34;
    double bitscore;

    int index[16];// store index of the sequences from database

    //------------------------------------- End initializing ----------------------------------------------------------

    for(int i=begin; i <end ; i++){

      if(freePosition != -1){//if there is a free place in analysedSequences
        index[freePosition] = i;
    		seqOffset = globalPIN->getSqOffset(i);//position in .psq file of the found sequence
        analysedSequences[freePosition] = seqOffset;
        len[freePosition] = globalPIN->getSqOffset(i+1)-seqOffset;
        analysedSeqIndex[freePosition] = 0;

        freePosition = -1;

        for(int t = 0; t < 16; t++){
          if(len[t] > analysedSeqIndex[t]){
            residue[t] = seqContainer[analysedSequences[t]+analysedSeqIndex[t]];
            done = false;//while something to analyse done = false
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

        for(int l = 0; l < len1; l++){//calculating scores for each residues with all query
          for(int j = 0; j < 16; j++){//before the SIMD need some initializes
            if(l == 0){//if it is at the first element of the query
              if(analysedSeqIndex[j] == 0){//if it is at the first element of the sequence from database
                for(int i = 0; i < len1+1; i++){
                  leftScore[i][j] = 0;
                  if(i != len1){
                    maxPosLine[i][j] = 0;
                    maxLineVal[i][j] = 0;
                  }
                }

                maxScore[j] = 0;
              }

              maxPosCol[j] = 0;
              maxColVal[j] = 0;
            }

            ColPos[j] = l;
            LinePos[j] = analysedSeqIndex[j];
            if(residue[j] != -1){
              diagGap[j] = blosumMatrix[query[l]][residue[j]];
            }
            else{
              diagGap[j] = 0;
            }
          }

          for(int w = 0; w < 2; w++){//do it 2 times because matrices are stored on 256 bits and the SSE takes 128 bits

            //--------------------------------- LOADING --------------------------------

            __m128i Z = _mm_load_si128((__m128i*) &zero[w*8]);
            __m128i RESC = _mm_load_si128((__m128i*) &colScore[w*8]);
            __m128i RESL = _mm_load_si128((__m128i*) &lineScore[w*8]);
            __m128i S = _mm_load_si128((__m128i*) &score[w*8]);

            __m128i OP = _mm_load_si128((__m128i*) &opGap[w*8]);
            __m128i COL = _mm_load_si128((__m128i*) &maxColVal[w*8]);
            __m128i PMCOL = _mm_load_si128((__m128i*) &maxPosCol[w*8]);
            __m128i PCOL = _mm_load_si128((__m128i*) &ColPos[w*8]);

            __m128i EX = _mm_load_si128((__m128i*) &exGap[w*8]);
            __m128i LIN = _mm_load_si128((__m128i*) &maxLineVal[l][w*8]);
            __m128i PMLIN = _mm_load_si128((__m128i*) &maxPosLine[l][w*8]);
            __m128i PLIN = _mm_load_si128((__m128i*) &LinePos[w*8]);

            __m128i D = _mm_load_si128((__m128i*) &leftScore[l][w*8]);
            __m128i P = _mm_load_si128((__m128i*) &diagGap[w*8]);

            __m128i M = _mm_load_si128((__m128i*) &maxScore[w*8]);
            __m128i UNIT = _mm_load_si128((__m128i*) &unit[w*8]);

            //------------------------------- calculating ---------------------------------------------


            __m128i SUBC2 = _mm_sub_epi16(PCOL,PMCOL);
            __m128i MULTC = _mm_mullo_epi16(EX,SUBC2);
            RESC = _mm_sub_epi16(COL,MULTC);// score Column = max score Column - extend Gap*(position Column - max position Column)

            __m128i SUBL2 = _mm_sub_epi16(PLIN,PMLIN);
            __m128i MULTL = _mm_mullo_epi16(EX,SUBL2);
            RESL = _mm_sub_epi16(LIN,MULTL);// score Line = max score Line - extend Gap*(position Line - max position Line)

            S = _mm_add_epi16(D,P);// score = diagonal Score + blosumMatrix[query(l)][sequence from databse(k)]

            S = _mm_max_epi16(S,RESC);// score = max(score, score Column)
            S = _mm_max_epi16(S,RESL);// score = max(score, score Line)
            S = _mm_max_epi16(S,Z);// score = max(score, 0)

            __m128i SUB1 = _mm_sub_epi16(S, OP);// sub1 = score - open Gap


            __m128i CMPLTC = _mm_add_epi16(UNIT,_mm_cmplt_epi16(SUB1, RESC));// =1 if score - open Gap >= score Column || =0 if score < score Column
            __m128i MULTC2 = _mm_mullo_epi16(PCOL,CMPLTC);
            __m128i NEWMPC = _mm_max_epi16(PMCOL, MULTC2);// new position = current column position if CMPLTC == 1 || new position does not change if CMPLTC == 0
            __m128i NEWMC = _mm_max_epi16(_mm_mullo_epi16(_mm_sub_epi16(UNIT,CMPLTC), COL),_mm_mullo_epi16(CMPLTC, SUB1));// save the new relative maximum value of the column

            //same for line
            __m128i CMPLTL = _mm_add_epi16(UNIT,_mm_cmplt_epi16(SUB1, RESL));
            __m128i MULTL2 = _mm_mullo_epi16(PLIN,CMPLTL);
            __m128i NEWMPL = _mm_max_epi16(PMLIN, MULTL2);
            __m128i NEWML = _mm_max_epi16(_mm_mullo_epi16(_mm_sub_epi16(UNIT,CMPLTL), LIN),_mm_mullo_epi16(CMPLTL, SUB1));

            //-------------------------------------- Storing the results in the arrays -------------------------------------------------------

            _mm_store_si128((__m128i*)&score[w*8], S);
            _mm_store_si128((__m128i*)&leftScore[l][w*8], S);
            _mm_store_si128((__m128i*)&maxPosCol[w*8], NEWMPC);
            _mm_store_si128((__m128i*)&maxColVal[w*8], NEWMC);
            _mm_store_si128((__m128i*)&maxPosLine[l][w*8], NEWMPL);
            _mm_store_si128((__m128i*)&maxLineVal[l][w*8], NEWML);
            _mm_store_si128((__m128i*)&maxScore[w*8], M = _mm_max_epi16(S,M));
          }


          if(l == len1-1){// if at the end of the column for the next colum the
            //leftScore has to be offseted because currently the diagonal score is the left score
            for(int f = len1-1; f >= 0; f--){
              if(f == 0)
                firstLine = true;
              offset(leftScore[f+1], leftScore[f], zero, firstLine);
              firstLine = false;
            }
          }
        }
        for(int k = 0 ; k < 16; k++){// updating the analysed sequences from database
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
              /*
              We have to use a mutex and lock this part
              because we are writing in a global variable
              */
              locker.lock();
              scoreList.push_back(bitscore);
              indexList.push_back(index[k]);
              locker.unlock();
            }
          }
        }
      }
    }
}

int findMax(int tableau[], int size){
  /*
  finding max in an array of a given size
  */
	int res = 0;
	for (int i = 1; i < size; i++){
		if (tableau[i]>tableau[res]){
			res = i;
		}
	}
	return res;
}

void merge(vector<int> &scorev, vector<int> &indexv, int left, int mid, int right){
	//Merge two sorted vectors
  //We are sorting according to score but we move the index with it
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
	for(int i=left; i<right; i++){
		int j=i;
		while ( j>left && scorev[j-1]>scorev[j] ){
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

void merge_sort(vector<int> &scorev, vector<int> &indexv, int left, int right){
	//Sort a vector with the merge_sort
  //We are sorting the two vectors according to score but we move index
  //to get them at the same position
	const int min_merge = 256;
  //If we have less than #min_merge elements, we use the insertion_sort
	if(right-left+1 < min_merge) { insertion_sortmerge(scorev, indexv,left, right); }
	if (left < right){
		int mid = ceil((float)(left+right)/2);
		merge_sort(scorev, indexv, mid, right);
		merge_sort(scorev, indexv, left, mid-1);
		merge(scorev, indexv,left,mid,right);
	}
}

int setupBlosumMatrix(string pathToBlosumMatrix){
	ifstream file(pathToBlosumMatrix);
	string line;
	bool first_line = true;
	int value;
	int n_line = 0;
	int n_column = 0;
  int first_number;
  vector<vector<int>> otherBlosumMatrix;
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
					otherBlosumMatrix.assign(charToInt.size(), vector<int> (charToInt.size(),0));
					//We initialize the matrix with full 0
				} else {
					line.erase(0,1);//remove the letter of the line
					for (int i = 0; i< line.length();i++){
						if (line.at(i) != ' ' && line.at(i) != '-'){//We skip the '-' character
            first_number = 0;
            if (i< line.length() -1 && line.at(i+1) != ' '){
              first_number = (int)line.at(i+1) -48;
              i++;
            }
            value = 10*(first_number) + (int)line.at(i) -48; //ASCII digits starts at 48
            if (i != 0 && line.at(i-1) == '-'){
								//If there is a '-' before the int, we have a negative value
								value = -value;
							}
							otherBlosumMatrix[n_line][n_column] = value; //Set the value
							n_column++;
						}
					}
					n_line++;
				}
			}
		}
		file.close();

    vector<int> tempv;
    tempv.assign(28,0);
    blosumMatrix.assign(28, tempv);
    int valueX;
    int valueY;
    int default_column;
    /*
    If we don't find a default character '*' AND
    we don't have all the letters,
    our blosumMatrix is not complete
    */
    if (charToInt.find('*') == charToInt.end() && charToInt.size() < 26){
      cout << "Problem while reading the blosumMatrix : missing scores for some amino acids combinations" << endl;
      return EXIT_FAILURE;
    }
    else {
      default_column = charToInt['*'];
    }
		for (int i = 1; i < 28; i++){
			for (int j = 1; j < 28; j++){
        //If the character is not in the table, we attribute it a default value
        if (charToInt.find(conversionTable[i]) == charToInt.end()){
          valueX = default_column;
        } else{
          valueX = charToInt[conversionTable[i]];
        }
        if (charToInt.find(conversionTable[j]) == charToInt.end()){
          valueY = default_column;
        } else{
          valueY = charToInt[conversionTable[j]];
        }
        blosumMatrix[i][j] = otherBlosumMatrix[valueX][valueY];
			}
		}
	} else {
		cout << "Problem while opening the blosum file" << endl;
    return EXIT_FAILURE;
	}
  return EXIT_SUCCESS;
}

void traceback(int maxX, int maxY, vector<vector<int>> rootAlignement){
  /*
  Building the alignment between 2 sequences
  */
	int x = maxX;
	int y = maxY;
	vector<int> alignement;//we start from the end
	int value = rootAlignement[x][y];
	while (value != 0){ //temp[0] = 0
		alignement.push_back(value);
		switch(value){
			case 1 	: y--;break;//gap in query
			case 2 	: x--;break;//gap in dbSeq
			case 3 	: x--; y--;break;//negative match
			case 4	: x--; y--;break;//perfect match
			case 5	: x--; y--; break; //positive match
		}
		value = rootAlignement[x][y];
	}
	//offset of the alignement
	//we push_back first the subject offset, then the query offset
	alignement.push_back(y); //subject offset
	alignement.push_back(x); //query offset

	alignementList.push_back(alignement);
}

void matching(int seq1[], int index, char db[], int len1, int len2){

  //We are working with two lines and not the entire matrix because
  //We don't need all the matrix to calculate a single case
	int *line1 = new int[len2+1];
	int *line2 = new int[len2+1];
	int *tempLine;

  //Define variables used in order to find the maxValue

  int maxValue = 0;
	int maxX = 0;
	int maxY = 0;

  //We set up the matrix for the traceback
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
		line1[j] = 0;
		line2[j] = 0;
  }
  int temp[4];
  temp[0] = 0;
	int tempIndex;
	int blosumGet;
	int aa1;
	int aa2;
  int possibleValue;

  for (int i = 1; i < len1+1; i++){
    //We switch line1 and line2
		tempLine = line1;
		line1 = line2;
		line2 = tempLine;
    for (int j = 1; j < len2+1; j++){
      temp[1] = maxColumn[i] - gap_ex*(j - posYmaxColumn[i]); //up
      temp[2] = maxLine[j] - gap_ex*(i - posXmaxLine[j]); //left
			aa1 = seq1[i-1];
			aa2 = db[index+j-1];
			blosumGet = blosumMatrix[aa1][aa2];
			temp[3] = line1[j-1] + blosumGet; //diag
			tempIndex = findMax(temp,4);
			line2[j] = temp[tempIndex];
      possibleValue = temp[tempIndex] -gap_op;
      if (possibleValue >= temp[1]){
        //We memoized the relative maximum of the column and its position
        maxColumn[i] = possibleValue;
        posYmaxColumn[i] = j;
      }
      if (possibleValue >= temp[2]){
        maxLine[j] = possibleValue;
        posXmaxLine[j] = i;
      }
      if (temp[tempIndex] >= maxValue){
				maxValue = temp[tempIndex];
				maxX = i;
				maxY = j;
      }
      if (tempIndex == 3){ //if negative match, we keep 3
				if(aa1==aa2){//perfect match
					tempIndex=4;
				}
				else if(blosumGet>0){ //positive match but aa1 != aa2
					tempIndex=5;
				}
			}
			rootAlignement[i][j] = tempIndex;
    }
  }
  double lambda = 0.267;
  double logk = -3.34;
  double bitscore = (lambda*(double)(maxValue) - logk)/log(2);
	traceback(maxX,maxY,rootAlignement);
	delete line1;
	delete line2;
}



vector<vector<int>> dbAlignment(string db, string Squery, PSQ* filePSQ,
  string smMatrix, int gpo, int gpe, int nbResults){

  vector<vector<int>> results;
  //This is the first function which is called in this file
  //We set up the global variables
  gap_op = gpo;
	gap_ex = gpe;
	globalPIN = filePSQ->getPIN();
  globalPSQ = filePSQ;
  PIN* filePIN = globalPIN;
	int dbSize = filePIN->getNumSeq();
  seqContainer = filePSQ->getDatabase();

  //If the nbResults is too high, we set a maximum value
  if(nbResults > dbSize){
    nbResults = dbSize;
  }

  //We set up the BlosumMatrix
	if(setupBlosumMatrix(smMatrix)){
    return results;
  }
	int len1 = Squery.size();
	map<int,char> charToInt{
		{'A',1},{'B',2},{'C',3},{'D',4},
		{'E',5},{'F',6},{'G',7},{'H',8},{'I',9},
		{'K',10},{'L',11},{'M',12},{'N',13},{'P',14},
		{'Q',15},{'R',16},{'S',17},{'T',18},{'V',19},
		{'W',20},{'X',21},{'Y',22},{'Z',23},{'U',24},
		{'*',25},{'O',26},{'J',27}
	};

  //We change the query into an int[]
  int query[len1];
	for (int i = 0; i < len1; i++){
		query[i] = charToInt[Squery.at(i)];//.push_back(conversionTable[query.at(i)]);
	}

  //--------------MULTI-THREADING----------

  //We look at the number of cores
  int n_cores = thread::hardware_concurrency();
  thread threads[n_cores];
  int seq_per_thread = dbSize/n_cores;

  struct thread_data td[n_cores];

  //If we have more than one core, we can split the work and
  //start an other thread
  for(int n = 1; n < n_cores; n++){
    td[n].begin = seq_per_thread*n;
    if (n == n_cores-1){
      td[n].end = dbSize;
    }
    else {
      td[n].end = seq_per_thread*(n+1);
    }
    td[n].query = query;
    td[n].len1 = len1;
    threads[n] = thread(job,td[n]);
  }

  //We start the main thread
  td[0].begin = 1;
  td[0].end = dbSize - (n_cores-1)*seq_per_thread;
  td[0].query = query;
  td[0].len1 = len1;
  job(td[0]);

  //We wait until the all job is done
  for (int n = 1; n < n_cores ; n++){
    threads[n].join();
  }

  //--------END MULTI-THREADING ---------

  int seqOffset;
  int size;
	merge_sort(scoreList, indexList, 0, scoreList.size()-1);
	vector<int> tempRes;
	int compt = 0;

  //The elements are sorted in ascending order
	for (int i = indexList.size()-1; i > indexList.size()-1 -nbResults; i--){
			tempRes.clear();
			tempRes.push_back(indexList[i]);
			tempRes.push_back(scoreList[i]);
			seqOffset = filePIN->getSqOffset(indexList[i]);
			size = filePIN->getSqOffset(indexList[i]+1)-seqOffset;
      //We redo the algorithm but only on the #nbResults best results to get the traceback
			matching(query,seqOffset,filePSQ->getDatabase(),len1,size);
			tempRes.insert(tempRes.end(), alignementList[compt].begin(), alignementList[compt].end());
			results.push_back(tempRes);
			compt++;
	}

	return results;
}

void offset(int16_t* leftScore, int16_t* diagonalScore, int16_t* zero, bool firstLine){
  for(int i = 0; i < 2; i++){
    /*
    Offseting 1 matrice with 2 submatrices of this leftScore and diagonalScore
    */
    __m128i L = _mm_load_si128((__m128i*) &leftScore[i*8]);
    __m128i D = _mm_load_si128((__m128i*) &diagonalScore[i*8]);
    _mm_store_si128((__m128i*)&leftScore[i*8], D);
    if(firstLine){
      __m128i Z = _mm_load_si128((__m128i*) &zero[i*8]);
      _mm_store_si128((__m128i*)&diagonalScore[i*8], Z);
    }
  }
}
