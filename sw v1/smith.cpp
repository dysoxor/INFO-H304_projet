
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "BlosumMatrix.h"
#include "ScoringMatrix.h"
#include "fasta.h"
#include <ctime>
#include <math.h>

using namespace std;


/*string findPath(Position* pos, string prot1, string prot2, BlosumMatrix* blosum){
		if (pos->getValue() == 0){ //If we fall to 0, it is the end of the alignement
			return "";
		}
		else {
			string returnValue;
			vector<Position*> rootPosVector = pos->getRoot();
			Position* rootPos = rootPosVector[0]; // we consider only the first root
			string res = findPath(rootPos, prot1, prot2, blosum);
			bool aa1 = true;
			bool aa2 = true;
		if (rootPos->getX() == pos->getX()){
			aa1=false; //if we go left, we don't take the amino acid of the first sequence
		} else if (rootPos->getY() == pos->getY()){
			aa2 = false; // same if we go up
		}
		if (res != ""){
			//If we don't get a 0, we can add the res and a / to separate
			//the local alignement between two amino acids of each sequence
			returnValue+=res;
			returnValue+= "/";
		}
		if (aa1){
			//If we go don't go left, we can add the corresponding amino acid to
			// alignement
			returnValue+=prot1[pos->getX()-1];
		} else {//Else, we put a gap
			returnValue+="-";
		}
		if (aa2){//Same if we don't go up
			returnValue+=prot2[pos->getY()-1];
		} else {
			returnValue+="-";
		}
		return returnValue;
	}
}
*/
int findMax(int tableau[], int size){
	int res = tableau[0];
	for (int i = 1; i < size; i++){
		if (tableau[i]>res){
			res = tableau[i];
		}
	}
	return res;
}

int findMax2(int t1[], int t2[], int s1, int s2){
	int res = t1[0];
	for (int i = 1; i < s1; i++){
		if (t1[i]> res){
			res = t1[i];
		}
	}
	for (int j = 0; j < s2; j++){
		if(t2[j]>res){
			res = t2[j];
		}
	}
	return res;
}

int main(int argc, char** argv){
	if (argc != 3){
		cout << "Error, not the right number of parameters" << endl;
		cout << "Number of parameters expected : 2" << endl;
		cout <<"Number of parameters given : " << (argc - 1) << endl;
		return EXIT_FAILURE;
	}

	List *listProtein1 = readFasta(argv[1]);
  string content1 = listProtein1->getHead()->getSequence();

	List *listProtein2 = readFasta(argv[2]);
  string content2 = listProtein2->getHead()->getSequence();

	//creation of the blosum matrix
	BlosumMatrix* blosum = new BlosumMatrix("blosum62");

	//Initialization of the ScoringMatrix
	//ScoringMatrix* matrix = new ScoringMatrix();

	//define the 2 constants : gap opening and gap expansion
	const int gap_op = 11;
	const int gap_ex = 1;

	//Define the 2 sequence and place them into 2 variables prot1 and prot2

	string prot1=content1;
	string prot2=content2;
	/*string prot1 = argv[1];
	string prot2 = argv[2];*/

	int len1 = prot1.size();
	int len2 = prot2.size();

	//pos is a object Position* we will use to build our matrix
	Position* pos;

	//First step of the algorithm is to place zeros
	//in the first line and first column
	/*for (int i = 0; i <= len1 ; i++){
		pos=new Position(i,0);
		matrix->addPosition(pos);
	}
	for (int j = 0; j <= len2; j++){
		pos=new Position(0,j);
		matrix->addPosition(pos);
	}

	matrix->setupMax(len1+1, len2+1);
	*/
	//Define variables used in order to fill the ScoringMatrix
	int match;
	int up = 0;
	int left = 0;
	int diag = 0;

	/*int upgap;
	int leftgap;*/


	/*int maxX = 0;
	int maxY = 0;*/
	int maxValue = 0;
	vector<vector<int>> matrice;//[len1+1][len2+1];
	vector<int> tempv;
	tempv.assign(len2+1,0);
	matrice.assign(len1+1,tempv);
	int maxLine[len2+1];
	int maxColumn[len1+1];
	int posXmaxLine[len1+1];
	int posYmaxColumn[len2+1];
	for (int i = 0; i < len1+1; i++){
		posXmaxLine[i] = 0;
		maxColumn[i] = 0;
		matrice[i][0] = 0;
	}
	for (int j = 0; j < len2+1; j++){
		posYmaxColumn[j] = 0;
		maxLine[j] = 0;
		matrice[0][j] = 0;
	}
	int temp[4];
	clock_t begin = clock();
	int tempVal;
	int mline;
	int mcolumn;
	for (int i = 1; i < len1+1; i++){
		for (int j = 1; j < len2+1; j++){
			mcolumn = maxColumn[i];
			mline = maxLine[j];
			temp[0] = mcolumn - gap_op - gap_ex*(j - posYmaxColumn[i]);
			temp[1] = mline - gap_op - gap_ex*(i - posXmaxLine[j]);
			temp[2] = matrice[i-1][j-1] + blosum->get(prot1.at(i-1), prot2.at(j-1));
			temp[3] = 0;
			tempVal = findMax(temp,4);
			matrice[i][j] = tempVal;
			if (tempVal>mline){
				maxLine[j] = tempVal;
				posXmaxLine[j] = i;
			}
			if (tempVal > mcolumn){
				maxColumn[i] = tempVal;
				posYmaxColumn[i] = j;
			}
			/*if (tempVal > maxValue){
				maxValue = matrice[i][j];
				maxX = i;
				maxY = j;
			}*/
		}
	}
	maxValue = findMax2(maxColumn, maxLine, len2+1, len1+1);
	//int root[len1+1][len2+1];

	/*Position* maxPos = new Position();
	//Position* allPos = (Position*)malloc(sizeof(Position)*(len1+1)*(len2+1));

	/*cout << (len1+1)*(len2+1) << endl;
	for (int i = 0; i < (len1+1)*(len2+1); i++){
		allPos[i].setX(0);
		allPos[i].setY(0);
		allPos[i].setValue(0);
	}
	cout << "or here" << endl;*/
	//We will go in every case to complete the matrix
	/*for (int i = 1; i <= len1 ; i++){
		for (int j = 1; j <= len2 ; j++){
			//pos = new Position(i,j);
			//cout << i << " " << j << " "<< (i-1)*len1+(j-1) << endl;
			/*pos = &allPos[(i-1)*len1+(j-1)];
			pos->setX(i);
			pos->setY(j);*/
			//matrix->addPosition(pos);
			/*upgap = gap_op;
			leftgap = gap_op;
			//We check if the root is from the upper case.
			//If yes, we check if the root of the upper case is also its own upper case
			//If yes, this is an expansion and so, the penality due to this gap fall to 1
			if (j>=2 && matrix->getPosition(i,j-1)->getRoot().size() >= 1 && matrix->getPosition(i,j-1)->getRoot()[0] == matrix->getPosition(i,j-2)){
				upgap = gap_ex;
			}
			//Same for the left case
			if (i>=2 && matrix->getPosition(i-1,j)->getRoot().size() >= 1 && matrix->getPosition(i-1,j)->getRoot()[0] == matrix->getPosition(i-2,j)){
				leftgap = gap_ex;
			}
			up = matrix->getValue(i,j-1) - upgap; // we substract the gap to the value of the upper case
			left = matrix->getValue(i-1,j) - leftgap; //Same for the left case
			diag = matrix->getValue(i-1,j-1); // We just take the value of the left up diagonal case
			*/

			//We find the value in the BlosumMatrix
			/*up = matrix->getMaxColumn(i)->getValue() - (gap_op + gap_ex*matrix->getDistYWithMax(pos));
			left = matrix->getMaxLine(j)->getValue() - (gap_op + gap_ex*matrix->getDistXWithMax(pos));
			match=blosum->get(prot1.at(i-1), prot2.at(j-1));
			diag = matrix->getValue(i-1,j-1);
			diag+=match; // Add the result of the blosum matrix to the diagonal value
*/
			//Now, we find the maximum value. If they are all negative, we take 0.
			/*int highValue = max(up,max(left,max(diag,0)));
			//We set it to the actual position
			pos->setValue(highValue);

			matrix->setMaxLine(j, pos);
			matrix->setMaxColumn(i, pos);*/
			//We set the links between the roots and targets
			/*if (highValue == diag){
				matrix->setRootTarget(matrix->getPosition(i-1,j-1), pos);
			}
			else if (highValue == up){ //Attention, on considère qu'il ne peut y avoir qu'une seule root
				matrix->setRootTarget(matrix->getPosition(i,j-1), pos);
			}
			else if (highValue == left){
				matrix->setRootTarget(matrix->getPosition(i-1,j), pos);
			}*/
			//If this value is greater than the previous maximum, we keep the position and the value
			/*if (highValue > maxPos->getValue()){
				maxPos = pos;
			}*/
		/*}
	}*/
	clock_t end = clock();
	double time = double(end - begin)/CLOCKS_PER_SEC;
	int test = 0;
	clock_t begin2 = clock();
	for (int i = 0; i <= len1; i++){
		for (int j = 0; j <= len2; j++){
			test++;
		}
	}
	clock_t end2 = clock();
	double time2 = double(end2-begin2)/CLOCKS_PER_SEC;
	cout << "Le test en double for met " << time2 << " secondes" << endl;
	//int score = maxPos->getValue();
	int score = maxValue;
	double lambda = 0.267;
	double k = 0.041;
	double bitscore = double(score);
	bitscore = (lambda*bitscore - log(k))/log(2);
	cout << "Score donne : "<< score << " en " << time << " secondes"<< endl;
	//free(allPos);
	//Pas oublier de delete !!!!

	//matrix->print();
	//string res = findPath(maxPos, prot1, prot2, blosum);
	//cout << "The alignement of " << argv[1] << " and " << argv[2] << " gives " << res << " and gets a score of "<< maxPos->getValue() << endl;
	return 0;
}
