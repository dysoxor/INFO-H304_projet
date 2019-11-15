
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "BlosumMatrix.h"
#include "ScoringMatrix.h"

using namespace std;


string findPath(Position* pos, string prot1, string prot2, int gap, BlosumMatrix* blosum){
		if (pos->getValue() == 0){ //If we fall to 0, it is the end of the alignement
			return "";
		}
		else {
			string returnValue;
			vector<Position*> rootPosVector = pos->getRoot();
			Position* rootPos = rootPosVector[0]; // we consider only the first root
			string res = findPath(rootPos, prot1, prot2, gap, blosum);
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




int main(int argc, char** argv){
	if (argc != 3){
		cout << "Error, not the right number of parameters" << endl;
		cout << "Number of parameters expected : 2" << endl;
		cout <<"Number of parameters given : " << (argc - 1) << endl;
		return EXIT_FAILURE;
	}
	//creation of the blosum matrix
	BlosumMatrix* blosum = new BlosumMatrix("blosum62");

	//Initialization of the ScoringMatrix
	ScoringMatrix* matrix = new ScoringMatrix();

	//define the 2 constants : gap opening and gap expansion
	const int gap_op = 11;
	const int gap_ex = 1;

	//Define the 2 sequence and place them into 2 variables prot1 and prot2

	string prot1=argv[1];
	string prot2=argv[2];

	int len1 = prot1.size();
	int len2 = prot2.size();

	//pos is a object Position* we will use to build our matrix
	Position* pos;

	//First step of the algorithm is to place zeros
	//in the first line and first column
	for (int i = 0; i <= len1 ; i++){
		pos=new Position(i,0);
		matrix->addPosition(pos);
	}
	for (int j = 0; j <= len2; j++){
		pos=new Position(0,j);
		matrix->addPosition(pos);
	}


	//Define variables used in order to fill the ScoringMatrix
	int match;
	int up = 0;
	int left = 0;
	int diag = 0;

	int upgap;
	int leftgap;
	Position* maxPos = new Position();

	//We will go in every case to complete the matrix
	for (int i = 1; i <= len1 ; i++){
		for (int j = 1; j <= len2 ; j++){
			pos = new Position(i,j);
			matrix->addPosition(pos);
			upgap = gap_op;
			leftgap = gap_op;
			//We check if the root is from the upper case.
			//If yes, we check if the root of the upper case is also its own upper case
			//If yes, this is an expansion and so, the penality due to this gap fall to 1
			if (j>=2 && matrix->getPosition(i,j-1)->getRoot().size() >= 1 && matrix->getPosition(i,j-1)->getRoot()[0] == matrix->getPosition(i,j-2)){
				upgap = gap_ex;
				cout << "Expansion up en " << i << " " << j << endl;
			}
			//Same for the left case
			if (i>=2 && matrix->getPosition(i-1,j)->getRoot().size() >= 1 && matrix->getPosition(i-1,j)->getRoot()[0] == matrix->getPosition(i-2,j)){
				leftgap = gap_ex;
				cout << "Expansion left en "  << i << " " << j << endl;
			}
			up = matrix->getValue(i,j-1) - upgap; // we substract the gap to the value of the upper case
			left = matrix->getValue(i-1,j) - leftgap; //Same for the left case
			diag = matrix->getValue(i-1,j-1); // We just take the value of the left up diagonal case
			//We find the value in the BlosumMatrix
			match=blosum->get(prot1.at(i-1), prot2.at(j-1));
			diag+=match; // Add the result of the blosum matrix to the diagonal value

			//Now, we find the maximum value. If they are all negative, we take 0.
			int highValue = max(up,max(left,max(diag,0)));
			//We set it to the actual position
			pos->setValue(highValue);
			//We set the links between the roots and targets
			if (highValue == diag){
				matrix->setRootTarget(matrix->getPosition(i-1,j-1), pos);
			}
			else if (highValue == up){ //Attention, on considère qu'il ne peut y avoir qu'une seule root
				matrix->setRootTarget(matrix->getPosition(i,j-1), pos);
			}
			else if (highValue == left){
				matrix->setRootTarget(matrix->getPosition(i-1,j), pos);
			}
			//If this value is greater than the previous maximum, we keep the position and the value
			if (highValue > maxPos->getValue()){
				maxPos = pos;
			}
		}
	}
	matrix->print();
	string res = findPath(maxPos, prot1, prot2, gap, blosum);
	cout << "The alignement of " << argv[1] << " and " << argv[2] << " gives " << res << " and gets a score of "<< maxPos->getValue() << endl;
	return 0;
}
