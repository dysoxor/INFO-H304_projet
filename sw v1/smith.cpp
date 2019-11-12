
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "BlosumMatrix.h"
#include "ScoringMatrix.h"

using namespace std;


string findPath(Position* pos, string prot1, string prot2, int gap, BlosumMatrix* blosum){

		cout << pos->getX() << " " << pos->getY() << endl;
		if (pos->getValue() == 0){
			return "";
		} else {
		string returnValue;
		vector<Position*> rootPosVector = pos->getRoot();
		Position* rootPos = rootPosVector[0];
		string res = findPath(rootPos, prot1, prot2, gap, blosum);
		bool aa1 = true;
		bool aa2 = true;
		if (rootPos->getX() == pos->getX()){
			aa1=false;
		} else if (rootPos->getY() == pos->getY()){
			aa2 = false;
		}
		if (res != ""){
			returnValue+=res;
			returnValue+= "/";
		}
		if (aa1){
			returnValue+=prot1[pos->getX()-1];
		} else {
			returnValue+="-";
		}
		if (aa2){
			returnValue+=prot2[pos->getY()-1];
		} else {
			returnValue+="-";
		}
		return returnValue;
	}
}




int main(int argc, char** argv){
	if (argc != 3){
		printf("Erreur, pas le bon nombre d'arguments rentrés\n");
		printf("Nombre d'arguments attendus : 2\n");
		printf("Nombre d'arguments donnés : %i\n", argc);
		return EXIT_FAILURE;

	}

	BlosumMatrix* blosum = new BlosumMatrix("blosum62");
	ScoringMatrix* matrice = new ScoringMatrix();

	string prot1=argv[1];
	string prot2=argv[2];

	int len1 = prot1.size();
	int len2 = prot2.size();
	Position* pos;
	for (int i = 0; i <= len1 ; i++){
		pos=new Position(i,0);
		matrice->addPosition(pos);
	}
	for (int j = 0; j <= len2; j++){
		pos=new Position(0,j);
		matrice->addPosition(pos);
	}


	int gap = 6;
	int match;
	//int mismatch;
	int up = 0;
	int left = 0;
	int diag = 0;

	int upgap;
	int leftgap;
	//int nb_of_max = 0; //pour un seul max pour l'instant
	Position* maxPos = new Position();
	for (int i = 1; i <= len1 ; i++){
		for (int j = 1; j <= len2 ; j++){
			pos = new Position(i,j);
			matrice->addPosition(pos);
			upgap = 11;
			leftgap = 11;
			if (j>=2 && matrice->getPosition(i,j-1)->getRoot().size() >= 1 && matrice->getPosition(i,j-1)->getRoot()[0] == matrice->getPosition(i,j-2)){
				upgap = 1;
				cout << "Expansion up en " << i << " " << j << endl;
			}
			if (i>=2 && matrice->getPosition(i-1,j)->getRoot().size() >= 1 && matrice->getPosition(i-1,j)->getRoot()[0] == matrice->getPosition(i-2,j)){
				leftgap = 1;
				cout << "Expansion left en "  << i << " " << j << endl;
			}
			/*up = tableau[i][j-1] + gap;
			left = tableau[i-1][j] + gap;
			diag = tableau[i-1][j-1];*/
			/*up = matrice->getValue(i,j-1) - gap;
			left = matrice->getValue(i-1,j) -gap;
			diag = matrice->getValue(i-1,j-1);*/
			up = matrice->getValue(i,j-1) - upgap;
			left = matrice->getValue(i-1,j) - leftgap;
			diag = matrice->getValue(i-1,j-1);
			/*match = 3;
			mismatch = -3;*/
			match=blosum->get(prot1.at(i-1), prot2.at(j-1));
			//cout << prot1.at(i-1) << prot2.at(j-1) << " " <<match << endl;
			/*if (prot1[i-1] == prot2[j-1]) //On a du ajouter une colonne et une ligne vide
				diag+=match;
			else
				diag+=mismatch;*/
			diag+=match;
			//tableau[i][j] = max(up,max(left,max(diag,0)));
			int highValue = max(up,max(left,max(diag,0)));
			pos->setValue(highValue);
			if (highValue == diag){
				matrice->setRootTarget(matrice->getPosition(i-1,j-1), pos);
			}
			else if (highValue == up){ //Attention, on considère qu'il ne peut y avoir qu'une seule root
				matrice->setRootTarget(matrice->getPosition(i,j-1), pos);
			}
			else if (highValue == left){
				matrice->setRootTarget(matrice->getPosition(i-1,j), pos);
			}

			if (highValue > maxPos->getValue()){
				maxPos = pos;
				cout<< maxPos->getValue()<< endl;

			}
		}
	}

	//Trouver les maxs (il faudra utiliser un meilleur algorithme pour le trouver)


	matrice->print();
	string res = findPath(maxPos, prot1, prot2, gap, blosum);
	cout << "L'alignement de " << argv[1] << " et " << argv[2] << " donne " << res << " et a un score de "<< maxPos->getValue() << endl;
	return 0;
}
