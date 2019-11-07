
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "BlosumMatrix.h"

using namespace std;

typedef struct{
	int x;
	int y;
} Pos;

string findPath(Pos pos, vector<vector<int>> tableau, string prot1, string prot2, int gap, BlosumMatrix* matrice){
	Pos nextPos;
	//printf("Je suis a la pos %i %i\n",pos.x, pos.y);
	if (tableau[pos.x][pos.y] == 0 ){// || (pos.x == 0) || (pos.y == 0)){
		//printf("J'arrive à un 0\n");
		return "";
	} else {
		//printf("Je suis pas encore à un 0\n");
		int diag = tableau[pos.x-1][pos.y-1];
		int up = tableau[pos.x][pos.y-1] + gap;
		int left = tableau[pos.x-1][pos.y] + gap;
		int match=matrice->get(prot1.at(pos.x-1), prot2.at(pos.y-1));
		/*if (prot1[pos.x-1] == prot2[pos.y-1])
			diag += match;
		else
			diag+=mismatch;*/
		diag+=match;
		nextPos.x = pos.x - 1;
		nextPos.y = pos.y - 1;
		bool aa1 = true;
		bool aa2 = true;
		if (tableau[pos.x][pos.y] == left){
			//printf("Je viens de la gauche\n");
			nextPos.y ++;
			aa2=false;
		}
		else if (tableau[pos.x][pos.y] == up){
			//printf("Je viens d'en haut\n");
			nextPos.x ++;
			aa1=false;
		}
		else if(tableau[pos.x][pos.y] == diag){
			//printf("Je viens de la diagonale\n");
		}
		else {
			cout << "Problème dans le findpath" << endl;
		}
		string res = findPath(nextPos, tableau, prot1, prot2, gap, matrice);
		string returnValue;
		if (res != ""){
			returnValue+=res;
			returnValue+= "/";
		}
		if (aa1){
			returnValue+=prot1[pos.x-1];
		} else {
			returnValue+="-";
		}
		if (aa2){
			returnValue+=prot2[pos.y-1];
		} else {
			returnValue+="-";
		}
		//std::cout << returnValue << std::endl;
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
	BlosumMatrix* matrice = new BlosumMatrix();
	string prot1=argv[1];
	string prot2=argv[2];
	int len1 = prot1.size();
	int len2 = prot2.size();
	vector<vector<int>> tableau(len1+1, vector<int>(len2+1,0));
	/*for (int i = 0; i <= len1 ; i++)
		tableau[i][0]= 0;
	for (int j = 0; j <= len2; j++)
		tableau[0][j]=0;*/

	int gap = 2;
	int match;
	int up = 0;
	int left = 0;
	int diag = 0;

	int vmax = 0;
	//int nb_of_max = 0; //pour un seul max pour l'instant
	Pos pos;
	pos.x=0;
	pos.y=0;
	for (int i = 1; i <= len1 ; i++){
		for (int j = 1; j <= len2 ; j++){
			up = tableau[i][j-1] + gap;
			left = tableau[i-1][j] + gap;
			diag = tableau[i-1][j-1];
			match=matrice->get(prot1.at(i-1), prot2.at(j-1));
			/*if (prot1[i-1] == prot2[j-1]) //On a du ajouter une colonne et une ligne vide
				diag+=match;
			else
				diag+=mismatch;*/
			diag+=match;
			tableau[i][j] = max(up,max(left,max(diag,0)));
			if (tableau[i][j] > vmax){
				vmax = tableau[i][j];
				//nb_of_max = 1;
				pos.x = i;
				pos.y = j;
			} /*else if (tableau[i][j] == max){
				nb_of_max++;
			}*/
		}
	}

	//Trouver les maxs (il faudra utiliser un meilleur algorithme pour le trouver)


	/*int max = 0;
	//int nb_of_max = 0; //pour un seul max pour l'instant
	Pos pos;
	pos.x=0;
	pos.y=0;
	for (int i = 0; i <= len1; i++){
		for (int j = 0; j <= len2; j++){
			if (tableau[i][j] > max){
				max = tableau[i][j];
				pos.x = i;
				pos.y = j;
			} *//*else if (tableau[i][j] == max){
				nb_of_max++;
			}*/
		/*}
	}*/
	//printf("Max en %i %i avec comme valeur %i\n",pos.x, pos.y, max);

	string res = findPath(pos, tableau, prot1, prot2, gap, matrice);
	std::cout << "L'alignement de " << argv[1] << " et " << argv[2] << " donne " << res << " et a un score de "<< vmax << std::endl;

}
