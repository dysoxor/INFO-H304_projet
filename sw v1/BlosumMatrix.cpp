#include "BlosumMatrix.h"
using namespace std;


BlosumMatrix::BlosumMatrix(){

};

void BlosumMatrix::setup(string pathToBlosumMatrix){
	ifstream file(pathToBlosumMatrix);
	string line;
	bool first_line = true;
	int value;
	int n_ligne = 0;
	int n_colonne = 0;
	if (file.is_open()){
		while(getline(file,line)){
			if (line != "#"){
				n_colonne=0;
				if (first_line){
					first_line=false;
					for (int i = 0; i < line.length();i++){
						if (line.at(i) != ' '){
							charToInt.insert(pair<char,int>(line.at(i),n_colonne));
							n_colonne++;
						}
					}
				} else {
					line.erase(0,1);//on enleve la lettre
					for (int i = 0; i< line.length();i++){
						if (line.at(i) != ' ' && line.at(i) != '-'){
							value = (int)line.at(i) -48; //ASCII digits commence à 48
							if (line.at(i-1) == '-'){
								value = -value;
							}
							matrice[n_ligne][n_colonne] = value;
							n_colonne++;
						}
					}
					n_ligne++;
				}
			}
		}
		file.close();
	} else {
		cout << "Problème dans l'ouverture du fichier blosum" << endl;
	}
};

const int BlosumMatrix::get(string aa1, string aa2){
	int line = charToInt[aa1.at(0)];
	int colonne = charToInt[aa2.at(0)];
	return matrice[line][colonne];
};

const void BlosumMatrix::print(){
	for (int i = 0; i < matrice.size(); i++){
		for (int j = 0; j < matrice[0].size(); j++){
			cout << matrice[i][j] << " ";
		}
		cout << endl;
	}
};
