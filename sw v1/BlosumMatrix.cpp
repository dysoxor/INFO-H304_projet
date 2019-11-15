#include "BlosumMatrix.h"
using namespace std;

//Create BlosumMatrix
BlosumMatrix::BlosumMatrix(string pathToBlosumMatrix){
	cout<<"Matrix created with " << pathToBlosumMatrix<< endl;
	BlosumMatrix::setup(pathToBlosumMatrix);
};

//Setup it
void BlosumMatrix::setup(string pathToBlosumMatrix){
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
							charToInt.insert(pair<char,int>(line.at(i),n_colonne));
							n_colonne++;
						}
					}
					matrix.assign(charToInt.size(), vector<int> (charToInt.size(),0));
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
							matrix[n_line][n_column] = value; //Set the value
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

const int BlosumMatrix::get(char aa1, char aa2){
	int line = charToInt[aa1];
	int column = charToInt[aa2];
	return matrix[line][column];
};

const int BlosumMatrix::get(string aa1, string aa2){
	return BlosumMatrix::get(aa1.at(0), aa2.at(0));
};

const void BlosumMatrix::print(){
	for (int i = 0; i < matrix.size(); i++){
		for (int j = 0; j < matrix[0].size(); j++){
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
};
