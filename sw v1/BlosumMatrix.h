#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

class BlosumMatrix{
	private:
		vector<vector<int>> matrice;
		map<char,int> charToInt;
	public:
		BlosumMatrix();
		void setup(string pathToBlosumMatrix);
		const int get(string aa1, string aa2);
		const void print();
};
