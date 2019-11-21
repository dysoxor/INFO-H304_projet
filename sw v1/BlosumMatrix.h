#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

class BlosumMatrix{
	private:
		vector<vector<int>> matrix;
		map<char,int> charToInt;//Map linking a amino acid to a column/line
	public:
		BlosumMatrix(string pathToBlosumMatrix);
		void setup(string pathToBlosumMatrix);
		const int get(string aa1, string aa2);
		const int get(char aa1, char aa2);
		const void print();
};
