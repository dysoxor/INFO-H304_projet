#include "BlosumMatrix.cpp"
//#include "BlosumMatrix.cpp"

int main(int argc, char** argv){
	BlosumMatrix* matrice = new BlosumMatrix();
	//string path = argv[1];
	matrice->setup("blosum62");
	int res1 = matrice->get("A","A");
	int res2 = matrice->get('A','A');
	cout << res1 << " " << res2 << endl;
	//matrice->print();
	//delete matrice;

}
