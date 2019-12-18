
#include <ctime>
#include <iostream>

using namespace std;
int main (int argc, char** argv){
  time_t actualTime = time(nullptr);
  cout << "Début : " << asctime(localtime(&actualTime)) << endl;
  clock_t begin = clock();


  for (int i = 0; i < 550000; i++){
    for (int j = 0; j < 5000; j++){
      for (int k = 0; k < 400; k++){
      }
    }
  }

  
  time_t endTime = time(nullptr);
  cout << "Fin : " << asctime(localtime(&endTime)) << endl;
  clock_t end = clock();
  cout << "#clocks par sec : " << (double)CLOCKS_PER_SEC << endl;
  cout << "Elapsed time : " << (double)(end-begin)/CLOCKS_PER_SEC << endl;
  cout << "# clocks : " << end-begin << endl;
}
