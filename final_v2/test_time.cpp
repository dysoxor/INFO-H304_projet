
#include <ctime>
#include <iostream>

using namespace std;
int main (int argc, char** argv){
  clock_t begin = clock();
  //int alpha = 0;
  for (int i = 0; i < 550000; i++){
    for (int j = 0; j < 5000; j++){
      for (int k = 0; k < 1000; k++){
        //alpha++;
      }
      //while (true);
    }
    //clock_t inter = clock();
    //cout << "estimated remaining time : " << (double)(inter-beg)*(550000-i)/(i)/CLOCKS_PER_SEC << endl;
    if (i == 300){
      clock_t inter = clock();
      cout << (double)(inter-begin)*(550000-300)/(300*CLOCKS_PER_SEC) << endl;
      while(true);
    }
  }
  clock_t end = clock();
  cout << "Elapsed time : " << (double)(end-begin)/CLOCKS_PER_SEC << endl;
  cout << "# clocks : " << end-begin << endl;
}
