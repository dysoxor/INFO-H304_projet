#include <vector>
#include <map>
#include "Position.h"
#include <iostream>
using namespace std;

class ScoringMatrix{
private:
  map<pair<int, int>,Position*> tableau;
  int sizeX;
  int sizeY;
public:
  ScoringMatrix();
  Position* getPosition(int x, int y);
  void setValue(int value, int x, int y);
  const int getValue(int x, int y);
  void addPosition(Position* pos);
  void setRootTarget(Position* rootPos, Position* targetPos);
  const void print();

};
