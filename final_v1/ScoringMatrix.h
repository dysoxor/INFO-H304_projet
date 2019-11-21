#include <vector>
#include <map>
#include <iostream>
#include "Position.h"

using namespace std;

class ScoringMatrix{
private:
  map<pair<int, int>,Position*> tableau;
  vector<Position*> maxLines;
  vector<Position*> maxColumns;
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
  void setMaxLine(int line, Position* pos);
  Position* getMaxLine(int line);
  void setMaxColumn(int column, Position* pos);
  Position* getMaxColumn(int coulumn);
  void setupMax(int len1, int len2);
  int getDistXWithMax(Position* pos);
  int getDistYWithMax(Position* pos);

};
