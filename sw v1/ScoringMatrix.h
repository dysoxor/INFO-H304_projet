#include <vector>
#include <map>
#include "Position.h"
using namespace std;

class ScoringMatrix{
private:
  map<pair<int, int>,Position*> tableau;
public:
  ScoringMatrix();
  Position* getPosition(int x, int y);
  void setValue(int value, int x, int y);
  const int getValue(int x, int y);
  void addPosition(Position* pos);
  void setRootTarget(Position* rootPos, Position* targetPos);

};
