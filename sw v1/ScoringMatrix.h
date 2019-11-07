#include <vector>
#include <map>
using namespace std;

class ScoringMatrix{
private:
  map<pair<int, int>,Position*> tableau;
public:
  ScoringMatrix();
  const Position* getPosition(int x, int y);
  void setValue(int value, int x, int y);
  const int getValue(int x, int y);
  void addPosition(int x, int y, Position* pos);
  void setRootTarget(Position* rootPos, Position* targetPos);

};
