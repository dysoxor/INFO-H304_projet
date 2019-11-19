
#include <vector>

using namespace std;

class Position{
private:
  int x;
  int y;
  int value;
  vector<Position*> root; //Vector of pointers pointing to its roots
  vector<Position*> target;//Vector of pointers pointing to its targets
public:
  Position();
  Position(int posX, int posY);
  const vector<Position*> getRoot();
  const vector<Position*> getTarget();
  void addRoot(Position* rootPos);
  void addTarget(Position* targetPos);
  const int getX();
  const int getY();
  void setX(int newX);
  void setY(int newY);
  const int getValue();
  void setValue(int newValue);
};
