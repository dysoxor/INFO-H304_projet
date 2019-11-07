#include "Position.h"

Position::Position(){
  x=0;
  y=0;
  value = 0;
};

Position::Position(int posX, int posY){
  x=posX;
  y=posY;
  value = 0;
}

const vector<Position*> Position::getRoot(){
  return root;
};

const vector<Position*> Position::getTarget(){
  return target;
};

void Position::addRoot(Position* rootPos){
  root.push_back(rootPos);
};

void Position::addTarget(Position* targetPos) {
  target.push_back(targetPos);
};

const int Position::getX(){
  return x;
};

const int Position::getY(){
  return y;
};

void Position::setX(int newX){
  x = newX;
};

void Position::setY(int newY){
  y = newY;
};

const int Position::getValue(){
  return value;
};

void Position::setValue(int newValue){
  value = newValue;
};
