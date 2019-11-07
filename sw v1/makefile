all: smith;

smith: smith.o BlosumMatrix.o;\
g++ -o smith.o
smith.o : smith.cpp;\
g++ -c smith.cpp
BlosumMatrix.o : BlosumMatrix.cpp;\
g++ -c BlosumMatrix.cpp
clean:\
rm *.o
