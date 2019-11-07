all: smith;

smith: smith.o BlosumMatrix.o;\
g++ -o smith smith.o BlosumMatrix.o
smith.o : smith.cpp BlosumMatrix.o;\
g++ -c smith.cpp
BlosumMatrix.o : BlosumMatrix.cpp;\
g++ -c BlosumMatrix.cpp
clean: ;\
rm *.o
