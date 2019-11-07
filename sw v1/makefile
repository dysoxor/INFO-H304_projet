all: smith;

smith: smith.o BlosumMatrix.o;\
g++ -o smith smith.o BlosumMatrix.o
smith.o : smith.cpp BlosumMatrix.o;\
g++ -c smith.cpp
BlosumMatrix.o : BlosumMatrix.cpp;\
g++ -c BlosumMatrix.cpp
test: test.cpp BlosumMatrix.o;\
	g++ -c test.cpp BlosumMatrix.o
clean: ;\
rm *.o
