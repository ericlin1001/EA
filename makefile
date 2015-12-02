all:*
	g++ -g *.cpp
	./a.out

param:
	mpic++ -DALGORITHM=4 *.cpp
	mpirun -n 5 ./a.out
