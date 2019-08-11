all: secular sa

PATH_TO_BOOST=./boost_1_70_0/
PATH_TO_SPACEHUB=./
CXX=g++
secular:
	${CXX} -std=c++17 -mavx -O3 -o secular main.cpp -I ${PATH_TO_BOOST} -pthread

sa:
	${CXX} -std=c++17 -mavx -O3 -o sa revision.cpp -I ${PATH_TO_BOOST} -pthread
	
clean:
	rm secular
