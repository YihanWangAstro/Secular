all: secular init_format

PATH_TO_BOOST=./boost_1_70_0/
PATH_TO_SPACEHUB=./
CXX=g++
secular:
	${CXX} -std=c++17 -march=native  -O3 -o secular main.cpp -I${PATH_TO_BOOST} -pthread

init_format:
	${CXX} -std=c++17 -march=native  -O3 -o format initial_format.cpp

clean:
	rm secular format
