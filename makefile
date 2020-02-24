all: secular dev

PATH_TO_BOOST=./boost_1_70_0/
PATH_TO_SPACEHUB=./
CXX=g++
secular:
	${CXX} -std=c++17 -march=native  -O3 -o secular main.cpp -I${PATH_TO_BOOST} -pthread

dev:
	${CXX} -std=c++17 -march=native  -O3 -o dev main-dev.cpp -I${PATH_TO_BOOST} -pthread

	
clean:
	rm secular
