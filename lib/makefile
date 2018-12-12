

all: engine.cpp
	g++ -std=c++11 -Wall -fPIC -O2 -c engine.cpp -o engine.o
	g++ -std=c++11 -shared -o libsim.so engine.o
