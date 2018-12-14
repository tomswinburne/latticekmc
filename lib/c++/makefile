

all: KMCBackend.cpp
	g++ -std=c++11 -Wall -fPIC -O2 -c KMCBackend.cpp -o KMCBackend.o
	g++ -std=c++11 -shared -o libkmcsim.so KMCBackend.o
