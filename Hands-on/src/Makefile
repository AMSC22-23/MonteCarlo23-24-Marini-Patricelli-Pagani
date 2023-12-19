.PHONY: all
all: main

main: main.cpp
	mpicxx main.cpp -I../muparser/include -L../muparser/build -lmuparser -Wall -fopenmp

.PHONY: clear
clear:
	rm -f main
