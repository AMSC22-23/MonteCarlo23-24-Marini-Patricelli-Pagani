.PHONY: all
all: main

# this is not even a good makefile command because does not exploit automatic rules
# but at least it works
main: main.cpp
	mpicxx main.cpp -I../muparser/include -L../muparser/build -lmuparser -Wall -fopenmp

.PHONY: clear
clear:
	rm -f main
