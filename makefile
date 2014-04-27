SOURCE = ./src/*.c
INC = ./inc
CFLAGS = -Wall -I$(INC) -lm -lgsl -lgslcblas
default:
	gcc -g $(SOURCE) -o transform $(CFLAGS)
