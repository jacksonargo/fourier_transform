SOURCE = ./src/*.c
INC = ./inc
CFLAGS = -Wall -I$(INC) -lm
default:
	gcc -g $(SOURCE) -o transform $(CFLAGS)
