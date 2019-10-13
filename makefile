CC = gcc
LDFLAGS = -lglog -lgflags

cstm:
	$(CC) -o cstm src/cstm.cpp src/model.cpp $(INCLUDE) $(LDFLAGS)