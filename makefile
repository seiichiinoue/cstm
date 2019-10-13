CC = clang++
BOOST = -lboost_serialization
FMATH = -fomit-frame-pointer -fno-operator-names -msse2 -mfpmath=sse -march=native
LDFLAGS = -lglog -lgflags

cstm:
	$(CC) -o cstm src/model.cpp $(INCLUDE) $(BOOST) $(FMATH) $(LDFLAGS)