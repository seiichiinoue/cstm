CC = clang++
LLDB = -g
BOOST = -lboost_serialization
FMATH = -fomit-frame-pointer -fno-operator-names -msse2 -mfpmath=sse -march=native
LDFLAGS = -lglog -lgflags

cstm:
	$(CC) -o cstm src/model.cpp $(LLDB) $(INCLUDE) $(BOOST) $(FMATH) $(LDFLAGS)

clean:
	rm -f cstm

.PHONY: clean