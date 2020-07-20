CC = g++
STD = -std=c++11
LLDB = -g
BOOST = -lboost_serialization
PYTHON = -lboost_python3-py36
FMATH = -fomit-frame-pointer -fno-operator-names -msse2 -mfpmath=sse -march=native
GFLAGS = -lglog -lgflags
#-I/root/boost -L/root/boost/stage/lib -I/usr/include/boost/ 
INCLUDE = -I/usr/include/
LDFLAGS = `python3-config --includes` `python3-config --ldflags`

cstm:
	$(CC) -O3 $(STD) -o cstm src/model.cpp $(BOOST) $(INCLUDE) $(FMATH) $(GFLAGS) 

install:
	$(CC) -O3 $(STD) -DPIC -shared -fPIC -o pycstm.so pycstm.cpp $(INCLUDE) $(LDFLAGS) $(PYTHON) $(BOOST) $(FMATH) $(GFLAGS)

test:
	$(CC) -O3 -Wall -o test src/model.cpp $(LLDB) $(BOOST) $(FMATH) $(GFLAGS)

clean:
	rm -f cstm

.PHONY: clean
