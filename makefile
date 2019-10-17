CC = clang++
LLDB = -g
BOOST = -lboost_serialization
PYTHON = -lboost_python37
FMATH = -fomit-frame-pointer -fno-operator-names -msse2 -mfpmath=sse -march=native
GFLAGS = -lglog -lgflags
INCLUDE = -I/usr/local/lib `/usr/local/Cellar/python/3.7.4_1/bin/python3.7-config --include`
LDFLAGS = `/usr/local/Cellar/python/3.7.4_1/bin/python3.7-config --ldflags`

cstm:
	$(CC) -o cstm src/model.cpp $(BOOST) $(FMATH) $(GFLAGS)

pycstm:
	$(CC) -Wall -DPIC -shared -fPIC -o pycstm.so pycstm.cpp $(PYTHON) $(INCLUDE) $(LDFLAGS) $(BOOST) $(FMATH) $(GFLAGS)

test:
	$(CC) -o cstm src/model.cpp $(LLDB) $(BOOST) $(FMATH) $(GFLAGS)

clean:
	rm -f cstm

.PHONY: clean