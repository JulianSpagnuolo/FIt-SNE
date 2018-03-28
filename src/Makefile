CXX = g++

.PHONY: all clean 

all: fast_tsne
	mkdir -p ../inst
	cp fast_tsne ../inst/

fast_tsne: sptree.cpp tsne.cpp nbodyfft.cpp
	$(CXX) -std=c++11 -O3 -pthread -lfftw3 -lm  -o fast_tsne sptree.cpp tsne.cpp nbodyfft.cpp

clean:
	rm fast_tsne
	rm ../inst/fast_tsne
