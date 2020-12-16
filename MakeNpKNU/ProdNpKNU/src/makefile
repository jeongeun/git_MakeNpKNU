ROOTCFLAGS = $(shell root-config --cflags)

CXXFLAGS = -O2 -Wall -fPIC

all: libNpKNU.so

libNpKNU.so: NpKNUDict.o NpKNU.o
	g++ -shared $(CXXFLAGS) $^ -o$@
	@rm -rf *Dict.h *Dict.cxx *Dict.h.gch *.o

NpKNU.o: NpKNU.cc
	g++ $(CXXFLAGS) $(ROOTCFLAGS) $^ -c -o$@

NpKNUDict.o: NpKNU.hh NpKNULinkDef.h
	@echo "Generating dictionary"
	#rootcint -f $@ -c $(CXXFLAGS) -p $^
	rootcling -f NpKNUDict.cxx -c $^
	g++ $(CXXFLAGS) $(ROOTCFLAGS) NpKNUDict.cxx -c -o$@

clean:
	rm -rf *.hh.gch  *Dict.h *Dict.cxx  *Dict.h.gch *.o *Dict_rdict.pcm lib*.so

