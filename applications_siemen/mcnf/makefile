include ../../mkfhosts/munari

ARCH := $(shell getconf LONG_BIT)
OBJECTS= mcnf-pdcgm.o pdcgm_env.o pdcgm_SMatrix.o 
MYLIBS = -L ../../ipx-master/lib -lipx
OUTPUT = mcnf-pdcgm
RPATH = -Wl,-rpath=../../ipx-master/lib
 
%.o: %.ccp
	$(CXX) $(CXXFLAGS) -c $(CFLAGS) $(DEFINE) $< -o $@

$(OUTPUT): $(OBJECTS)
	$(CXX)   $(CFLAGS) $(OBJECTS) -o $(OUTPUT) $(MYLIBS) $(RPATH) -g

clean:
	\rm $(OBJECTS) $(OUTPUT)
