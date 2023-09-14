include ../../mkfhosts/munari

OBJECTS= qpexample1.o pdcgm_SMatrix.o pdcgm_env.o
MYLIBS = -L ../../ipx-master/lib -lipx
OUTPUT = example-pdcgm
RPATH = -Wl,-rpath=../../ipx-master/lib
 
.c.o:
	$(CXX) $(CXXFLAGS) -c $(CFLAGS) $(DEFINE) $< -o $@

$(OUTPUT): $(OBJECTS)
	$(CXX)  $(OBJECTS) -o $(OUTPUT) $(MYLIBS) $(RPATH)

clean:
	\rm $(OBJECTS) $(OUTPUT)
