# by R.Lie, Nov 01, 2002

include /usr/local/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS += -L/usr/local/AFEPack/library/lib -lAFEPack

#all : 
%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)
	mv $@ ../obj/
main : $(object)
	$(CXX) -o $@ ../obj/*.o $(LDFLAGS) $(LIBS) 
debug : main
	$(CXX) -o $@ ../obj/*.o $(LDFLAGS) $(LIBS) -g -I ../obj/*.o

clean : 
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx
	-rm -f *.m
	-rm -f debug

.PHONY : default clean

