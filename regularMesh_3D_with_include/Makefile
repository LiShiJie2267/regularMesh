# by R.Lie, Nov 01, 2002

include /usr/local/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS += -L/usr/local/AFEPack/library/lib -lAFEPack

#all : 

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ possion_equation.o ../obj/RectangleDomain.o $(LDFLAGS) $(LIBS)

debug :main
	$(CXX) -o $@ possion_equation.o ../obj/RectangleDomain.o $(LDFLAGS) $(LIBS) -g
clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx
	-rm -f *.m
	-rm -f *.o
	-rm -f debug

.PHONY : default clean
