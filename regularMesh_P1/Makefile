# by R.Lie, Nov 01, 2002

include /usr/local/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS += -L/usr/local/AFEPack/library/lib -lAFEPack

#all : 

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ possion_equation.o $(LDFLAGS) $(LIBS)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx
	-rm -f step1
	-rm -f step2
	-rm -f *.m

.PHONY : default clean
