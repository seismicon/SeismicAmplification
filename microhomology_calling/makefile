CCC     = g++
CCFLAGS = -O2 -w -std=c++11

%.o: %.cc makefile installdir.h
	@echo [Compiling] $@ "\c"
	@$(CCC) $(CCFLAGS) -o $@ -c $<
	@echo " ...done."



MODULES = 	nibtools.o


all:	mh-caller 

mh-caller: mh-caller.o mh-caller.cc $(MODULES) makefile
	@echo [Linking] $@ "\c"
	@$(CCC) $(CCFLAGS)  -o mh-caller mh-caller.o $(MODULES)
	@echo " ...done."

installdir.h:
	@echo "Creating installdir.h..."
	@echo "#define INSTALLDIR \""`pwd`"\"" > installdir.h
	@echo "...done."

clean:
	@rm -f *.o installdir.h
