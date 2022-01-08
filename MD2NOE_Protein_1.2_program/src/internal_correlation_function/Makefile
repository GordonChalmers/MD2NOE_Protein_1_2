GLYLIB_BASE = /home3/chalmers/Desktop/GLYLIB
GLYLIB=$(GLYLIB_BASE)/lib
GLYINC=$(GLYLIB_BASE)/inc
CFLAGS=-l glylib -Wall -lm -fpermissive -lpthread


EXEFILE=-o internal_correlation_function
SOURCE=internal_correlation_function.cc

default:
	g++ $(SOURCE) -I $(GLYINC) -L $(GLYLIB) $(CFLAGS) $(EXEFILE)

gdb:
	g++ $(SOURCE) -I$(GLYINC) -L$(GLYLIB) $(CFLAGS) $(EXEFILE) -g


