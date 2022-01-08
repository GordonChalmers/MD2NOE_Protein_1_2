GLYLIB_BASE = /home3/chalmers/Desktop/GLYLIB
GLYLIB=$(GLYLIB_BASE)/lib
GLYINC=$(GLYLIB_BASE)/inc
CFLAGS=-l glylib -Wall -lm -fpermissive -lpthread


EXEFILE=-o noe_curve
SOURCE=noe_curve.cc

default:
	g++ $(SOURCE) -I $(GLYINC) -L $(GLYLIB) $(CFLAGS) $(EXEFILE)

gdb:
	g++ $(SOURCE) -I$(GLYINC) -L$(GLYLIB) $(CFLAGS) $(EXEFILE) -g


