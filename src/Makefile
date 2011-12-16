FC	= gfortran
CC	= gcc
LD	= $(FC)
F77	= $(FC)
SWP	= swplist
RM	= /bin/rm -f
PROC	= ip27
OLEVEL	= -O0
FOPTS   = 
FFLAGS	= $(OLEVEL) $(FOPTS) -I/usr/include
CFLAGS	= $(COPTS)
LIBS	= 
LDFLAGS	= -lcgns $(FFLAGS)

FOBJS   = mesh_var.o read_cgns.o write_cgns.o mod_coord.o\
          sort3.o 

OBJS	= $(FOBJS)

DEFAULT:	random_blade setbc

random_blade:	$(OBJS) random_blade.o
	$(LD) -o $@ $(LDFLAGS) random_blade.o $(OBJS) $(LIBS)

setbc:	$(OBJS) setbc.o
	$(LD) -o $@ $(LDFLAGS) setbc.o $(OBJS) $(LIBS)

clean:
	$(RM) $(EXEC) $(OBJS) random_blade.o mesh_var.mod random_blade setbc

.SUFFIXES:
.SUFFIXES: .o .F .c .f .swp

.F.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.f.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.c.o:
	$(CC)  -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<

.F.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.f.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.c.swp:
	$(SWP) -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<

