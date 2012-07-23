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
          sort3.o read_cgns_stator.o write_cgns_stator.o mod_coord_stator.o\
          read_cgns_2D.o write_cgns_2D.o mod_coord_2D.o

OBJS	= $(FOBJS)

DEFAULT:	random_blade random_blade_stator random_blade_2D setbc

random_blade:	$(OBJS) random_blade.o
	$(LD) -o $@ $(LDFLAGS) random_blade.o $(OBJS) $(LIBS)

random_blade_stator:	$(OBJS) random_blade_stator.o
	$(LD) -o $@ $(LDFLAGS) random_blade_stator.o $(OBJS) $(LIBS)

random_blade_2D:	$(OBJS) random_blade_2D.o
	$(LD) -o $@ $(LDFLAGS) random_blade_2D.o $(OBJS) $(LIBS)

setbc:	$(OBJS) setbc.o
	$(LD) -o $@ $(LDFLAGS) setbc.o $(OBJS) $(LIBS)

clean:
	$(RM) $(EXEC) $(OBJS) random_blade.o mesh_var.mod random_blade setbc.o setbc\
                          random_blade_stator.o random_blade_stator random_blade_2D.o\
                          random_blade_2D

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

