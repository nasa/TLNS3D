SOURCE_MPI=	\
	main.f	        \
	tlns3d.f	\
	bcflow.f	\
	bcturb.f	\
	bcvol.f		\
	disdr.f		\
	errexit.f	\
	getconv.f	\
	savsol2.f	\
        savsolavg.f     \
	savdata.f	\
	savtim.f	\
	sendpch.f	\
	force.f         \
	turbdr.f        \
	solve.f	

OBJECT_MPI = tlns3d.o bcflow.o bcturb.o bcvol.o	disdr.o	errexit.o \
getconv.o savsol2.o savdata.o savtim.o sendpch.o force.o turbdr.o \
solve.o	main.o savsolavg.o

SOURCE_MB =	\
	main.f	\
	tlns3d.f	\
	bcflow.f	\
	bcturb.f	\
	bcvol.f	 	\
	disdr.f		\
	errexit.f	\
	force.f 	\
	turbdr.f	\
	savsol2.f	\
        savsolavg.f     \
	savdata.f	\
	savtim.f	\
	solve.f	

OBJECT_MB = main.o	\
	tlns3d.o	\
	bcflow.o	\
	bcturb.o	\
	bcvol.o		\
	disdr.o		\
	errexit.o	\
	force.o 	\
	turbdr.o	\
	savsol2.o	\
        savsolavg.o     \
	savdata.o	\
	savtim.o	\
	solve.o	

DOUBLE = -DTLN_DBL -r8 
#DOUBLE =
FLAGS_MB = $(DOUBLE)
FLAGS_MPI = -DBUILD_MPI $(DOUBLE)

LIBDIR = lib_obj
MPIDIR = mpi_obj
MBDIR = mb_obj

MV = -mv -f
RM = -rm -f
CP = -cp -f
MKDIR = -mkdir

all: lib mb mpi

lib:
	$(MV) $(LIBDIR)/*.o .
	make_lib_cpp
	make -f makelib
	$(MKDIR) $(LIBDIR)
	$(MV) *.o $(LIBDIR)
	$(MV) *.f $(LIBDIR)

mpi:
	$(MV) $(MPIDIR)/*.o .
	make_mpi_cpp
	make -f makempi "OBJ=$(OBJECT_MPI)"  "F90=mpif77"
	$(MKDIR) $(MPIDIR)
	$(MV) *.o $(MPIDIR)
	$(MV) *.f $(MPIDIR)

mb:
	$(MV) $(MBDIR)/*.o .
	make_mb_cpp
	make -f makempi "OBJ=$(OBJECT_MB)" \
	 "PROG=tlns3dmb"
	$(MKDIR) $(MBDIR)
	$(MV) *.o $(MBDIR)
	$(MV) *.f $(MBDIR)

cleanmpi:
	$(RM) $(MPIDIR)/*.o

cleanmb:
	$(RM) $(MBDIR)/*.o

clean:
	$(RM) *.o *.mod $(MPIDIR)/*.o $(MBDIR)/*.o $(LIBDIR)/*.o

