#########################################################################
# Makefile to create library for the eigensolver
#
# ================================
# 1 - set the compiler to be used
# 2 - create the program typing:
#          make           to create the dynamic version
#          make dynamic   to create the dynamic version
#          make static    to create the static version
#          make all       to create both versions
# 3 - other choices:
#          make clean     removes object files
#          make build     creates a file containing all sources
#########################################################################
#     Marco       Wed Aug 14 11:45:37 MET DST 1996
#########################################################################

#
# Set the desired compiler
#
################### uncomment the desired option ########################
#FORTRAN      = f77
#FORTRAN      = g77
FORTRAN      = gfortran
#FORTRAN      = ifort
#########################################################################


##################### don't change anything below #######################

#
# Library name
#
LIBRARY      =  libeigen

#
# name of makefile and lists of sources
#

MAKEFILE     = makefile

FILES1 = dcopy  \
         dlaneg \
         dlarrc \
         dlarrv \
         dlasq4 \
         dstegr \
         lsame  \
         disnan \
         dlanst \
         dlarrd \
         dlartg \
         dlasq5 \
         dstemr \
         xerbla \
         dlae2  \
         dlapy2 \
         dlarre \
         dlaruv \
         dlasq6 \
         dsteqr \
         dlaebz \
         dlar1v \
         dlarrf \
         dlascl \
         dlasr  \
         dswap  \
         dlaev2 \
         dlarnv \
         dlarrj \
         dlaset \
         dlasrt \
         ieeeck \
         dlaisnan \
         dlarra \
         dlarrk \
         dlasq2 \
         dlassq \
         ilaenv \
         dlamch \
         dlarrb \
         dlarrr \
         dlasq3 \
         dscal  \
         iparmq \
	 vandermonde \
         dgesv \
         dgetrf \
         dgetf2 \
	 dlaswp \
         dgetrs \
         dgemm \
         idamax \
         dger \
         dtrsm 

#
# working directories
#
SRCDIR       = src
OBJDIRSTA    = obj/sta
OBJDIRDYN    = obj/dyn
LIBDIR       = ./

#
# options for preprocessor, compiler and linker
#
FSTA         = -O5
FDYN         = -O5 -fPIC
FFLAGSSTA    = $(FSTA)
FFLAGSDYN    = $(FDYN)
#LDFLAGSFAST  = -dalign -libmil -native -xlic_lib=sunperf
#LDFLAGSDEBG  = -g -C -dalign -xlic_lib=sunperf

#
# basic commands for preprocessing, compiling and linking
#
FC           = $(FORTRAN) -c
ARCHI        = ar rcs
LINK         = $(FORTRAN) -shared

#
# define name of files
#
SOURCES      = $(FILES1:%=$(SRCDIR)/%.f )
OBJECTSSTA   = $(FILES1:%=$(OBJDIRSTA)/%.o )
OBJECTSDYN   = $(FILES1:%=$(OBJDIRDYN)/%.o )
STA          = $(LIBRARY:%=$(LIBDIR)/%.a)
DYN          = $(LIBRARY:%=$(LIBDIR)/%.so)


#
# commands
#
default : dynamic
all     : dynamic static
dynamic : $(DYN)
static  : $(STA)

#
# generate executable fast/debug version:
#
$(STA)  : $(OBJECTSSTA) $(MAKEFILE)
	$(ARCHI) $@ $(OBJECTSSTA)
$(DYN)  : $(OBJECTSDYN) $(MAKEFILE)
	$(LINK)  -o $@ $(OBJECTSDYN)

#
# generate object files for fast version:
#
$(OBJDIRSTA)/%.o : $(SRCDIR)/%.f
	$(FC) $(FFLAGSSTA) -o $@ $<

#
# generate object files for debug version:
#
$(OBJDIRDYN)/%.o : $(SRCDIR)/%.f
	$(FC) $(FFLAGSDYN) -o $@ $<

#
# rm o-files(debug+fast) if requested
#
clean :
	rm $(OBJDIRDYN)/*.o #$(OBJDIRSTA)/*.o 

