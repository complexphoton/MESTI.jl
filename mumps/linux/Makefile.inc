#
#  This file is part of MUMPS 5.6.2, released
#  on Wed Oct 11 09:36:25 UTC 2023
#
################################################################################
#
#   Makefile.inc.generic
#
#   This defines some parameters dependent on your platform; you should
#   look for the approriate file in the directory ./Make.inc/ and copy it
#   into a file called Makefile.inc. For example, from the MUMPS root
#   directory, use
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc"
#   (see the main README file for details)
#
#   If you do not find any suitable Makefile in Makefile.inc, use this file:
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc" and modify it according
#   to the comments given below. If you manage to build MUMPS on a new platform,
#   and think that this could be useful to others, you may want to send us
#   the corresponding Makefile.inc file.
#
################################################################################


########################################################################
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. It is recommended to
# install other orderings. For that, you need to obtain the corresponding package
# and modify the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 6.0 and later) orderings are recommended.
#

#SCOTCHDIR  = ${HOME}/scotch_6.0
#ISCOTCH    = -I$(SCOTCHDIR)/include
#
# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

# The following is the path of the METIS on Hsu's group in USC Discovery cluster.
# Please modify the path based on your machine.
# It is strongly recommended to install METIS for 3D calculations. The 3D calculations can be
# extremely slow and time consuming without METIS. You should specify the library and include 
# directory for METIS in the following.
LMETISDIR = /project/cwhsu_38/shared/software/metis-5.1.0/build/Linux-x86_64/libmetis
IMETIS    = -I/project/cwhsu_38/shared/software/metis-5.1.0/include

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

LMETIS    = -L$(LMETISDIR) -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
# If you want to use Metis 4.X or an older version, you should use -Dmetis4 instead of -Dmetis
# or in addition with -Dparmetis (if you are using parmetis 3.X or older).
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis

ORDERINGSF  = -Dpord -Dmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################

########################################################################

# DEFINE HERE SOME COMMON COMMANDS, THE COMPILER NAMES, ETC...

# PLAT : use it to add a default suffix to the generated libraries
PLAT    =
# Library extension, + C and Fortran "-o" option
# may be different under Windows
LIBEXT  = .a
LIBEXT_SHARED  = .so
FPIC_OPT = -fPIC
SONAME = -soname
# Adapt/uncomment RPATH_OPT to avoid modifying
# LD_LIBRARY_PATH in case of shared libraries
# The following is the path of the MUMPS shared libraries on Hsu's group in USC Discovery cluster:
RPATH_OPT = -Wl,-rpath,/project/cwhsu_38/shared/software/mumps-5.6.2-par/lib/
OUTC    = -o 
OUTF    = -o 
# RM : remove files
RM = /bin/rm -f
# CC : C compiler
CC = mpiicc
# FC : Fortran 90 compiler
FC = mpiifort
# FL : Fortran linker
FL = mpiifort
# AR : Archive object in a library
#      keep a space at the end if options have to be separated from lib name
AR = ar vr 
#RANLIB = ranlib
RANLIB  = echo

# DEFINE HERE YOUR LAPACK and ScaLAPACK LIBRARIES

# Make this variable point to the path where the Intel MKL library is
# installed. It is set to the default install directory for Intel MKL.

# Default install directory for Intel MKL on a local machine is following:
# MKLROOT=/opt/intel/oneapi/mkl
# However, in different Linux clusters, the Intel MKL may be put in different
# path. If you are using in Lmod module system and Intel MKL are installed,
# you can use "module load intel", "module load intel-mkl", and "echo $MKLROOT".
# MKLROOT path should be printed out, if Intel MKL is available.
# The MKLROOT should be under the path ".../mkl", where "..." is the Intel oneapi path.
#
# If Lmod module system is not used in your cluster, then you can also try to
# find the path like ".../mkl", where "..." is the Intel oneapi directory path.
# You can specify the MKLROOT by "source .../mkl/bin/mklvars.sh intel64"
# and type "echo $MKLROOT". MKLROOT path should be printed out.
#
# If you still cannot figure out the path of the libraries, you may need to
# contact the maintenance team of your cluster.

# LAPACK and ScaLAPACK libraries are under the folder of $(MKLROOT)/lib/intel64

LAPACK = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
SCALAP = -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

# The next line should not be modified. They concern
# the parallel library provided by MUMPS.
LIBPAR = $(SCALAP) $(LAPACK)

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

# DEFINE HERE YOUR BLAS LIBRARY
# BLAS library is under the folder of $(MKLROOT)/lib/intel64
LIBBLAS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O -nofor-main -DBLR_MT -qopenmp -DGEMMT_AVAILABLE
OPTL    = -O -nofor-main -qopenmp
OPTC    = -O -qopenmp
#End Optimized options

INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =
