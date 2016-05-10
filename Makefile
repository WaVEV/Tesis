# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/umfpack/sym directory makefile.

#include Makefile.inc
# ARPACK++ v1.2 2/20/2000
# c++ interface to ARPACK code.
# examples/umfpack/sym directory makefile.

# including other makefiles.

include Makefile.inc

# defining CSCMAT directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: pozo-2

# compiling and linking each symmetric problem.

pozo-2:	pozo-2.o
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -o pozo-2 pozo-2.o $(UMFPACK_LIB) $(ALL_LIBS) -fopenmp -O3

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o core usymreg usymshf usymgreg usymgshf usymgbkl usymgcay

# defining pattern rules.

%.o:	%.cc
	$(CPP) $(CPP_FLAGS) -I$(CSCMAT_DIR) -I$(UMFPACK_DIR) -c $< -fopenmp -O3
