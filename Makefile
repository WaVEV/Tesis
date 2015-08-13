# Binary file
BIN			= pozo-1

# Flags
# CFLAGS			= -std=gnu99 -Wall -Wextra
CFLAGS			= -std=gnu99 -Wall -Wextra -fopenmp  -g -O3
LDFLAGS			= -llapack -lm 

# Default Values
Rmin			= 0.0
Rmax			= 50.0f
#l			= 2048
#l			= 510
kord			= 5
r1			= 5.0f
r2			= 10.0f
me			= 1.f
intg			= 500
nev			= 15
lmax			= 0
lambda_in		= 0.0f
lambda_fin		= 20.0f
numero_puntos_lambda	= 200
NUM_THREADS             = 8
OFILE			= [CPU,$(Q),$(Rmin),$(Rmax),$(l),$(kord),$(r1),$(r2),$(me),$(intg), \
			   $(nev),$(lmax),$(lambda_in),$(lambda_fin),$(numero_puntos_lambda)].dat

# Simulation Parameters
#PARAMETERS		= -DQ=$(Q) -DRmin=$(Rmin) -DRmax=$(Rmax) -Dl=$(l) -Dkord=$(kord) \
#			  -Dr1=$(r1) -Dr2=$(r2) -Dme=$(me) -Dintg=$(intg) -Dnev=$(nev) -Dlamx=$(lmax) \
#			  -Dlambda_in=$(lambda_in) -Dlambda_fin=$(lambda_fin) -Dnumero_puntos_lambda=$(numero_puntos_lambda) -DNUM_THREADS=$(NUM_THREADS)

# Compilers
CC			= gcc
LINKER			= gcc

# Files
C_SOURCES		= $(BIN).c
HEADERS			=
C_OBJS			= $(patsubst %.c, %.o, $(C_SOURCES))


# Rules
$(BIN): clean $(C_OBJS) $(HEADERS)
	$(LINKER) $(CFLAGS) -o $(BIN) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS) 
	$(CC) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES) $(PARAMETERS)

run: $(BIN)
	./$(BIN) # > $(OFILE) &

clean:
	rm -f $(BIN) *.o
