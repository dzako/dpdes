# a list of all the programs in your project 
#PROGS = 2DIntegratorTest IndicesTests ExactFFTTest FFTTests PolyBd2DTests KSPOTest BurgersFPTest SHTest CHTest
PROGS = DBCPModelHetConProof
#PROGS = DBCP2DNonrigorous CHNonrigorous
#PROGS = IndicesTests 2DIntegratorTest

# a list of all your units to be linked with your programs (space separated)
OTHERS =

# directory where capd scripts are (e.g. capd-config)
CAPDBINDIR =/usr/local/bin/

# setting compiler and linker flags
CAPDFLAGS = `${CAPDBINDIR}capd-config --cflags`
CAPDLIBS = `${CAPDBINDIR}capd-config --libs`
CXXFLAGS += ${CAPDFLAGS} -O2 -Wall -frounding-math -ggdb

#CXXFLAGS += ${CAPDFLAGS} -O2 -Wall -D__USE_FILIB__ -frounding-math -ffloat-store
#THE FFLOAT-STORE ABOVE MAKES PROGRAM 3 x SLOWER WHEN FILIB INTERVALS WITH ROUNDING ARE USED

#CXXFLAGS += ${CAPDFLAGS} -O0 -Wall -D__USE_FILIB__ -frounding-math
#THE -O0 FLAG MAKES PROGRAM 5 X SLOWER COMPARED TO -O2

# directory where object and dependancy files will be created
OBJDIR = .obj/

#============ the following should not be changed =========

OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
OBJ_FILES = ${OTHERS_OBJ} ${PROGS:%=${OBJDIR}%.o}

.PHONY: all
all: ${PROGS}

# rule to link executables
${PROGS}: % : ${OBJDIR}%.o ${OTHERS_OBJ}
	${CXX} -o $@ $< ${OTHERS_OBJ} ${CAPDLIBS}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}

#rule to compile .cpp files and generate corresponding files with dependencies
${OBJ_FILES}: ${OBJDIR}%.o : %.cpp
	@mkdir -p ${OBJDIR}
	$(CXX) ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<

# rule to clean all object files, dependencies and executables
.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o ${OBJDIR}*.o.d ${PROGS}


