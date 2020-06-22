defaults: main

#CFLAGS          = -Wall -O -Wuninitialized 
CFLAGS           = 
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =

BIN_DIR=.

#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# obj files
OBJS = \
${BIN_DIR}/main.o \
${BIN_DIR}/myformfunction.o \
${BIN_DIR}/initial.o \

${BIN_DIR}/%.o : %.c def.h
	-${PETSC_COMPILE_SINGLE} $< -o $@
main: ${OBJS} chkopts
	-${CLINKER} -o ${BIN_DIR}/main ${OBJS} ${PETSC_LIB}
clear:
	-${RM} ${ONJS} ${BIN_DIR}/main
run:
	-@${MPIEXEC} -n 6 ./main -da_grid_x 200 -da_grid_y 200 -ep 1
run32:
	-@${MPIEXEC} -n 32 ./main -da_grid_x 200 -da_grid_y 200 -ep 0.1 
run128:
	-@${MPIEXEC} -n 128 ./main -da_grid_x 400 -da_grid_y 400 -ep 0.1
run512:
	-@${MPIEXEC} -n 512 ./main -da_grid_x 800 -da_grid_y 800 -ep 0.1
run2048:
	-@${MPIEXEC} -n 2048 ./main -da_grid_x 1600 -da_grid_y 1600 -ep 0.1
