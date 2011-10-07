CC=cc
#OPTS=-lm -O3 -march=pentium4 -mcpu=pentium4 -Wall -m32
#OPTS=-lm -Wall -m32
OPTS=-lm -Wall -O3 -m32

# top-level rule to compile the whole program
all: netclustwacnmdadiff

# program is made of several source files
netclustwacnmdadiff: netclustwacnmdadiff.o clust_cn.o
	${CC} ${OPTS} netclustwacnmdadiff.o clust_cn.o -o netclustwacnmdadiff
	/bin/rm -f netclustwacnmdadiff.o clust_cn.o
# rule for netclustwacnmdadiff.o
netclustwacnmdadiff.o: netclustwacnmdadiff.c lnet.h 
	${CC} ${OPTS} -c netclustwacnmdadiff.c

# rule for clust_cn.o
clust_cn.o: clust_cn.c clust_cn.h
	${CC} ${OPTS} -c clust_cn.c

# rule for cleaning 
clean:
	/bin/rm -f netclustwacnmdadiff netclustwacnmdadiff.o clust_cn.o


