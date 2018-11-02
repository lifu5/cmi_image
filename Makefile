BIN = ./bin
SRC = ./src
INCLUDE = ./include
OBJ = ./obj

# basic grammer:
# $@ obj;  $^ all requst;  $@< the first requst

src = $(wildcard ${SRC}/*.c ${INCLUDE}/*.c )
objects = $(patsubst %.c, ${OBJ}/%.o,$(notdir ${src}))

TARGET = myexe
BIN_TARGET = ${BIN}/${TARGET}

CXX = clang
CFLAGS = -g -Wall -I${INCLUDE}

# using .o file to generate
${BIN_TARGET}: ${objects}
	$(CXX) $(objects) -lpthread -O2 -o $@

# middle file
${OBJ}/%.o: ${SRC}/%.c
	$(CXX) $(CFLAGS) -c $< -o $@

${OBJ}/%.o: ${INCLUDE}/%.c
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -rf ${objects}
	rm -rf ${BIN_TARGET}
