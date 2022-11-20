# OBJS = opt_vec.o

# main: $(OBJS)
# 	g++ -o $@ $(OBJS)

# %.o: %.cpp Vector3d.h opt_vec.h
# 	g++ -c $< -Wall -O3 

# OBJS = opt_vec2.o

# main: $(OBJS)
# 	g++ -o $@ $(OBJS)

# %.o: %.cpp Vector3d.h opt_vec.h
# 	g++ -c $< -Wall -O3 



# OBJS = opt_vec2.o opt_vertical.o

# HEADERS = opt_vec.h

# PREFIX = $(HOME)
# LIB_DIR = $(PREFIX)/lib
# INC_DIR = $(PREFIX)/include

# LIB = -lAndoLab_20
# LINKER_OPTS = -Wl,-R$(LIB_DIR) -L$(LIB_DIR)

# OPTS = -Wall -O3 -I$(INC_DIR)

# .PHONY: all clean

# all: main

# main: $(OBJS)
# 	g++ -o $@ $(OBJS) $(LINKER_OPTS) $(LIB)

# %.o: %.cpp $(HEADERS)
# 	g++ -c $< $(OPTS)

# clean:
# 	rm -rf *.o main



OBJS = opt_vec.o opt_vertical.o

HEADERS = opt_vec.h

PREFIX = $(HOME)
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include

LIB = -lAndoLab_20
LINKER_OPTS = -Wl,-R$(LIB_DIR) -L$(LIB_DIR)

OPTS = -Wall -O3 -I$(INC_DIR)

.PHONY: all clean

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS) $(LINKER_OPTS) $(LIB)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

clean:
	rm -rf *.o main