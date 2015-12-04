#------------------------------------------------------------------------------
# Makefile for revised_cloud_mask
#------------------------------------------------------------------------------

# Set up compile options
CC = gcc
RM = rm -f
EXTRA = -Wall

# Define the include files
INC = input.h output.h revised_cloud_mask.h
INCDIR  = -I. -I$(XML2INC) -I$(ESPAINC) -I$(OPENCVINC) -I$(BOOST_INC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC = rule_based_model.c \
      buffer.c            \
      get_args.c          \
      input.c             \
      make_index.c        \
      output.c            \
      variance.c          \
      revised_cloud_mask.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
LIB   = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -L$(XML2LIB) -lxml2 \
        -L$(BOOST_LIB) -lboost_program_options \
        -L$(OPENCVLIB) -lopencv_imgproc -lopencv_core \
        -lz -lrt -lpthread -lm

# Define the executable
EXE = revised_cloud_mask

# Target for the executable
all: $(EXE)

revised_cloud_mask: $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LIB)

install:
	install -d $(BIN)
	install -m 755 $(EXE) $(BIN)
#	install -m 755 ../scripts/*.py $(BIN)

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

