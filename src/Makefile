#------------------------------------------------------------------------------
# Makefile for sample HDF code
#------------------------------------------------------------------------------

# Set up compile options
CC = gcc
RM = rm -f
EXTRA = -Wall -g

# Define the include files
INC = bool.h const.h date.h error_handler.h input.h myhdf.h mystring.h \
output.h space.h sca.h
INCDIR  = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC) -I$(JPEGINC) -I$(SZIPINC) -I$(ZLIBINC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC = cloud_cover_class.c \
      combine_qa.c        \
      date.c              \
      error_handler.c     \
      get_args.c          \
      input.c             \
      myhdf.c             \
      mystring.c          \
      output.c            \
      qa_mask.c           \
      shaded_relief.c     \
      snow_cover_class.c  \
      space.c             \
      write_envi_hdr.c    \
      scene_based_sca.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
LIB   = -L$(HDFLIB) -lmfhdf -ldf -lxdr -L$(JPEGLIB) -ljpeg \
        -L$(ZLIBLIB) -lz -L$(SZIPLIB) -lsz -lm
##LIB   = -L$(HDFLIB) -lmfhdf -ldf -lxdr -L$(JPEGLIB) -ljpeg \
##        -L$(ZLIBLIB) -lz -L$(SZIPLIB)/libsz.a -lm
EOSLIB = -L$(HDFEOS_LIB) -lhdfeos -L$(HDFEOS_GCTPLIB) -lGctp

# Define the executable
EXE = scene_based_sca

# Target for the executable
all: $(EXE)

scene_based_sca: $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(EOSLIB) $(LIB)

install:
	cp $(EXE) $(BIN)

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

