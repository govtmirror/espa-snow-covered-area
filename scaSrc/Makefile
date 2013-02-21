#------------------------------------------------------------------------------
# Makefile for sample HDF code
#------------------------------------------------------------------------------

# Set the HDF4 environment variables
#HDFINC=/home/gschmidt/hdf-4.2.5/include
#HDFLIB=/home/gschmidt/hdf-4.2.5/lib

# Set the HDF_EOS environment variables
#HDFEOS_INC=/home/gschmidt/hdfeos2.17v1.00/hdfeos/include
#HDFEOS_LIB=/home/gschmidt/hdfeos2.17v1.00/hdfeos/lib

# Set the GCTP environment variables (as delivered with HDF-EOS)
#HDFEOS_GCTPINC=/home/gschmidt/hdfeos2.17v1.00/hdfeos/gctp/include
#HDFEOS_GCTPLIB=/home/gschmidt/hdfeos2.17v1.00/hdfeos/lib

# Set the SZIP environment variables
#SZIPINC=/home/gschmidt/szip-2.1/include
#SZIPLIB=/home/gschmidt/szip-2.1/lib

# Set the ZLIB environment variables
#ZLIBINC=/home/gschmidt/zlib-1.2.7/include
#ZLIBLIB=/home/gschmidt/zlib-1.2.7/lib

# Set the JPEG environment variables
#JPEGINC=/home/gschmidt/jpeg-6b/include
#JPEGLIB=/home/gschmidt/jpeg-6b/lib

# Set up compile options
CC = gcc
RM = rm -f
EXTRA = -g -O -Wall

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
LIB   = -static -L$(HDFLIB) -lmfhdf -ldf -L$(ZLIBLIB) -lz -L$(SZIPLIB) -lsz \
        -L$(JPEGLIB) -ljpeg -L$(HDFLIB) -lxdr
EOSLIB = -static -L$(HDFEOS_LIB) -lhdfeos -L$(HDFEOS_GCTPLIB) -lGctp -lm

# Define the executable
EXE = scene_based_sca

# Target for the executable
all: $(EXE)

scene_based_sca: $(OBJ) $(INC)
	$(CC) -o $(EXE) $(OBJ) $(LIB) $(EOSLIB)

install:
	cp $(EXE) $(BIN)

clean:
	$(RM) *.o

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

