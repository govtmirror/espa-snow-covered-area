#-----------------------------------------------------------------------------
# Makefile
#
# For building the revised cloud mask
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC = gcc
RM = rm
EXTRA_OPTIONS = $(static_option) $(threading_options) $(profiling_options)
EXTRA   = -Wall -O2 $(EXTRA_OPTIONS)

# Define the include files
INC = input.h output.h revised_cloud_mask.h

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

# Define include paths
INCDIR  = -I. -I$(XML2INC) -I$(ESPAINC) -I$(OPENCVINC) -I$(BOOST_INC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the object libraries and paths
EXLIB   = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
          -L$(XML2LIB) -lxml2 \
          -L$(JBIGLIB) -ljbig \
          -L$(BOOST_LIB) -lboost_program_options \
          -L$(OPENCVLIB) -lopencv_imgproc -lopencv_core \
          -lz -lrt -lpthread -lstdc++
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define C executable
EXE = revised_cloud_mask

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

#-----------------------------------------------------------------------------
install: $(EXE)
	install -d $(link_path)
	install -d $(bin_install_path)
	install -m 755 $(EXE) $(bin_install_path)
	ln -sf $(link_source_path)/$(EXE) $(link_path)/$(EXE)

#-----------------------------------------------------------------------------
clean:
	$(RM) -f *.o $(EXE)

#-----------------------------------------------------------------------------
$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

