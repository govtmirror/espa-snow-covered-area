#-----------------------------------------------------------------------------
# Makefile
#
# Project Name: snow-covered-area
# Makefile that will invoke subsequent Makefiles in subdirectories.
#-----------------------------------------------------------------------------
.PHONY: check-environment all install clean

include make.config

#-----------------------------------------------------------------------------
all:
	echo "make all in revised_cm"; \
        (cd revised_cm; $(MAKE) all);

#-----------------------------------------------------------------------------
install: check-environment
	echo "make install in revised_cm"; \
        (cd revised_cm; $(MAKE) install);

#-----------------------------------------------------------------------------
clean:
	echo "make clean in revised_cm"; \
        (cd revised_cm; $(MAKE) clean);

#-----------------------------------------------------------------------------
check-environment:
ifndef PREFIX
    $(error Environment variable PREFIX is not defined)
endif

