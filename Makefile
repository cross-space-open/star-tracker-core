################################################################################
# CROSS Core Makefile
################################################################################
PROJECT_NAME = star_tracker_core

################################################################################
# Key paths and settings
################################################################################
CC = gcc -g3 -O0 -DDEBUG
ODIR = obj
SDIR = src
LSDIR = libsrc/ica
BDIR = bin
HDIR = hdr

CPPFLAGS = -Ihdr -MMD -MP -I$(LSDIR)
CFLAGS	 = -Wall -D_FILE_OFFSET_BITS=64
LDFLAGS  = -Llib
LDLIBS   = -Wl,-Bdynamic -lm

################################################################################
# Master inc/lib/obj/dep settings
################################################################################
EXE  = $(BDIR)/$(PROJECT_NAME)
SRC  = $(wildcard $(SDIR)/*.c)
LSRC = $(wildcard $(LSDIR)/*.c)
OBJ  = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(SRC))
OBJ += $(patsubst $(LSDIR)/%.c, $(ODIR)/%.o, $(LSRC))

################################################################################
# Rules/recipes
################################################################################
.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(BDIR)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

$(ODIR)/%.o: $(SDIR)/%.c | $(ODIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(ODIR)/%.o: $(LSDIR)/%.c | $(ODIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BDIR) $(ODIR):
	mkdir -p $@

clean:
	-del -f $(BDIR) $(ODIR)

-include $(OBJ:.o=.d)
