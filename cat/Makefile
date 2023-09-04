################################################################################
# CROSS Makefile
################################################################################
PROJECT_NAME = gencat

################################################################################
# Key paths and settings
################################################################################
CC = gcc -g3 -O0 -DDEBUG
ODIR = obj
SDIR = src
BDIR = bin

CPPFLAGS = -Ihdr -MMD -MP
CFLAGS	 = -Wall
LDFLAGS  = -Llib
LDLIBS   = -Wl,-Bdynamic -lm

################################################################################
# Master inc/lib/obj/dep settings
################################################################################
EXE  = $(BDIR)/$(PROJECT_NAME)
SRC  = $(wildcard $(SDIR)/*.c)
OBJ  = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(SRC))

################################################################################
# Rules/recipes
################################################################################
.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(BDIR)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

$(ODIR)/%.o: $(SDIR)/%.c | $(ODIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BDIR) $(ODIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BDIR) $(ODIR)

-include $(OBJ:.o=.d)
