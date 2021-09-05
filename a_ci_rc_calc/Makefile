# Opsti makefile za shemu direktoriuma
# /--src/
#  |
#  *-obj/
#  |
#  --bin/
#
###############

SUBMITNAME=rc_with_aci_curve_V1.tar.gz

SRCDIR=src
BINDIR=bin
OBJDIR=obj

executable=$(BINDIR)/aci

SOURCES=$(wildcard $(SRCDIR)/*.c)
OBJECTS=$(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SOURCES))

CC=gcc
CFLAGS=-g -Wall
LDFLAGS=-lm

all: $(executable)

$(executable): $(OBJECTS)
	$(CC) -o $@ $(CFLAGS) $(OBJECTS) $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(LDFLAGS)

archive: clean
	$(eval DATENOW:=$(shell date +'%Y%m%d%H%M%S'))
	$(eval TARGZNAME:=$(subst tar.gz,$(DATENOW).tar.gz,$(SUBMITNAME)))
	@echo $(TARGZNAME)
	tar czvf ../$(TARGZNAME) .

clean:
	rm -fr $(BINDIR)/* $(OBJDIR)/*
