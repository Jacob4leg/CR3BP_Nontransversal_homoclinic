# Executables
PROGS = main

# Directory with capd-config; empty means capd-config is taken from PATH.
# Set it on the command line: make CAPDBINDIR=/path/to/CAPD/build/bin
CAPDBINDIR ?=
CAPDCONFIG = $(if $(CAPDBINDIR),$(CAPDBINDIR)/)capd-config

# Directories
SRCDIR = src
INCLUDEDIR = include
OBJDIR = .obj

# Source and header files
SOURCES = $(filter-out $(SRCDIR)/main.cpp, $(wildcard $(SRCDIR)/*.cpp))
HEADERS = $(wildcard $(INCLUDEDIR)/*.h)

# Object files (without main.cpp)
OBJ_FILES = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Compiler and flags
CAPDCXX := $(shell $(CAPDCONFIG) --variable=capd_cxx)
CAPDFLAGS = `$(CAPDCONFIG) --cflags`
CAPDLIBS = `$(CAPDCONFIG) --libs`
CXXFLAGS += $(CAPDFLAGS) -I$(INCLUDEDIR)

# Main rule
.PHONY: all
all: $(PROGS)

# Linking executables
$(PROGS): % : $(OBJDIR)/%.o $(OBJ_FILES)
	$(CAPDCXX) -o $@ $< $(OBJ_FILES) $(CAPDLIBS)

# Compiling sources to objects
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CAPDCXX) $(CXXFLAGS) -MT $@ -MD -MP -MF $(@:.o=.d) -c -o $@ $<

# Dependency files
-include $(OBJ_FILES:.o=.d)

# Cleaning
.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(PROGS)
