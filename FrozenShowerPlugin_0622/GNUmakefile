name := FrozenShowerPlugin
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib 

CPPFLAGS += -I../PluginManager/include -I../../Kernel/include -I../../Geometry/CGA/include
CPPFLAGS += $(shell root-config --cflags)

ifdef G4ANALYSIS_USE
  CPPFLAGS += $(shell aida-config --include)
endif

include ../../Kernel/GNUmakefile.common
include $(G4INSTALL)/config/binmake.gmk
