#    Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it)
LIBRARY_DIR=$(HOME)/itensor

APP1=brickwork
APP2=cluster

HEADERS=src/clustermap.h

CCFILES=src/$(APP1).cc src/$(APP2).cc src/clustermap.cc

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	@$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	@$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

all: $(APP1) $(APP2)
debug: $(APP1)-g $(APP2)-g

$(APP1): $(OBJECTS) $(ITENSOR_LIBS)
	@$(CCCOM) $(CCFLAGS) src/$(APP1).o src/clustermap.o -o $(APP1) $(LIBFLAGS)

#$(APP1)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
#	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP1)-g $(LIBGFLAGS)

$(APP2): $(OBJECTS) $(ITENSOR_LIBS)
	@$(CCCOM) $(CCFLAGS) src/$(APP2).o src/clustermap.o -o $(APP2) $(LIBFLAGS)

#$(APP2)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
#	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP2)-g $(LIBGFLAGS)

clean:
	@rm -fr src/.debug_objs src/*.o $(APP1) $(APP1)-g $(APP2) $(APP2)-g

mkdebugdir:
	@mkdir -p .debug_objs

