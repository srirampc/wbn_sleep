release_exec = build/generate_network
debug_exec = build.debug/generate_network
objects = build/format.o build/generate_network.o
debug_objects = build.debug/format.o build.debug/generate_network.o
src_files = ext/fmt/src/format.cc

# Compilers
ifeq ($(origin CC),default)
CC = gcc
endif
ifeq ($(origin CXX),default)
CXX = g++
endif

#
DEFINES=-DGEN_NETWORK_VERSION=1.0
INCLUDES=-I./include -I./ext/json/include -I./ext/fmt/include -I./ext/mio/include -I./ext/CLI11/include 
CPPFLAGS=-Wall -Wno-maybe-uninitialized -fPIC -O2 -g -DNDEBUG -ffast-math -funroll-loops -std=gnu++17 -fopenmp
LINK_ARGS=-Wall -Wno-maybe-uninitialized -fPIC -O2 -g -DNDEBUG -ffast-math -funroll-loops 
LINK_LIBS=-ldl -lpthread -lgomp
BUILD_DIR=build
DEBUG_BUILD_DIR=build.debug

#
debug: CPPFLAGS=-Wall -Wno-maybe-uninitialized -fPIC -O2 -g -DNDEBUG -ffast-math -funroll-loops -std=gnu++17 -fopenmp
debug: BUILD_DIR=$(DEBUG_BUILD_DIR)
debug: objects=$(debug_objects)


all: generate_network
# Note: PHONY is important here. Without it, implicit rules will try to build the executable "all", since the prereqs are ".o" files.
.PHONY: all 

# Define a pattern rule that compiles every .c file into a .o file
$(objects) : $(BUILD_DIR)/%.o : src/%.cpp
		$(CXX) -c $(DEFINES) $(INCLUDES) $(CPPFLAGS) $< -o $@

doxy:
	doxygen ./docs/Doxyfile

clean: 
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/generate_network
	rm -rf $(DEBUG_BUILD_DIR)/*.o $(DEBUG_BUILD_DIR)/generate_network

$(BUILD_DIR):
	ln -s ../ext/fmt/src/format.cc src/format.cpp
	mkdir -p $(BUILD_DIR)/

generate_network: $(BUILD_DIR) $(objects)
#	mkdir -p $(BUILD_DIR)/
#	$(CXX) $(DEFINES) $(INCLUDES) $(CPPFLAGS) -o $(BUILD_DIR)/format.o -c ./ext/fmt/src/format.cc
#	$(CXX) $(DEFINES) $(INCLUDES) $(CPPFLAGS) -o $(BUILD_DIR)/generate_network.o -c ./src/generate_network.cpp
	$(CXX) $(LINK_ARGS) $(BUILD_DIR)/generate_network.o $(BUILD_DIR)/format.o -o $(BUILD_DIR)/generate_network $(LINK_LIBS)

release: generate_network

debug: generate_network

# constructs random network based on config/mriNetworkNew.cfg file parameters
network: release
	mkdir -p output
	./build/generate_network -c ./config/mriNetworkNew.cfg -d ./config/data_config.json -s output/connect_summary.txt -o output/random_connection_info.txt
