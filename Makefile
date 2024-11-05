release_exec = build/generate_network
debug_exec = build.debug/generate_network
release_obj_files = build/format.o build/generate_network.o
debug_obj_files = build.debug/format.o build.debug/generate_network.o
src_files = ext/fmt/src/format.cc 
DEFINES=-DGEN_NETWORK_VERSION=1.0
INCLUDES=-I./include -I./ext/json/include -I./ext/fmt/include -I./ext/mio/include -I./ext/CLI11/include 
CARGS=-Wall -Wno-maybe-uninitialized -fPIC -O2 -g -DNDEBUG -ffast-math -funroll-loops -std=gnu++17 -fopenmp
DEBUG_CARGS=-Wall -Wno-maybe-uninitialized -fPIC -O2 -g -DNDEBUG -ffast-math -funroll-loops -std=gnu++17 -fopenmp
LINK_ARGS=-Wall -Wno-maybe-uninitialized -fPIC -O2 -g -DNDEBUG -ffast-math -funroll-loops 
LINK_LIBS=-ldl -lpthread
BUILD_DIR=build
DEBUG_BUILD_DIR=build.debug

all: release

doxy:
	doxygen ./docs/Doxyfile

clean: 
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/generate_network
	rm -rf $(DEBUG_BUILD_DIR)/*.o $(DEBUG_BUILD_DIR)/generate_network

generate_network_debug: src/generate_network.cpp ext/fmt/src/format.cc include/*.hpp
	mkdir -p $(DEBUG_BUILD_DIR)/
	g++ $(DEFINES) $(INCLUDES) $(DEBUG_CARGS) -o $(DEBUG_BUILD_DIR)/format.o -c ./ext/fmt/src/format.cc
	g++ $(DEFINES) $(INCLUDES) $(DEBUG_CARGS) -o $(DEBUG_BUILD_DIR)/generate_network.o -c ./src/generate_network.cpp
	g++ $(LINK_ARGS) $(DEBUG_BUILD_DIR)/generate_network.o $(DEBUG_BUILD_DIR)/format.o -o $(DEBUG_BUILD_DIR)/generate_network $(LINK_LIBS)

generate_network: src/generate_network.cpp ext/fmt/src/format.cc include/*.hpp
	mkdir -p $(BUILD_DIR)/
	g++ $(DEFINES) $(INCLUDES) $(CARGS) -o $(BUILD_DIR)/format.o -c ./ext/fmt/src/format.cc
	g++ $(DEFINES) $(INCLUDES) $(CARGS) -o $(BUILD_DIR)/generate_network.o -c ./src/generate_network.cpp
	g++ $(LINK_ARGS) $(BUILD_DIR)/generate_network.o $(BUILD_DIR)/format.o -o $(BUILD_DIR)/generate_network $(LINK_LIBS)


release: generate_network

debug: generate_network_debug

# constructs random network based on network.cfg file parameters
network: release
	mkdir -p output
	./build/generate_network -c ./config/mriNetworkNew.cfg -o output/random_connection_info.txt
