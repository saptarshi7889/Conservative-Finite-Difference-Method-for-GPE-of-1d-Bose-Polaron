CXX ?= c++
EIGEN_INCLUDE ?= /usr/include/eigen3
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
PREFIX ?= /usr/local

.PHONY: all test install clean
all: build/gpe
build/gpe: src/main.cpp src/solver.cpp src/continuation.cpp src/numerics.hpp src/snapshot_archive.hpp include/gpe/solver.hpp include/gpe/GPE_solver_contd.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) -I"$(EIGEN_INCLUDE)" $(CXXFLAGS) -Iinclude src/main.cpp src/solver.cpp src/continuation.cpp $(LDFLAGS) -lz -larchive -o $@
build/test_numerics: tests/numerics.cpp src/numerics.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) -I"$(EIGEN_INCLUDE)" $(CXXFLAGS) tests/numerics.cpp $(LDFLAGS) -o $@
build/test_library: tests/library.cpp src/solver.cpp src/continuation.cpp src/numerics.hpp src/snapshot_archive.hpp include/gpe/solver.hpp include/gpe/GPE_solver_contd.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) -I"$(EIGEN_INCLUDE)" -Iinclude $(CXXFLAGS) tests/library.cpp src/solver.cpp src/continuation.cpp $(LDFLAGS) -lz -larchive -o $@
build/test_continuation: tests/continuation.cpp src/solver.cpp src/continuation.cpp src/numerics.hpp src/snapshot_archive.hpp include/gpe/solver.hpp include/gpe/GPE_solver_contd.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) -I"$(EIGEN_INCLUDE)" -Iinclude $(CXXFLAGS) tests/continuation.cpp src/solver.cpp src/continuation.cpp $(LDFLAGS) -lz -larchive -o $@
test: build/gpe build/test_numerics build/test_library build/test_continuation
	./build/test_numerics
	python3 tests/library.py ./build/test_library ./build/gpe
	python3 tests/integration.py ./build/gpe
	python3 tests/continuation.py ./build/test_continuation
install: build/gpe
	install -Dm755 build/gpe "$(DESTDIR)$(PREFIX)/bin/gpe"
clean:
	rm -rf build
