CXX ?= c++
EIGEN_INCLUDE ?= /usr/include/eigen3
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
PREFIX ?= /usr/local

.PHONY: all test install clean
all: build/gpe
build/gpe: src/main.cpp src/numerics.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) -I"$(EIGEN_INCLUDE)" $(CXXFLAGS) src/main.cpp $(LDFLAGS) -o $@
test: build/gpe
	python3 tests/integration.py ./build/gpe
install: build/gpe
	install -Dm755 build/gpe "$(DESTDIR)$(PREFIX)/bin/gpe"
clean:
	rm -rf build
