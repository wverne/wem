# ---------------------------------------------------------------------
# wem makefile
# Author: Wesley Verne
# ---------------------------------------------------------------------

# Macros
SRC = src

# Variables

executables = wemEG wemEP wemMR
debug_executables = wemEG_debug wemEP_debug wemMR_debug
test_executables = planetTest
library = src/stdafx.h src/files.h src/constants.h src/constants.cpp src/EOS.h src/EOS.cpp src/Planet.h src/Planet.cpp src/PlanetComp.h src/PlanetComp.cpp

# Non File Targets

all: $(executables)

debug: $(debug_executables)

test: $(test_executables)

clean: 
	rm -f $(executables) $(debug_executables) $(test_executables)

# File Targets

wemEG: src/wemEG.cpp $(library)
	g++ -O3 src/wemEG.cpp src/constants.cpp src/EOS.cpp -o wemEG

wemEP: src/wemEP.cpp $(library)
	g++ -O3 src/wemEP.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp src/PlanetComp.cpp -o wemEP

wemMR: src/wemMR.cpp $(library)
	g++ -O3 src/wemMR.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp src/PlanetComp.cpp -o wemMR

wemEG_debug: src/wemEG.cpp $(library)
	g++ -g src/wemEG.cpp src/constants.cpp src/EOS.cpp -o wemEG_debug

wemEP_debug: src/wemEP.cpp $(library)
	g++ -g src/wemEP.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp src/PlanetComp.cpp -o wemEP_debug

wemMR_debug: src/wemMR.cpp $(library)
	g++ -g src/wemMR.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp src/PlanetComp.cpp -o wemMR_debug

planetTest: src/PlanetTest.cpp $(library)
	g++ -g src/PlanetTest.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp -o planetTest
