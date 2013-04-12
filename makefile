# ---------------------------------------------------------------------
# wem makefile
# Author: Wesley Verne
# ---------------------------------------------------------------------

# Macros
SRC = src

# Variables

executables = wemEG wemEP
debug_executables = wemEG_debug wemEP_debug
library = src/files.h src/constants.h src/constants.cpp src/EOS.h src/EOS.cpp src/Planet.h src/Planet.cpp src/stdafx.h

# Non File Targets

all: $(executables)

debug: $(debug_executables)

clean: 
	rm -f $(executables) *_debug *~

# File Targets

wemEG: src/wemEG.cpp $(library)
	g++ -O3 src/wemEG.cpp src/constants.cpp src/EOS.cpp -o wemEG

wemEP: src/wemEP.cpp $(library)
	g++ -O3 src/wemEP.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp -o wemEP

wemEG_debug: src/wemEG.cpp $(library)
	g++ -g src/wemEG.cpp src/constants.cpp src/EOS.cpp -o wemEG_debug

wemEP_debug: src/wemEP.cpp $(library)
	g++ -g src/wemEP.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp -o wemEP_debug