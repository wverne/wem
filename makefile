# ---------------------------------------------------------------------
# wem makefile
# Author: Wesley Verne
# ---------------------------------------------------------------------

# Macros
SRC = src

# Variables

executables = wemEG wemEP
library = src/files.h src/constants.h src/constants.cpp src/EOS.h src/EOS.cpp src/Planet.h src/Planet.cpp src/stdafx.h

# Non File Targets

all: $(executables)

clean: 
	rm -f $(executables) *~

# File Targets

wemEG: src/wemEG.cpp $(library)
	g++ src/wemEG.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp -o wemEG

wemEP: src/wemEP.cpp $(library)
	g++ src/wemEP.cpp src/constants.cpp src/EOS.cpp src/Planet.cpp -o wemEP