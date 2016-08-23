################################
#           Makefile           #
################################

WaveAnalyzer: Analysis_40Gs.h Analysis_40Gs.cpp 
	g++ -O3 -Wall -Wextra -o Analysis_40Gs Analysis_40Gs.cpp   `root-config --cflags --glibs `   
	#
