EX_CC_FLAGS=-std=c++11 -lstdc++

fft_corr:
	gcc signal_xcorr.c \
	./kissfft/tools/kiss_fftr.c \
	./kissfft/kiss_fft.c \
	-I kissfft/ \
	-I kissfft/tools \
	-o signal_xcorr -lm

