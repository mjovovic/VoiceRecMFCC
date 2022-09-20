
default: compile run

compile:
	#mora export DISPLAY=:0 u terminalu
	g++ main.c -o main -lm -lsndfile -I/usr/include/python2.7 -lpython2.7  -I kissfft kissfft/kiss_fft.c -I DTW_cpp/include/DTW_hpp
run:
	./main

clean:
	rm main