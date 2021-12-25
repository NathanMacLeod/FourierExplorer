#pragma once
#include <vector>
#include "ComplexNumber.h"

struct FourierTerm {
	FourierTerm(Cmplx coeff, int freq) {
		this->coeff = coeff;
		this->freq = freq;
	}
	FourierTerm() {
		coeff = Cmplx(0, 0);
		freq = 0;
	}
	Cmplx coeff;
	int freq;
};

int getNextPow2(int n, int* pow=nullptr);

std::vector<FourierTerm> PointFourierTransform(std::vector<Cmplx>& samples, int startK, int endK);
std::vector<FourierTerm> LerpFourierTransform(std::vector<Cmplx>& path, int startK, int endK);

std::vector<Cmplx> DFT(std::vector<Cmplx>& samples);
std::vector<Cmplx> InvDFT(std::vector<Cmplx>& series);
std::vector<Cmplx> FFT(std::vector<Cmplx> samples);
std::vector<Cmplx> InvFFT(std::vector<Cmplx> samples);

Cmplx InterrogateSeries(std::vector<Cmplx>& series, int startIndx, double x, double period);