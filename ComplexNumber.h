#pragma once
#include <math.h>

class Cmplx {
public:
	Cmplx() {
		a = 0;
		b = 0;
	}
	Cmplx(double a, double b) {
		this->a = a;
		this->b = b;
	}
	void print();
	Cmplx operator+(const Cmplx& c);
	Cmplx operator-(const Cmplx& c);
	Cmplx operator*(const Cmplx& c);
	Cmplx operator/(const Cmplx& c);
	Cmplx& operator+=(const Cmplx& c);
	bool operator==(const Cmplx& c);
	bool operator!=(const Cmplx& c);
	
	static double abs(const Cmplx& c) {
		return sqrt(c.a * c.a + c.b * c.b);
	}

	double a;
	double b;
};

Cmplx cmplxExp(double x); //e^ix

Cmplx operator+(const Cmplx& c, const double& d);
Cmplx operator+(const double& d, const Cmplx& c);
Cmplx operator-(const Cmplx& c, const double& d);
Cmplx operator-(const double& d, const Cmplx& c);
Cmplx operator*(const Cmplx& c, const double& d);
Cmplx operator*(const double& d, const Cmplx& c);
Cmplx operator/(const Cmplx& c, const double& d);