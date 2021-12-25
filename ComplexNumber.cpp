#include "ComplexNumber.h"
#include <math.h>
#include <stdio.h>

//e^ix
Cmplx cmplxExp(double x) {
	return Cmplx(cos(x), sin(x));
}

void Cmplx::print() {
	if (b >= 0.0) {
		printf("%.3f + %.3fi\n", a, b);
	}
	else {
		printf("%.3f - %.3fi\n", a, -b);
	}
}

Cmplx Cmplx::operator+(const Cmplx& c) {
	return Cmplx(a + c.a, b + c.b);
}

Cmplx Cmplx::operator-(const Cmplx& c) {
	return Cmplx(a - c.a, b - c.b);
}
Cmplx Cmplx::operator*(const Cmplx& c) {
	return Cmplx(a * c.a - b * c.b, a * c.b + b * c.a);
}
Cmplx Cmplx::operator/(const Cmplx& c) {
	double denom = 1 / (c.a * c.a + c.b * c.b);
	return this->operator*(Cmplx(c.a * denom, -c.b * denom));
}

Cmplx& Cmplx::operator+=(const Cmplx& c) {
	a += c.a;
	b += c.b;
	return *this;
}

bool Cmplx::operator==(const Cmplx& c) {
	return a == c.a && b == c.b;
}

bool Cmplx::operator!=(const Cmplx& c) {
	return a != c.a || b != c.b;
}

Cmplx operator+(const Cmplx& c, const double& d) {
	return Cmplx(c.a + d, c.b);
}
Cmplx operator+(const double& d, const Cmplx& c) {
	return Cmplx(c.a + d, c.b);
}
Cmplx operator-(const Cmplx& c, const double& d) {
	return Cmplx(c.a - d, c.b);
}
Cmplx operator-(const double& d, const Cmplx& c) {
	return Cmplx(-c.a + d, -c.b);
}
Cmplx operator*(const Cmplx& c, const double& d) {
	return Cmplx(c.a * d, c.b * d);
}
Cmplx operator*(const double& d, const Cmplx& c) {
	return Cmplx(c.a * d, c.b * d);
}
Cmplx operator/(const Cmplx& c, const double& d) {
	return Cmplx(c.a / d, c.b / d);
}