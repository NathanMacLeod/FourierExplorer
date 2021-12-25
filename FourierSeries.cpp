#include "FourierSeries.h"
#define PI 3.1415926535897932

//calcualtes fourier series on discrete data, number of terms is parameter
std::vector<FourierTerm> LerpFourierTransform(std::vector<Cmplx>& path, int startK, int endK) {
	std::vector<FourierTerm> outVector = std::vector<FourierTerm>(endK - startK + 1);
	if (path.size() <= 1) {
		for (int i = startK; i <= endK; i++) {
			outVector.at(i - startK) = FourierTerm(Cmplx(0, 0), i);
		}
	}

	double total_length = 0;
	for (int i = 0; i < path.size(); i++) {
		Cmplx p1 = path.at(i);
		Cmplx p2 = (i == path.size() - 1) ? path.at(0) : path.at(i + 1);

		total_length += Cmplx::abs(p2 - p1);
	}

	for (int i = startK; i <= endK; i++) {

		Cmplx sum(0, 0);
		double t = 0;

		
		for (int j = 0; j < path.size(); j++) {
			Cmplx p1 = path.at(j);
			Cmplx p2 = (j == path.size() - 1) ? path.at(0) : path.at(j + 1);
			double dt = Cmplx::abs(p2 - p1) / total_length;

			if (i == 0) {
				sum += (p2 + p1) * dt/2.0;
			}
			else {
				sum += Cmplx(0, 1 / (2 * PI * i)) * (p2 * cmplxExp(-2 * PI * i * (t + dt)) - p1 * cmplxExp(-2 * PI * i * t));
				sum += (p2 - p1) * Cmplx(1 / (dt * 4 * PI * PI * i * i), 0) * (cmplxExp(-2 * PI * i * (t + dt)) - cmplxExp(-2 * PI * i * t));
				//sum += Cmplx(1 / (dt * 4 * PI * PI * i * i), 1 / (2 * PI * i)) * (p2 * cmplxExp(-2 * PI * i * (t + dt)) - p1 * cmplxExp(-2 * PI * i * t));
			}
			t += dt;
		}

		/*Cmplx sum(0, 0);
		double dt = 1/ (float)fidelity;
		for (double t = 0; t <= 1.0; t += dt) {

			int data_curr = ((int)(t * samples.size())) % samples.size();
			int data_next = (data_curr + 1) % samples.size();
			double interp_t = fmod(t * samples.size(), 1);

			Cmplx data_point = dt * samples.size() * (interp_t * samples.at(data_next) + (1 - interp_t) * samples.at(data_curr));
			sum += data_point * cmplxExp(-2 * PI * i * t);
		}*/

		outVector.at(i - startK) = FourierTerm(sum, i);
	}
	return outVector;
}

//calcualtes fourier series on discrete data, number of terms is parameter
std::vector<FourierTerm> PointFourierTransform(std::vector<Cmplx>& samples, int startK, int endK) {
	std::vector<FourierTerm> outVector = std::vector<FourierTerm>(endK - startK  + 1); 
	
	for (int i = startK; i <= endK; i++) {
		Cmplx sum(0, 0);
		for (int j = 0; j < samples.size(); j++) {
			if (i == 0) {
				sum += samples[j];
			}
			else {
				sum += samples[j] * Cmplx(0, samples.size() / (2 * i * PI)) * (cmplxExp(-2 * PI * i * (j + 1) / samples.size()) - cmplxExp(-2 * PI * i * (j) / samples.size()));
			}
		}

		/*Cmplx sum(0, 0);
		double dt = 1/ (float)fidelity;
		for (double t = 0; t <= 1.0; t += dt) {

			int data_curr = ((int)(t * samples.size())) % samples.size();
			int data_next = (data_curr + 1) % samples.size();
			double interp_t = fmod(t * samples.size(), 1);

			Cmplx data_point = dt * samples.size() * (interp_t * samples.at(data_next) + (1 - interp_t) * samples.at(data_curr));
			sum += data_point * cmplxExp(-2 * PI * i * t);
		}*/
		
		outVector.at(i - startK) = FourierTerm(sum / samples.size(), i);
	}
	return outVector;
}

int getNextPow2(int n, int* pow) {
	int pow2 = 1;
	int exp = 0;
	while (n > pow2) {
		pow2 <<= 0x1;
		exp++;
	}
	if (pow != nullptr) {
		*pow = exp;
	}
	return pow2;
}

std::vector<Cmplx> FFT(std::vector<Cmplx> samples) {
	int pow;
	int N = getNextPow2(samples.size(), &pow);
	while (samples.size() < N) {
		samples.push_back(Cmplx(0, 0));
	}

	std::vector<Cmplx> v1 = std::vector<Cmplx>(N);
	std::vector<Cmplx> v2 = std::vector<Cmplx>(N);

	//partition into set of even indexes, odd indexes (even on top, odd on bottom), then repartition those (reoccours logn times)
	for (int i = 0; i < N; i++) {
		int set_i = 0;
		int rel_indx = i;
		for (int j = 0; j < pow - 1; j++) {
			set_i *= 2;
			if (rel_indx % 2 != 0) {
				set_i++;
			}
			rel_indx /= 2;
		}
		v1[2 * set_i + rel_indx] = samples[i];
	}

	std::vector<Cmplx>* curr = &v1;
	std::vector<Cmplx>* next = &v2;

	//F'[K] = F_e[K] + (w_n)^k * F_o[K]
	int partition_size = 1;
	while (partition_size < N) {
		for (int offst = 0; offst < N; offst += partition_size * 2) {
			for (int k = 0; k < partition_size * 2; k++) {
				next->at(offst + k) = curr->at(offst + (k % partition_size)) + curr->at(offst + partition_size + (k % partition_size)) * cmplxExp(-2 * PI * k / (partition_size * 2));
			}
		}
		partition_size *= 2;
		std::vector<Cmplx>* tmp = curr;
		curr = next;
		next = tmp; //(old values will be overwritten);
	}

	std::vector<Cmplx> outVector = std::vector<Cmplx>(samples.size());
	for (int i = 0; i < outVector.size(); i++) {
		outVector[i] = curr->at(i);
	}

	return outVector;
}

std::vector<Cmplx> InvFFT(std::vector<Cmplx> samples) {
	int pow;
	int N = getNextPow2(samples.size(), &pow);
	while (samples.size() < N) {
		samples.push_back(Cmplx(0, 0));
	}

	std::vector<Cmplx> v1 = std::vector<Cmplx>(N);
	std::vector<Cmplx> v2 = std::vector<Cmplx>(N);

	//partition into set of even indexes, odd indexes (even on top, odd on bottom), then repartition those (reoccours logn times)
	for (int i = 0; i < N; i++) {
		int set_i = 0;
		int rel_indx = i;
		for (int j = 0; j < pow - 1; j++) {
			set_i *= 2;
			if (rel_indx % 2 != 0) {
				set_i++;
			}
			rel_indx /= 2;
		}
		v1[2 * set_i + rel_indx] = samples[i];
	}

	std::vector<Cmplx>* curr = &v1;
	std::vector<Cmplx>* next = &v2;

	//F'[K] = F_e[K] + (w_n)^k * F_o[K]
	int partition_size = 1;
	while (partition_size < N) {
		for (int offst = 0; offst < N; offst += partition_size * 2) {
			for (int k = 0; k < partition_size * 2; k++) {
				next->at(offst + k) = curr->at(offst + (k % partition_size)) + curr->at(offst + partition_size + (k % partition_size)) * cmplxExp(2 * PI * k / (partition_size * 2));
			}
		}
		partition_size *= 2;
		std::vector<Cmplx>* tmp = curr;
		curr = next;
		next = tmp; //(old values will be overwritten);
	}

	std::vector<Cmplx> outVector = std::vector<Cmplx>(samples.size());
	double scale = 1.0 / N;
	for (int i = 0; i < outVector.size(); i++) {
		outVector[i] = curr->at(i) * scale;
	}

	return outVector;
}

std::vector<Cmplx> DFT(std::vector<Cmplx>& samples) {
	std::vector<Cmplx> outVector = std::vector<Cmplx>(samples.size());
	int N = samples.size();
	for (int i = 0; i < N; i++) {
		Cmplx sum(0, 0);
		for (int j = 0; j < N; j++) {
			sum += samples[j] * cmplxExp(-2 * PI * i * j / N);
		}
		outVector.at(i) = sum;
	}
	return outVector;
}

std::vector<Cmplx> InvDFT(std::vector<Cmplx>& series) {
	std::vector<Cmplx> outVector = std::vector<Cmplx>(series.size());
	int N = series.size();
	for (int i = 0; i < N; i++) {
		Cmplx sum(0, 0);
		for (int j = 0; j < N; j++) {
			sum += series[j] * cmplxExp(2 * PI * i * j / N);
		}
		outVector.at(i) = sum / N;
	}
	return outVector;
}

Cmplx InterrogateSeries(std::vector<Cmplx>& series, int startIndx, double x, double period) {
	int N = series.size();
	Cmplx sum(0, 0);
	for (int j = 0; j < N; j++) {
		sum += series[j] * cmplxExp(2 * PI * x * (startIndx + j) / period);
	}
	return sum / period;
}