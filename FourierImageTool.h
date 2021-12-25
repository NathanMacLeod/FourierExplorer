#pragma once
#include "olcPixelGameEngine.h"
#include "cstdio"
#include <iostream>    
#include <fstream>
#include <stdio.h>
#include <chrono>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <filesystem>
#include "FourierSeries.h"
#include <ctime>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <chrono>


class FourierImageTool : public olc::PixelGameEngine
{
private:
	int imageSize;
	std::vector<Cmplx> partitionPath;
	std::vector<std::vector<Cmplx>> series;
	olc::Sprite* input_image = nullptr;
	olc::Sprite* imageOG = nullptr;
	olc::Sprite* imageReconstructed = nullptr;
	olc::Sprite* fourierImage = nullptr;

	olc::Sprite* getGrayscaledAndPadded(olc::Sprite& imageIn) {
		int N = (imageIn.height > imageIn.width) ? getNextPow2(imageIn.height) : getNextPow2(imageIn.width);
		olc::Sprite* imageOut = new olc::Sprite(N, N);
		for (int x = 0; x < N; x++) {
			for (int y = 0; y < N; y++) {
				if (x >= imageIn.width || y >= imageIn.height) {
					imageOut->SetPixel(x, y, olc::BLACK);
				}
				else {
					olc::Pixel original = imageIn.GetPixel(x, y);
					int value = 0.3 * original.r + 0.59 * original.g + 0.11 * original.b;
					imageOut->SetPixel(x, y, olc::Pixel(value, value, value));
				}
			}
		}
		return imageOut;
	}

	olc::Sprite* makeFourierImage(std::vector<std::vector<Cmplx>>& series) {
		double magnitude_max = 0;
		for (int i = 0; i < series.size(); i++) {
			for (int j = 0; j < series[i].size(); j++) {
				double mag = Cmplx::abs(series[i][j]);
				if (mag > magnitude_max) {
					magnitude_max = mag;
				}
			}
		}

		magnitude_max = log(magnitude_max);

		int offset = imageSize / 2;
		olc::Sprite* image = new olc::Sprite(imageSize, imageSize);
		for (int x = 0; x < series.size(); x++) {
			for (int y = 0; y < series[x].size(); y++) {
				int value = 255 * log(Cmplx::abs(series[x][y])) / magnitude_max;
				olc::Pixel color(value, value, value);
				image->SetPixel((x + offset) % imageSize, (y + offset) % imageSize, color);
			}
		}
		return image;
	}

	//assumes image is already NxN dimensions, N power of 2, and grayscale
	void FFT_Image(olc::Sprite& imageIn) {
		int N = imageIn.width;
		if (imageIn.height != imageIn.width || imageIn.height != getNextPow2(imageIn.height)) {
			printf("INVALID CONDITIONS TO FFT_IMAGE\n");
			return;
		}

		std::vector<std::vector<Cmplx>> coeff_matrix = std::vector<std::vector<Cmplx>>(imageIn.height);
		//first FFT on each row
		for (int i = 0; i < imageIn.height; i++) {
			std::vector<Cmplx> rowData = std::vector<Cmplx>(imageIn.width);
			for (int j = 0; j < imageIn.width; j++) {
				int value = imageIn.GetPixel(j, i).r;
				rowData[j] = Cmplx(value, 0);
			}
			coeff_matrix[i] = FFT(rowData);
		}

		series = std::vector<std::vector<Cmplx>>(imageIn.width);
		int offset = imageIn.width / 2;

		std::vector<FourierTerm> significant_frequencies;

		//finaly FFT on each Column
		for (int i = 0; i < imageIn.width; i++) {
			std::vector<Cmplx> column = std::vector<Cmplx>(coeff_matrix.size());
			for (int j = 0; j < coeff_matrix.size(); j++) {
				column[j] = coeff_matrix[j][i];
			}
			series[i] = FFT(column);
			for (int j = 0; j < series[i].size(); j++) {
				significant_frequencies.push_back(FourierTerm(series[i][j], N * i + j));
			}
		}

	}

	olc::Sprite* Inverse_FFT_Image(std::vector<std::vector<Cmplx>> series) {
		int N = series.size();
		for (int i = 0; i < N; i++) {

			series[i] = InvFFT(series[i]);
		}

		olc::Sprite* outImage = new olc::Sprite(N, N);

		//finaly FFT on each Column
		for (int i = 0; i < N; i++) {
			std::vector<Cmplx> row = std::vector<Cmplx>(N);
			for (int j = 0; j < N; j++) {
				row[j] = series[j][i];
			}
			std::vector<Cmplx> values = InvFFT(row);
			for (int j = 0; j < N; j++) {
				int value = values[j].a;
				if (value > 255) {
					value = 255;
				}
				if (value < 0) {
					value = 0;
				}
				outImage->SetPixel(j, i, olc::Pixel(value, value, value));
			}
		}

		return outImage;
	}

	void deleteCoefficients(bool deleteOutside) {
		for (int y = 0; y < imageSize; y++) {
			std::vector<int> intersections;
			for (int i = 0; i < partitionPath.size(); i++) {
				Cmplx c1 = partitionPath[i];
				Cmplx c2 = (i == partitionPath.size() - 1) ? partitionPath[0] : partitionPath[i + 1];
				if ((c1.b > y && c2.b > y) || (c1.b <= y && c2.b <= y)) {
					continue; //no intersection possible
				}
				else if (c1.b == c2.b) {
					continue; //dont bother with calculating intersection if its parallel to the row
				}
				else if (c1.a == c2.a) {
					intersections.push_back(c1.a); //special case for verticle to avoid div by 0
				}
				else {
					double slope = (c2.a - c1.a) / (c2.b - c1.b);
					intersections.push_back(c1.a + (y - c1.b) * slope);
				}
			}
			std::sort(intersections.begin(), intersections.end(), [](double a, double b) {return a < b; });
			int intrsct_indx = 0; //track the number of intersections crossed
			bool outside = true;

			for (int x = 0; x < imageSize; x++) {
				while (intrsct_indx < intersections.size() && intersections[intrsct_indx] < x) {
					intrsct_indx++;
					outside = !outside;
				}

				if (outside && deleteOutside) {
					series[(x + imageSize / 2) % imageSize][(y + imageSize / 2) % imageSize] = Cmplx(0, 0);
				}
				else if (!outside && !deleteOutside) {
					series[(x + imageSize / 2) % imageSize][(y + imageSize / 2) % imageSize] = Cmplx(0, 0);
				}
			}
		}
	}

public:
	FourierImageTool(olc::Sprite* input_image)
	{
		sAppName = "FourierImageTool";
		this->input_image = input_image;
	}

	~FourierImageTool()
	{
		if (input_image != nullptr) {
			delete input_image;
		}
		if (imageOG != nullptr) {
			delete imageOG;
		}
		if (imageReconstructed != nullptr) {
			delete imageReconstructed;
		}
		if (fourierImage != nullptr) {
			delete fourierImage;
		}
	}


	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		imageOG = getGrayscaledAndPadded(*input_image);
		imageSize = imageOG->height;

		FFT_Image(*imageOG);
		fourierImage = makeFourierImage(series);
		imageReconstructed = Inverse_FFT_Image(series);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{

		Clear(olc::DARK_BLUE);
		DrawSprite(0, 0, imageOG);
		DrawSprite(imageOG->width + 10, 0, fourierImage);
		DrawSprite(imageOG->width + 10 + fourierImage->width + 10, 0, imageReconstructed);

		const int fourierImageX = imageSize + 10;
		if (GetMouse(0).bPressed) {
			int mouseX = GetMouseX();
			int mouseY = GetMouseY();

			//check in bounds
			if (mouseX >= fourierImageX && mouseX < fourierImageX + imageSize && mouseY >= 0 && mouseY < imageSize) {
				partitionPath.push_back(Cmplx(mouseX - fourierImageX, mouseY));
			}

		}

		for (int i = 0; i + 1 < partitionPath.size(); i++) {
			Cmplx c1 = partitionPath.at(i);
			Cmplx c2 = partitionPath.at(i + 1);
			DrawLine(fourierImageX + c1.a, c1.b, fourierImageX + c2.a, c2.b, olc::RED);
		}

		//Reset
		if (GetKey(olc::R).bPressed) {
			delete(fourierImage);
			delete(imageReconstructed);

			FFT_Image(*imageOG);
			fourierImage = makeFourierImage(series);
			imageReconstructed = Inverse_FFT_Image(series);
		}
		else if (partitionPath.size() >= 3 && (GetKey(olc::I).bPressed || GetKey(olc::O).bPressed)) { //partition
			deleteCoefficients(GetKey(olc::I).bPressed); //remove either outside or inside
			partitionPath.clear();

			delete(fourierImage);
			delete(imageReconstructed);

			fourierImage = makeFourierImage(series);
			imageReconstructed = Inverse_FFT_Image(series);
		}

		return true;
	}
};