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
#include "FourierVisualizer.h"
#include <ctime>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <chrono>

#define PI 3.1415926535897932

class FourierVisualizer : public olc::PixelGameEngine
{
private:
	std::vector<FourierTerm> series;
	std::vector<Cmplx> data;
	std::vector<Cmplx> path;

	enum Mode { DisplayDFT, DisplaySeries, DrawFunction, InputDFTData };
	Mode mode;

	float t = 0;
	double timeFactor = 1;
	bool followTip = false;
	bool drawCircles = false;
	bool drawFunction = true;
	bool drawPath = true;
	double scale = 400;
	int N_terms = 1;
	double cameraX = 0;
	double cameraY = 0;
	int mouseX;
	int mouseY;

public:
	FourierVisualizer(bool DFT)
	{
		sAppName = "FourierVisualizer";
		if (DFT) {
			mode = InputDFTData;
		}
		else {
			mode = DrawFunction;
		}
	}

	void calculatePath(int resolution) {
		double t = 0;
		double dt = 1.0 / resolution;
		path = std::vector<Cmplx>(resolution);
		for (int i = 0; i < resolution; i++) {
			Cmplx seriesTip(0, 0);
			for (int j = 0; j < series.size(); j++) {
				seriesTip += series.at(j).coeff * cmplxExp(2 * 3.141592 * series.at(j).freq * t);
			}
			path[i] = seriesTip;
			t += dt;
		}
	}

	void renderPath(int offsetX, int offsetY, double scale, double currT, int fadeExp, olc::Pixel color) {
		//int startIndx = (int) (fmod(currT, 1.0) * path.size());
		for (int i = 0; i < path.size(); i++) {
			Cmplx p1 = path[i];
			Cmplx p2 = (i == path.size() - 1) ? path[0] : path[i + 1];
			DrawLine(offsetX + p1.a * scale, offsetY + p1.b * scale, offsetX + p2.a * scale, offsetY + p2.b * scale, color);
		}
	}

	void updateSeries(int N) {
		series.clear();
		series = LerpFourierTransform(data, -N, N);
		std::sort(series.begin(), series.end(), [](FourierTerm f1, FourierTerm f2) {
			if (abs(f1.freq) < abs(f2.freq)) {
				return true;
			}
			else if (abs(f1.freq) > abs(f2.freq)) {
				return false;
			}
			else {
				return f1.freq > f2.freq;
			}
			});
		calculatePath(10000);
	}

	bool OnUserCreate() override
	{
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{

		t += fElapsedTime * timeFactor;

		int centerX = ScreenWidth() / 2;
		int centerY = ScreenHeight() / 2;

		///////////////////////
		//////USER INPUT///////
		///////////////////////
		if (mode == DrawFunction || mode == InputDFTData) {
			if (GetMouse(0).bPressed) {
				int mouseX = GetMouseX();
				int mouseY = GetMouseY();

				Cmplx mousePos((mouseX - centerX) / scale, (mouseY - centerY) / scale);
				data.push_back(mousePos);
			}

			if (GetKey(olc::SPACE).bPressed && data.size() > 0) {
				timeFactor = 2.0 / data.size();


				if (mode == DrawFunction) {
					mode = DisplaySeries;
					N_terms = 1;
					updateSeries(N_terms);
				}
				else {
					mode = DisplayDFT;
					N_terms = data.size();
					std::vector<Cmplx> dft_out = DFT(data);
					for (int i = 0; i < dft_out.size(); i++) {
						series.push_back(FourierTerm(dft_out[i] / dft_out.size(), i));
					}
					calculatePath(10000);
				}
			}
		}

		if (mode == DisplaySeries || mode == DisplayDFT) {
			if (GetKey(olc::F).bPressed) {
				followTip = !followTip;
			}
			if (GetKey(olc::C).bPressed) {
				drawCircles = !drawCircles;
			}
			if (GetKey(olc::R).bPressed) {
				timeFactor /= 1.5;
			}
			if (GetKey(olc::T).bPressed) {
				timeFactor *= 1.5;
			}
			if (GetKey(olc::D).bPressed) {
				drawPath = !drawPath;
			}
			if (GetKey(olc::X).bPressed) {
				drawFunction = !drawFunction;
			}
			if (GetKey(olc::P).bPressed && mode == DisplaySeries) {
				mode = DisplaySeries;
				N_terms++;
				updateSeries(N_terms);
			}
			if (GetKey(olc::O).bPressed && mode == DisplaySeries) {
				mode = DisplaySeries;
				N_terms = (N_terms <= 0) ? 0 : N_terms - 1;
				updateSeries(N_terms);
			}
		}


		//Camera control
		if (mode != DrawFunction && mode != InputDFTData) {
			if (GetMouse(0).bPressed) {
				mouseX = GetMouseX();
				mouseY = GetMouseY();
			}
			else if (GetMouse(0).bHeld) {
				int newMouseX = GetMouseX();
				int newMouseY = GetMouseY();

				cameraX -= (newMouseX - mouseX);
				cameraY -= (newMouseY - mouseY);

				mouseX = newMouseX;
				mouseY = newMouseY;
			}
		}
		int scroll = GetMouseWheel();
		if (scroll > 0) {
			scale *= 1.1;
			cameraX *= 1.1; //very dumb but I dont want to redo stuff in the more sensible way
			cameraY *= 1.1;
		}
		else if (scroll < 0) {
			scale /= 1.1;
			cameraY /= 1.1;
			cameraX /= 1.1;
		}

		///////////////////////
		//////RENDERING////////
		///////////////////////
		Clear(olc::BLACK);

		if (mode == DrawFunction) {

			for (int i = 0; i + 1 < data.size(); i++) {
				Cmplx data_curr = data.at(i) * scale;
				Cmplx data_next = data.at((i + 1) % data.size()) * scale;
				DrawLine(centerX - cameraX + data_curr.a, centerY - cameraY + data_curr.b, centerX - cameraX + data_next.a, centerY - cameraY + data_next.b, olc::WHITE);
			}
		}
		else if (mode == InputDFTData) {
			for (int i = 0; i < data.size(); i++) {
				Cmplx data_i = data.at(i) * scale;
				FillCircle(centerX - cameraX + data_i.a, centerY - cameraY + data_i.b, 5, olc::WHITE);
			}
		}
		else if (mode == DisplayDFT || mode == DisplaySeries) {
			Cmplx currPos(0, 0);
			Cmplx finalPos(0, 0);

			if (drawFunction) {
				if (mode == DisplayDFT) {
					for (int i = 0; i < data.size(); i++) {
						Cmplx data_i = data.at(i) * scale;
						olc::Pixel color = ((int)(t * data.size()) % data.size() == i) ? olc::DARK_RED : olc::DARK_YELLOW;
						FillCircle(centerX - cameraX + data_i.a, centerY - cameraY + data_i.b, 5, color);
					}
				}
				else {
					for (int i = 0; i < data.size(); i++) {
						Cmplx data_curr = data.at(i) * scale;
						Cmplx data_next = data.at((i + 1) % data.size()) * scale;
						DrawLine(centerX - cameraX + data_curr.a, centerY - cameraY + data_curr.b, centerX - cameraX + data_next.a, centerY - cameraY + data_next.b, olc::DARK_YELLOW);
					}
				}
			}

			if (followTip) {
				//lazy
				for (int i = 0; i < series.size(); i++) {
					finalPos = finalPos + series.at(i).coeff * cmplxExp(2 * 3.141592 * series.at(i).freq * t) * scale;
				}
				cameraX = finalPos.a;
				cameraY = finalPos.b;
			}

			if (drawPath) {
				renderPath(centerX - cameraX, centerY - cameraY, scale, t, 1, olc::DARK_CYAN);
			}

			//drawing the series itself
			for (int i = 0; i < series.size(); i++) {
				Cmplx newPos = currPos + series.at(i).coeff * cmplxExp(2 * 3.141592 * series.at(i).freq * t) * scale;

				//draw arrowhead
				Cmplx arrow_height = (currPos - newPos) / 8.0;
				Cmplx perp = Cmplx(-arrow_height.b, arrow_height.a) / 2;

				Cmplx p1 = newPos;
				Cmplx p2 = newPos + arrow_height + perp;
				Cmplx p3 = newPos + arrow_height - perp;

				DrawTriangle(centerX - cameraX + p1.a, centerY - cameraY + p1.b, centerX - cameraX + p2.a, centerY - cameraY + p2.b, centerX - cameraX + p3.a, centerY - cameraY + p3.b, olc::WHITE);

				Cmplx headBase = newPos + arrow_height;
				DrawLine(centerX - cameraX + currPos.a, centerY - cameraY + currPos.b, centerX - cameraX + headBase.a, centerY - cameraY + headBase.b, olc::WHITE);

				if (drawCircles) {
					DrawCircle(centerX - cameraX + currPos.a, centerY - cameraY + currPos.b, Cmplx::abs(newPos - currPos), olc::WHITE);
				}

				currPos = newPos;
			}

			FillCircle(centerX - cameraX + currPos.a, centerY - cameraY + currPos.b, 2, olc::CYAN);
		}

		return true;
	}
};