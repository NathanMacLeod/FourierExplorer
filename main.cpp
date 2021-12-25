#define OLC_PGE_APPLICATION
#include <stdio.h>
#include <string.h>
#include "FourierVisualizer.h"
#include "FourierImageTool.h"



int main(int argc, char* argv[])
{
	 
	/*FourierImageTool demo;
	if (demo.Construct(1920, 1080, 1, 1))
		demo.Start();*/

	if (argc > 1) {
		const char* arg = argv[1];
		if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
			printf("Supported arguments are:\n-dft: Run DFT visualization program\n-fs: Run Fourier Series visualization program\n-image [file]: Run image fourier transform program on image [file]\n");
			printf("\nDFT Controls: Input data with the mouse, then press space bar with the mouse to generate and visualize DFT\n");
			printf("\nFS Controls: Input data with the mouse, then press space bar with the mouse to generate and visualize the fourier series. By default there is only one term\n");
			printf("press the 'P' key to increase the number of terms to improve accuracy, or use 'O' to decrease the number of terms.\n");
			printf("\nOther Controls: Press T to speed up time, press R to slow down. D To toggle path visualization, X to toggle drawing the original function.\n");
			printf("Press F to toggle camera following the tip of the series. Zoom with mouse wheel, click and drag to move camera.\n");
			printf("\nImage controls: Click with the mouse on the frequency image to plot a path. Press 'O' to delete all coefficients inside the enclosed area, or 'I' to delete all outside. To Reset, press 'R'.\n");
		}
		else if (strcmp(arg, "-dft") == 0 || strcmp(arg, "-DFT") == 0) {
			FourierVisualizer window(true);
			if (window.Construct(1920, 1080, 1, 1))
				window.Start();
		}
		else if (strcmp(arg, "-fs") == 0 || strcmp(arg, "-FS") == 0) {
			FourierVisualizer window(false);
			if (window.Construct(1920, 1080, 1, 1))
				window.Start();
		}
		else if (strcmp(arg, "-image") == 0 || strcmp(arg, "-IMAGE") == 0) {
			if (argc != 3) {
				printf("Missing arugment for image file. Run with arguments -image [path_to_image]\n");
			}
			else {
				olc::Sprite* image = new olc::Sprite(argv[2]);
				int main_dimension = (image->width > image->height) ? image->width : image->height;
				int size = getNextPow2(main_dimension);
				FourierImageTool window(image);
				if (window.Construct((size + 10) * 3, size + 10, 1, 1))
					window.Start();
			}
		}
		else {
			printf("Unrecognized argument. Run with the -h flag to see options\n");
		}
	}
	else {
		printf("Program arguments required. Run with -h flag to see options\n");
	}
}
