// James Rogers 8-26-2023
// Takes MCPL input file and backpropagates selected tracks to desired z values, 
// then writes histogram in McStas file format for plotting/ determining distribution characteristics in different projections
// Backpropagates from target_pos to zvb
// Selection is in the form (--rectangle x0 y0 x1 y1), and it will backpropagate only neutrons within the defined x, y bounds.
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include "NeutronData.cpp"
#include "NeutronData.h"

int main(int argc, char** argv) {
	/* Parse arguments */
	if (argc < 6) {
		std::cerr << "Usage: " << argv[0] << " <zdet> <zvb> <inFile> <outFile> (--rectangle x0 y0 x1 y1 OR --circle x0 y0 r)" << std::endl;
		return 1;
	}

	double target_pos = std::stod(argv[1]);	// z position from where tracks collected [m]
	double zvb = std::stod(argv[2]);		// z position where tracks are proagated to [m]

	std::string inFile(argv[3]);
	std::string outFile(argv[4]);
	std::string specifier(argv[5]);

	std::vector<NeutronData> filtered_data;

	if (specifier == "--rectangle" && argc == 10) {
		double x0 = std::stod(argv[6]);
		double y0 = std::stod(argv[7]);
		double x1 = std::stod(argv[8]);
		double y1 = std::stod(argv[9]);

		filtered_data = select_rectangle(inFile, x0, y0, x1, y1);
	}
	else if (specifier == "--circle" && argc == 9) {
		double x0 = std::stod(argv[6]);
		double y0 = std::stod(argv[7]);
		double r = std::stod(argv[8]);

		filtered_data = select_circle(inFile, x0, y0, r);
	}
	else {
		std::cerr << "Invalid specifier or incorrect number of arguments." << std::endl;
		return 1;
	}

	/* Define constants */
	const double vb_x = 0.0;			// x position of vb center [m] 
	const double vb_y = 0.0;			// y position of vb center [m] 

	const double binSize = 0.01;		// bin Size [m]
	const double zvals[] = {zvb};		// z values for which histogram is recorded [m]
		
	std::cout << zvals[0] << std::endl;

	// bounds on histogram plane in x-y [m]
	double xwidth = 1.00;
	double ywidth = 1.00;

	double xmax =  vb_x + xwidth/2.;
	double xmin =  vb_x - xwidth/2.; 
	double ymax =  vb_y + ywidth/2.;
	double ymin =  vb_y - ywidth/2.;

	int numBinsX = static_cast<int>(std::ceil((xmax - xmin) / binSize));
	int numBinsY = static_cast<int>(std::ceil((ymax - ymin) / binSize));

	// this code loops though every z plane value to calculate backpropagation (inefficient)
	for (double zpos : zvals) {
		std::vector<std::vector<double>> I_histogram(numBinsY, std::vector<double>(numBinsX, 0));
		std::vector<std::vector<int>> N_histogram(numBinsY, std::vector<int>(numBinsX, 0));

		// for each track in filtered_data:
		for (const auto& data : filtered_data) {
			// Backpropagate neutrons based on initial position and velocity
			// position [m]
			double dt = (zpos - target_pos) / data.vz;
			//std::cout << dt << std::endl;
			double xAtZ = 0.01*data.x + data.vx * dt;
			double yAtZ = 0.01*data.y + data.vy * dt;
			//std::cout << xAtZ << " " << yAtZ << std::endl;

			if (xAtZ >= xmin && yAtZ >= ymin && xAtZ <= xmax && yAtZ <= ymax) {
				int xBin = static_cast<int>((xAtZ - xmin) / binSize);
				int yBin = static_cast<int>((yAtZ - ymin) / binSize);
				I_histogram[yBin][xBin] += data.weight;
				N_histogram[yBin][xBin] ++;
			}
		}

		std::ofstream outFile(argv[4]);

		if (!outFile.is_open()) {
			std::cerr << "Error opening output file!" << std::endl;
			return 1;
		}

		outFile << "# Format: McCode with text headers" << std::endl;
		outFile << "# Instrument: HIBEAM" << std::endl;
		outFile << "# type: array_2d(" << numBinsX << ", " << numBinsY << ")" << std::endl;
		outFile << "# component: BackTracing Histogram" << std::endl;
		outFile << "# position: " << vb_x << " " << vb_y << " " << zpos << std::endl;
		outFile << "# title: PSD monitor" << std::endl;
		outFile << "# xvar: X" << std::endl;
		outFile << "# yvar: Y" << std::endl;
		outFile << "# xlabel: X position [cm]" << std::endl;
		outFile << "# ylabel: Y position [cm]" << std::endl;
		outFile << "# zvar: I" << std::endl;
		outFile << "# xylimits: " << 100*xmin << " " << 100*xmax << " " << 100*ymin << " " << 100*ymax << std::endl;

		// output I array
		outFile << "# Data [] I:" << std::endl;
		for (int yBin = 0; yBin < numBinsY; ++yBin) {
			for (int xBin = 0; xBin < numBinsX; ++xBin) {
				outFile << I_histogram[yBin][xBin] << " ";
			}
			outFile << std::endl;
		}

		// output dummy Ierr array
		outFile << "# Errors [] I_err:" << std::endl;
		for (int yBin = 0; yBin < numBinsY; ++yBin) {
			for (int xBin = 0; xBin < numBinsX; ++xBin) {
				outFile << 0 << " ";
			}
			outFile << std::endl;
		}

		// output N array
		outFile << "# Events [] N:" << std::endl;
		for (int yBin = 0; yBin < numBinsY; ++yBin) {
			for (int xBin = 0; xBin < numBinsX; ++xBin) {
				outFile << N_histogram[yBin][xBin] << " ";
			}
			outFile << std::endl;
		}

		outFile.close();
	}

	return 0;
}

