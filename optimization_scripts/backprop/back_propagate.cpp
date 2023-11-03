// James Rogers 8-26-2023
// Takes MCPL input file and backpropagates tracks to desired z values, 
// then writes histogram in McStas file format for plotting/ determining min fwhm in different projections
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

struct NeutronData {
	// position [cm], u velocity from sin() cos() of angle, t [ms], v velocity [m/s]
	double index, pdgcode, ekin, x, y, z, vx, vy, vz, t, weight, pol_x, pol_y, pol_z, userflags;
};

int main(int argc, char** argv) {
	if (argc != 4) {
		std::cerr << "Usage: " << argv[0] << " <zvb> <inFile> <outFile>" << std::endl;
		return 1;
	}   
	double zvb = std::stod(argv[1]);
	if (zvb <= 0) {
		std::cerr << "Error: zvb must be positive." << std::endl;
		return 1;
	} 
	std::string inFile(argv[2]);
	std::string outFile(argv[3]);

	const double vb_x = 0.0;			// x position of vb center [m] 
	const double vb_y = 0.0;			// y position of vb center [m] 
	const double target_pos = 65.00;	// z position from where tracks collected [m]

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

		std::ifstream inputFile(inFile);

		if (!inputFile.is_open()) {
			std::cerr << "Error opening input file!" << std::endl;
			return 1;
		}

		std::cout << zpos << std::endl;

		NeutronData data;
		std::string line;

		// Read neutron data from the input file
		while (std::getline(inputFile, line)) {
			std::istringstream iss(line);
			if (!(iss >> data.index >> data.pdgcode >> data.ekin >> data.x >> data.y >> data.z >> data.vx >> data.vy >> data.vz >> data.t >> data.weight >> data.pol_x >> data.pol_y >> data.pol_z >> data.userflags)) {
				std::cerr << "Error reading line: " << line << std::endl;
				continue;
			}

			// Calculate vx, vy, and vz based on ux, uy, uz, and ekin
			// u velocity from sin() cos() of angle, v velocity [m/s]
			double scale = 1.382e10 * std::sqrt(data.ekin); // v [m/s] ~= 1.382e10 * sqrt(E [Mev])
			data.vx = data.vx * scale;
			data.vy = data.vy * scale;
			data.vz = data.vz * scale;
			//std::cout << data.vx << " " << data.vy << " " << data.vz << std::endl;

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

		inputFile.close();

		//std::ofstream outFile("histogram_output_" + std::to_string(zpos) + ".dat");
		//std::ofstream outFile;
		std::ofstream outFile(argv[3]);

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

