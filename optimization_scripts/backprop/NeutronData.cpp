// NeutronData.cpp

#include "NeutronData.h"
#include <iostream>
#include <fstream>
#include <sstream>

std::vector<NeutronData> select_circle(const std::string& inFile, double x0, double y0, double r) {
	std::ifstream inputFileStream(inFile);
	std::vector<NeutronData> selectedData;
	std::string line;
	NeutronData data;

	while (std::getline(inputFileStream, line)) {
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


		// Check if the data point is inside the circle
		double distanceSquared = (data.x - x0) * (data.x - x0) + (data.y - y0) * (data.y - y0);
		if (distanceSquared <= r * r) {
			// Data point is inside the circle, add it to the selected data
			selectedData.push_back(data);
		}
	}

	return selectedData;
}

std::vector<NeutronData> select_rectangle(const std::string& inFile, double x0, double y0, double x1, double y1) {
	std::ifstream inputFileStream(inFile);
	std::vector<NeutronData> selectedData;
	std::string line;
	NeutronData data;

	while (std::getline(inputFileStream, line)) {
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

		// Check if the data point is inside the rectangle 
		if (x0 < data.x && data.x < x1 && y0 < data.y && data.y < y1) {
			// Data point is inside the rectangle, add it to the selected data
			selectedData.push_back(data);
		}
	}

	return selectedData;
}
