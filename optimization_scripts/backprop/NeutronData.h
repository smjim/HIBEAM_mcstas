// NeutronData.h

#ifndef NEUTRONDATA_H
#define NEUTRONDATA_H

#include <vector>

struct NeutronData {
    // position [cm], u velocity from sin() cos() of angle, t [ms], v velocity [m/s]
    double index, pdgcode, ekin, x, y, z, vx, vy, vz, t, weight, pol_x, pol_y, pol_z, userflags;
};

std::vector<NeutronData> select_circle(const std::string& inFile, double x0, double y0, double r);

#endif // NEUTRONDATA_H

