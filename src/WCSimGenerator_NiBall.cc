/*********************************************************************************/
/**     WCSimGenerator_NiBall                                                   **/
/**     Author: Pablo                                                           **/
/**     Date: March 28th 2021                                                   **/
/**     Desc: Get the energy, mode, branching ration and multiplicity of gammas **/
/**     produced by Ni+n-->Ni+gamma reaction of the 252Cf+Ni calibration source **/
/*********************************************************************************/

#include "WCSimGenerator_NiBall.hh"
#include <array>

// Setting parameters

double WCSimGenerator_NiBall::NiBallBR[37]   = {}; // Cumulative branching ratio for the Ni ball source
int    WCSimGenerator_NiBall::NiBallMulti[37] = {}; // Multiplicities for each mode
double WCSimGenerator_NiBall::NiBallEnergy[4][37] = {}; // Gamma energies for each mode

WCSimGenerator_NiBall::WCSimGenerator_NiBall(WCSimDetectorConstruction* myDC) {
        myDetector = myDC;
        this->Initialize();
}

WCSimGenerator_NiBall::~WCSimGenerator_NiBall() {
}

void WCSimGenerator_NiBall::Initialize() {
	std::array<double, 37> ni_branching_ratio = {0.3761,0.5587,0.5919,0.6010,0.6032,0.6228,0.6259,0.6398,0.6434,0.6460,0.6516,0.6724,0.6736,0.6757,0.6787,0.6821,0.6887,0.6911,0.6932,0.6951,0.6976,0.7060,0.7086,0.7118,0.7136,0.8038,0.8565,0.8581,0.8593,0.8609,0.8713,0.9845,0.9864,0.9911,0.9955,0.9989,1.0000};
	std::array<int, 37> ni_multiplicity = {1, 2, 2, 2, 3, 3, 2, 3, 2, 3, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 1, 2, 2, 3, 2, 2, 1, 2, 3, 3, 4, 3};
	std::array<double, 37> ni_energy_1stphoton = {9.000, 8.534, 8.122, 7.698, 7.266, 6.584, 6.106, 6.106, 5.974, 5.974, 5.817, 5.817, 5.621, 5.437, 5.437, 5.313, 5.313, 5.313, 5.270, 4.977, 4.859, 4.859, 4.284, 4.031, 3.930, 7.820, 7.537, 6.720, 6.720, 6.635, 5.696, 6.838, 6.682, 6.320, 5.837, 5.837, 5.515};
	std::array<double, 37> ni_energy_2ndphoton = {0., 0.466, 0.878, 1.302, 0.856, 1.950, 2.894, 2.555, 3.026, 2.686, 3.182, 2.842, 3.039, 3.563, 2.686, 3.347, 3.221, 3.687, 3.266, 3.555, 4.141, 3.675, 4.716, 4.969, 4.192, 0., 0.283, 1.031, 0.816, 1.185, 2.124, 0., 0.156, 0.362, 0.846, 0.483, 1.168};
	std::array<double, 37> ni_energy_3rdphoton = {0., 0., 0., 0., 0.878, 0.466, 0., 0.340, 0., 0.340, 0., 0.340, 0., 0., 0.878, 0.340, 0.466, 0.878, 0.466, 0.466, 0., 0.466, 0.,0., 0.878,0., 0., 0., 0.284,0., 0.,0. , 0.,0.156 , 0.156,0.362, 0.156};
	std::array<double, 37> ni_energy_4thphoton = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.156,0.};
	for (int i=0; i<37; i++) {
		NiBallBR[i] = ni_branching_ratio[i];
		NiBallMulti[i] = ni_multiplicity[i];
		NiBallEnergy[0][i] = ni_energy_1stphoton[i];
		NiBallEnergy[1][i] = ni_energy_2ndphoton[i];
		NiBallEnergy[2][i] = ni_energy_3rdphoton[i];
		NiBallEnergy[3][i] = ni_energy_4thphoton[i];
  	}
}

void WCSimGenerator_NiBall::SettingNiBall(double fNiBallPosition[3], double rn[4]) {
	int i = 0;
	while (rn[0] > NiBallBR[i]) i++;
	fNiGammaMode = i;
	fNiGammaMultiplicity = NiBallMulti[fNiGammaMode];
	for (int j=0; j<4; j++) fNiGammaEnergy[j] = NiBallEnergy[j][fNiGammaMode];
        double radius, theta, z; // Match the geometry of Ni ball (taken from SK)
        radius = rn[1] * 2.5;//in cm
        theta  = rn[2] * 2.0*3.14159265;
        z      = rn[3] * 12.0 + 5.0;//in cm
        fNiGammaPosition[0] = fNiBallPosition[0] + radius*cos(theta);
        fNiGammaPosition[1] = fNiBallPosition[1] + radius*sin(theta);
        fNiGammaPosition[2] = fNiBallPosition[2] + z - 8.5;

}

