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
	std::array<double, 37> ni_branching_ratio = {1.,0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	std::array<int, 37> ni_multiplicity = {1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	std::array<double, 37> ni_energy_1stphoton = {9.000, 0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	std::array<double, 37> ni_energy_2ndphoton = {0., 0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	std::array<double, 37> ni_energy_3rdphoton = {0., 0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	std::array<double, 37> ni_energy_4thphoton = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
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

