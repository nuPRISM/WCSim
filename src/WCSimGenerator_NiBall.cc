/*********************************************************************************/
/**     WCSimGenerator_NiBall                                                   **/
/**     Author: Pablo                                                           **/
/**     Date: March 28th 2021                                                   **/
/**     Desc: Get the energy, mode, branching ration and multiplicity of gammas **/
/**     produced by Ni+n-->Ni+gamma reaction of the 252Cf+Ni calibration source **/
/*********************************************************************************/

#include "WCSimGenerator_NiBall.hh"
#include "TRandom3.h"

// Setting parameters

double WCSimGenerator_NiBall::NiBallBR[37]   = {}; // Cumulative branching ratio for the Ni ball source
double WCSimGenerator_NiBall::NiBallMulti[37] = {}; // Multiplicities for each mode
double WCSimGenerator_NiBall::NiBallEnergy[4][37] = {}; // Gamma energies for each mode

WCSimGenerator_NiBall::WCSimGenerator_NiBall(WCSimDetectorConstruction* myDC) {
        myDetector = myDC;
        this->Initialize();
}

WCSimGenerator_NiBall::~WCSimGenerator_NiBall() {
}

void WCSimGenerator_Radioactivity::Initialize() {
	NiBallBR   = {0.3761,0.5587,0.5919,0.6010,0.6032,0.6228,0.6259,0.6398,0.6434,0.6460,0.6516,0.6724,0.6736,0.6757,0.6787,0.6821,0.6887,0.6911,0.6932,0.6951,0.6976,0.7060,0.7086,0.7118,0.7136,0.8038,0.8565,0.8581,0.8593,0.8609,0.8713,0.9845,0.9864,0.9911,0.9955,0.9989,1.0000};
	NiBallMulti = {1, 2, 2, 2, 3, 3, 2, 3, 2, 3, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 1, 2, 2, 3, 2, 2, 1, 2, 3, 3, 4, 3};
	NiBallEnergy = {{9.000, 0., 0., 0.,}, {8.534, 0.466, 0.  ,0.}, {8.122, 0.878 , 0., 0.}, {7.698, 1.302 ,0.  ,0.}, {7.266, 0.856, 0.878, 0.}, {6.584, 1.950 ,0.466,0.}, {6.106, 2.894 , 0., 0.}, {6.106, 2.555 ,0.340,0.}, {5.974, 3.026 , 0., 0.}, {5.974, 2.686 ,0.340,0.}, {5.817, 3.182 , 0., 0.}, {5.817, 2.842 ,0.340,0.}, {5.437, 2.686, 0.878,0.}, {5.313, 3.347 ,0.340, 0.}, {5.313, 3.221, 0.466,0.}, {5.313, 3.687 ,0.878, 0.}, {5.270, 3.266, 0.466,0.}, {4.977, 3.555 ,0.466 ,0.}, {4.859, 4.141, 0.   ,0.}, {4.859, 3.675 ,0.466 ,0.}, {4.284, 4.716, 0. , 0.}, {4.031, 4.969 ,0.    ,0.}, {3.930, 4.192, 0.878,0.}, {7.820, 0.    ,0.    ,0.}, {7.537, 0.283, 0.   ,0.}, {6.720, 1.031 ,0.    ,0.}, {6.720, 0.816, 0.284,0.}, {6.635, 1.185 ,0.    ,0.}, {5.696, 2.124, 0.   ,0.}, {6.838, 0.    ,0.    ,0.}, {6.682, 0.156, 0.   ,0.}, {6.320, 0.362 ,0.156 ,0.}, {5.837, 0.846, 0.156,0.}, {5.837,0.483,0.362, 0.156}, {5.515, 1.168 , 0.156 , 0.}};
}

int WCSimGenerator_Radioactivity::GetNiGammaMode() {
	randGen = new TRandom3();
        double rn = randGen->Rndm();
	while (rn > NiBallBR[i]) i++;
	return i;
}

double WCSimGenerator_Radioactivity::GetNiGammaMultiplicity(int mode) {
	return NiBallMulti[mode];
}

double WCSimGenerator_Radioactivity::GetNiGammaEnergy(int mode){
	double eg[4];
	for (int i=0; i<4, i++) eg[i] = NiBallEnergy[i][mode];
	return eg;
}
