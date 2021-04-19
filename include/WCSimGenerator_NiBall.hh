#ifndef WCSimGenerator_NiBall_hh
#define WCSimGenerator_NiBall_hh 1

#include "G4ReactionProductVector.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "WCSimDetectorConstruction.hh"
#include "TF2.h"
#include "TGraph.h"
#include "TMath.h"


using namespace std;

class WCSimDetectorConstruction;

class WCSimGenerator_NiBall
{
	public:
		WCSimGenerator_NiBall(WCSimDetectorConstruction* myDC);
		~WCSimGenerator_NiBall();
		void Initialize();

		G4ThreeVector GetRandomVertex(G4int tSymNumber);
		int GetNiGammaMode()			{ return fNiGammaMode; }
		int GetNiGammaMultiplicity() 	 	{ return fNiGammaMultiplicity; }
		double * GetNiGammaEnergy()		{ return fNiGammaEnergy; }
		double * GetNiGammaPosition() 		{ return fNiGammaPosition;}

		void SettingNiBall(double fNiBallPosition[3], double rn[6]);

	private:
	
		int fNiGammaMode;
		int fNiGammaMultiplicity;
		double fNiGammaEnergy[4];
		double fNiGammaPosition[3];
		WCSimDetectorConstruction*      myDetector;

		// Constants of gamma production as taken from SK
		static G4double NiBallBR[37];
		static G4int NiBallMulti[37];
		static G4double NiBallEnergy[4][37];
};
#endif
