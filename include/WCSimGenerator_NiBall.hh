#ifndef WCSimGenerator_NiBall_hh
#define WCSimGenerator_NiBall_hh 1

#include "G4ReactionProductVector.hh"
#include "G4ThreeVector.hh"
#include <vector>

#include "WCSimDetectorConstruction.hh"
#include "TF2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"


using namespace std;

class WCSimDetectorConstruction;

class WCSimGenerator_NiBall
{
	public:
		WCSimGenerator_NiBall(WCSimDetectorConstruction* myDC);
		~WCSimGenerator_NiBall();
		void Initialize();

		//void Configuration(G4int iScenario, G4double dLifeTime=0);
		
		G4ThreeVector GetRandomVertex(G4int tSymNumber);
		G4double GetNiGammaMode()		{ return fNiGammaMode; }
		G4double GetNiGammaMultiplicity()	{ return fNiGammaMultiplicity; }
		G4double GetNiGammaEnergy()		{ return fNiGammaEnergy; }

		TRandom3 * randGen;

		
		
	private:
	
/*		// Diffusion functions
		static double DiffusionZ(double x, double origin, double factor);
		static double DiffusionR(double x, double origin, double factor);
		
		// Concentration functions
		static double RadonFormulaR(double* val, double* par);
		static double RadonFormulaZ(double* val, double* par);
		
		// Function
		TF2* thRnFunction;
		G4double fIntegral;
		G4double fIntegralReco;
		
		G4double fConcentrationID;
		G4double fConcentrationReco;
		
		void SetScenario(G4int iScenario);
		static double RadonFormula(double *val, double *par);

		WCSimDetectorConstruction*      myDetector;
		
		// Constant	
		G4int    fScenario;
		static G4double fRnDiffusion_Coef;
		static G4double fRnLambda_Global;
		G4double fRnLambda;
		
		static G4double fRn_PerPMT;
		
		// Parameter
		G4double CONC_BOTTOM;
		G4double CONC_CENTER;
		G4double CONC_INTERMEDIATE;
	
		G4double CONC_MIDDLE;
		G4double CONC_INT_R2_6;
		G4double CONC_INT_R2_7;
		G4double CONC_INT_R2_8;
		G4double CONC_INT_R2_9;
		
		static G4double fRn_Border;
		
		static G4double fHK_Z_min;
		static G4double fHK_Z_max;
		static G4double fHK_R_max;
		static G4double fHK_R2_max;
		
		static G4double fSK_Z_min;
		static G4double fSK_Z_max;
		static G4double fSK_R_max;
		static G4double fSK_R2_max;
		
		static G4double fZ_min_reco;
		static G4double fZ_max_reco;
		static G4double fR_max_reco;
		static G4double fR2_max_reco;
		
		static G4double fZ_min, fZ_max, fR_max , fR2_max;
		static G4double fZ_reco, fR_reco, fR2_reco;
		
		// Parameter array declaration:
		static G4double vParam_Z [RnModel_Bin_Rmax][7];	
		static G4double vParam_R2[RnModel_Bin_Zmax][7];
	
		static G4double vLayer_MinR2_Z[RnModel_Bin_Rmax];
		static G4double vLayer_MaxR2_Z[RnModel_Bin_Rmax];

		static G4double vLayer_MinZ_R2[RnModel_Bin_Zmax];
		static G4double vLayer_MaxZ_R2[RnModel_Bin_Zmax];
*/
		// Constants of gamma production
		static G4double NiBallBR[37];
};
#endif
