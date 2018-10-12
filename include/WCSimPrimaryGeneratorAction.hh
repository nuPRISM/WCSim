#ifndef WCSimPrimaryGeneratorAction_h
#define WCSimPrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "jhfNtuple.h"
#include <vector>

#include "WCSimEnumerations.hh"

#include <fstream>

#include "WCSimRootOptions.hh"
#include "TFile.h"
#include "TTree.h"
#include "TNRooTrackerVtx.hh"
#include "TClonesArray.h"

class WCSimDetectorConstruction;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class WCSimPrimaryGeneratorMessenger;
class G4Generator;

class WCSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
struct radioactive_source {
    G4String IsotopeName;
    G4String IsotopeLocation;
    G4double IsotopeActivity;
  };

    public:
        WCSimPrimaryGeneratorAction(WCSimDetectorConstruction*);
        ~WCSimPrimaryGeneratorAction();

    public:
        void GeneratePrimaries(G4Event* anEvent);
        void SetupBranchAddresses(NRooTrackerVtx* nrootrackervtx);
        void OpenRootrackerFile(G4String fileName);
        void CopyRootrackerVertex(NRooTrackerVtx* nrootrackervtx);
        bool GetIsRooTrackerFileFinished(){return (fEvNum==fNEntries);}

        // Gun, laser & gps setting calls these functions to fill jhfNtuple and Root tree
        void SetVtx(G4ThreeVector i)     { vtxs[0] = i; nvtxs = 1; };
        void SetBeamEnergy(G4double i, G4int n = 0)   { beamenergies[n] = i; };
        void SetBeamDir(G4ThreeVector i, G4int n = 0) { beamdirs[n] = i; };
        void SetBeamPDG(G4int i, G4int n = 0)         { beampdgs[n] = i; };
        void SetNvtxs(G4int i)     { nvtxs = i; };
        void SetVtxs(G4int i, G4ThreeVector v)     { vtxs[i] = v; };

        // These go with jhfNtuple
        G4int GetVecRecNumber(){return vecRecNumber;}
        //G4int GetMode() {return mode;};
        InteractionType_t GetMode() {return mode;};
        G4int GetNvtxs() {return nvtxs;};
  G4int GetVtxVol(G4int n = 0) {return vtxsvol[n];};
        G4ThreeVector GetVtx(G4int n = 0) {return vtxs[n];}
        G4int GetNpar() {return npar;};
        G4int GetBeamPDG(G4int n = 0) {return beampdgs[n];};
        G4double GetBeamEnergy(G4int n = 0) {return beamenergies[n];};
        G4ThreeVector GetBeamDir(G4int n = 0) {return beamdirs[n];};
        G4int GetTargetPDG(G4int n = 0) {return targetpdgs[n];};
        G4double GetTargetEnergy(G4int n = 0) {return targetenergies[n];};
        G4ThreeVector GetTargetDir(G4int n = 0) {return targetdirs[n];};

        // older ...
        G4double GetNuEnergy() {return nuEnergy;};
        G4double GetEnergy() {return energy;};
        G4double GetXPos() {return xPos;};
        G4double GetYPos() {return yPos;};
        G4double GetZPos() {return zPos;};
        G4double GetXDir() {return xDir;};
        G4double GetYDir() {return yDir;};
        G4double GetZDir() {return zDir;};

  G4bool GetStorePhotons() {return storephotons;}
  void SetStorePhotons(G4double tparam) {storephotons=tparam;}

        G4String GetGeneratorTypeString();

        void SaveOptionsToOutput(WCSimRootOptions * wcopt);
    
  private:
        WCSimDetectorConstruction*      myDetector;
        G4ParticleGun*                  particleGun;
        G4GeneralParticleSource*        MyGPS;  //T. Akiri: GPS to run Laser
        WCSimPrimaryGeneratorMessenger* messenger;

        // Variables set by the messenger
        G4bool   useMulineEvt;
        G4bool   useRootrackerEvt;
        G4bool   useGunEvt;
        G4bool   useLaserEvt;  //T. Akiri: Laser flag
        G4bool   useGPSEvt;
        G4bool useRadioactiveEvt;
        std::vector<struct radioactive_source> radioactive_sources;
        G4double radioactive_time_window;
        std::fstream inputFile;
        G4String vectorFileName;
        G4bool   GenerateVertexInRock;
        G4bool   usePoissonPMT;
        G4double poissonPMTMean;

        // These go with jhfNtuple
        //G4int mode;
        InteractionType_t mode;
        G4int nvtxs;
  G4int vtxsvol[MAX_N_PRIMARIES];
        G4ThreeVector vtxs[MAX_N_PRIMARIES];
        G4int npar;
        G4int beampdgs[MAX_N_PRIMARIES], targetpdgs[MAX_N_PRIMARIES];
        G4ThreeVector beamdirs[MAX_N_PRIMARIES], targetdirs[MAX_N_PRIMARIES];
        G4double beamenergies[MAX_N_PRIMARIES], targetenergies[MAX_N_PRIMARIES];
        G4int vecRecNumber;

        G4double nuEnergy;
        G4double energy;
        G4double xPos, yPos, zPos;
        G4double xDir, yDir, zDir;

        G4int    _counterRock; 
        G4int    _counterCublic;

        // Counters to read Rootracker event file
        int fEvNum;
        int fNEntries;
        TFile* fInputRootrackerFile;

        // Pointers to Rootracker vertex objects
        // Temporary vertex that is saved if desired, according to WCSimIO macro option
        TTree* fRooTrackerTree;
        TTree* fSettingsTree;
        NRooTrackerVtx* fTmpRootrackerVtx;
        float fNuPrismRadius;
        float fNuBeamAng;
        float fNuPlanePos[3];

        G4bool storephotons;

    public:

        inline TFile* GetInputRootrackerFile(){ return fInputRootrackerFile;}

        inline void SetMulineEvtGenerator(G4bool choice) { useMulineEvt = choice; }
        inline G4bool IsUsingMulineEvtGenerator() { return useMulineEvt; }

        inline void SetRootrackerEvtGenerator(G4bool choice) { useRootrackerEvt = choice; }
        inline G4bool IsUsingRootrackerEvtGenerator() { return useRootrackerEvt; }

        inline void SetGunEvtGenerator(G4bool choice) { useGunEvt = choice; }
        inline G4bool IsUsingGunEvtGenerator()  { return useGunEvt; }

        //T. Akiri: Addition of function for the laser flag
        inline void SetLaserEvtGenerator(G4bool choice) { useLaserEvt = choice; }
        inline G4bool IsUsingLaserEvtGenerator()  { return useLaserEvt; }
  
        inline void SetGPSEvtGenerator(G4bool choice) { useGPSEvt = choice; }
        inline G4bool IsUsingGPSEvtGenerator()  { return useGPSEvt; }

        inline void AddRadioactiveSource(G4String IsotopeName, G4String IsotopeLocation, G4double IsotopeActivity){
            struct radioactive_source r;
            r.IsotopeName = IsotopeName;
            r.IsotopeLocation = IsotopeLocation;
            r.IsotopeActivity = IsotopeActivity;
            radioactive_sources.push_back(r);
        }
        inline std::vector<struct radioactive_source> Radioactive_Sources()  { return radioactive_sources; }
  
        inline void SetRadioactiveEvtGenerator(G4bool choice) { useRadioactiveEvt = choice; }
        inline G4bool IsUsingRadioactiveEvtGenerator()  { return useRadioactiveEvt; }
        inline void SetRadioactiveTimeWindow(G4double choice) { radioactive_time_window = choice; }
        inline G4double GetRadioactiveTimeWindow()  { return radioactive_time_window; }
        inline void OpenVectorFile(G4String fileName) 
        {
            if ( inputFile.is_open() ) 
                inputFile.close();

            vectorFileName = fileName;
            inputFile.open(vectorFileName, std::fstream::in);
	    if ( !inputFile.is_open() ) {
	      G4cout << "Vector file " << vectorFileName << " not found" << G4endl;
	      exit(-1);
	    }

        }
        inline G4bool IsGeneratingVertexInRock() { return GenerateVertexInRock; }
        inline void SetGenerateVertexInRock(G4bool choice) { GenerateVertexInRock = choice; }

        inline void SetPoissonPMT(G4bool choice) { usePoissonPMT = choice; }
        inline G4bool IsUsingPoissonPMT(){ return usePoissonPMT; }
  
        inline void SetPoissonPMTMean(G4double val){ poissonPMTMean = val; }
        inline G4double GetPoissonPMTMean(){ return poissonPMTMean; }
};

#endif


