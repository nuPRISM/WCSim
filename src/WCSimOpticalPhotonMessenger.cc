#include "WCSimOpticalPhotonMessenger.hh"
#include "WCSimOpticalPhotonTrackInfo.hh"


WCSimOpticalPhotonMessenger::WCSimOpticalPhotonMessenger(WCSimOpticalPhotonTrackInfo* info)
:opInfo(info)
{
  opInfoDir = new G4UIdirectory("/opInfo/");
  
  opEnabled = new G4UIcmdWithABool("/opInfo/Enabled", this);
  opEnabled->SetGuidance("Enables Optical Photon track info tree");
  opEnabled->SetParameterName("opEnabled", true);
  opEnabled->SetDefaultValue(false);
  opEnabled->AvailableForStates(G4State_PreInit, G4State_Idle);

  opKillScatterRef = new G4UIcmdWithABool("/opInfo/KillScatterRef", this);
  opKillScatterRef->SetGuidance("Kill any Optical Photon that scatters or reflects");
  opKillScatterRef->SetParameterName("opKillScatterRef", true);
  opKillScatterRef->SetDefaultValue(false);
  opKillScatterRef->AvailableForStates(G4State_PreInit, G4State_Idle);
}

WCSimOpticalPhotonMessenger::~WCSimOpticalPhotonMessenger(){
  delete opInfoDir;
  delete opEnabled;
  delete opKillScatterRef;
}

void WCSimOpticalPhotonMessenger::SetNewValue(G4UIcommand* command,G4String newValue){

  if(command == opEnabled){
    bool value = opEnabled->GetNewBoolValue(newValue);
    std::cout << "WCSimOpticalPhotonInfo enable Value is " << value << ". New value is " << newValue << std::endl; 
    if(value){
      opInfo->enable();
    } else {
      opInfo->disable();
    }
  } 

  if(command == opKillScatterRef){
    bool value = opKillScatterRef->GetNewBoolValue(newValue);
    std::cout << "Kill Optical Photon on Scatter/Reflect Value is " << value << ". New value is " << newValue << std::endl; 
    opInfo->setKillScatRef( value );
  }

}
