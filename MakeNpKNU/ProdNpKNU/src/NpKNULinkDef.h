#include "NpKNU.hh"

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ namespace npknu;
#pragma link C++ class npknu::Evt+;
#pragma link C++ class npknu::GenInfo+;
#pragma link C++ class npknu::PDFxfx+;
#pragma link C++ class npknu::GenParticle+;
#pragma link C++ class npknu::P4+;
#pragma link C++ class npknu::Cut+;
#pragma link C++ class npknu::Photon+;
#pragma link C++ class npknu::Vertex+;
#pragma link C++ class npknu::Pileup+;
#pragma link C++ class npknu::Electron+;
#pragma link C++ class npknu::Trigger+;
#pragma link C++ class npknu::TriggerObject+;
#pragma link C++ class npknu::MET+;
#pragma link C++ class npknu::Jet+;
#pragma link C++ class npknu::Muon+;
#pragma link C++ class npknu::PUReweight+;
#pragma link C++ class npknu::P4Hist+;
#pragma link C++ class npknu::P4DHist+;
#pragma link C++ class npknu::ElectronHist+;
#pragma link C++ class npknu::MuonHist+;
#pragma link C++ class npknu::METHist+;
#pragma link C++ class npknu::PhotonHist+;
#pragma link C++ class npknu::EtcHist+;

#pragma link C++ function npknu::EffAreaElectron_Spring15_25ns(double);
#pragma link C++ function npknu::EffAreaElectron_Spring15_50ns(double);
#pragma link C++ function npknu::EffAreaPhotonCha(double, bool);
#pragma link C++ function npknu::EffAreaPhotonNeu(double, bool);
#pragma link C++ function npknu::EffAreaPhotonPho(double, bool);

#pragma link C++ function npknu::CompareWild(const char*, const char*);
#pragma link C++ function npknu::GetElectronIdTCA(TClonesArray*, TClonesArray*, int, double, double);
#pragma link C++ function npknu::GetElectronIdTCA(TClonesArray*, int, double, double);
#pragma link C++ function npknu::NumElectronMVATCA(TClonesArray*, unsigned int, double, double, bool);
#pragma link C++ function npknu::GetElectronMVATCA(TClonesArray*, TClonesArray*, unsigned int, double, double, bool);

#pragma link C++ function npknu::GetPhotonIdTCA(TClonesArray*, TClonesArray*, int, double, double, bool);
#pragma link C++ function npknu::GetPhotonIdTCA(TClonesArray*, int, double, double, bool);
#pragma link C++ function npknu::GetMuonCutTCA(TClonesArray*, TClonesArray*, double,double, double);
#pragma link C++ function npknu::GetMuonCutTCA(TClonesArray*, double,double, double);
#pragma link C++ function npknu::GetMuonPreselection(TClonesArray*, double,double, double);
#pragma link C++ function npknu::GetElectronPreselection(TClonesArray*, double,double, double);
#pragma link C++ function npknu::GetMuonVeto2nd(TClonesArray*,  double);
#pragma link C++ function npknu::GetMuonIdTCA(TClonesArray*, TClonesArray*, int, double,double, double);
#pragma link C++ function npknu::GetMuonIdTCA(TClonesArray*, int, double, double, double);
#pragma link C++ function npknu::GetPassTriggerTCA(TClonesArray*, TClonesArray*);
#pragma link C++ function npknu::GetPassTriggerTCA(TClonesArray*);

#pragma link C++ function npknu::GetPassTriggerTCA(TClonesArray*, TClonesArray*, std::string);
#pragma link C++ function npknu::GetPassTriggerTCA(TClonesArray*, std::string);
#pragma link C++ function npknu::GetPassTriggerTCA(TClonesArray*, TClonesArray*, std::vector<std::string>);
#pragma link C++ function npknu::isPassTriggerTCA(TClonesArray*, TClonesArray*, std::vector<std::string>);
#pragma link C++ function npknu::GetPassTriggerTCA(TClonesArray*, std::vector<std::string>);
#pragma link C++ function npknu::GetPassTriggerObjTCA(TClonesArray*, TClonesArray*, std::string);
#pragma link C++ function npknu::isPassTriggerObjTCA(TClonesArray*, TClonesArray*, std::string);
#pragma link C++ function npknu::GetPassTriggerObjTCA(TClonesArray*, std::string);
#pragma link C++ function npknu::GetGoodVertex(TClonesArray*, TClonesArray*);
#pragma link C++ function npknu::GetGoodVertex(TClonesArray*);

#pragma link C++ function npknu::PrintBit(unsigned int);
#pragma link C++ function npknu::GetDirEvents(TString, TString);
#pragma link C++ function npknu::DiLepMass(TClonesArray*);
#pragma link C++ function npknu::GetPileupPU(TClonesArray*);
#pragma link C++ function npknu::GetPileupTrue(TClonesArray*);
#pragma link C++ function npknu::PrintCMD(int, char**);
#pragma link C++ function npknu::PrintElectronTCA(TClonesArray*);

#pragma link C++ function WelcomeNpKNU();


#endif

