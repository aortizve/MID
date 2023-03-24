//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "TrackerSD.hh"
#include "Analysis.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String &name, const G4String &hitsCollectionName)
    : G4VSensitiveDetector(name), fHitsCollection(NULL) {
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent *hce) {
  // Create hits collection

  // G4cout << "  Sensitive detector name=" << SensitiveDetectorName << G4endl;
  fHitsCollection =
      new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.)
    return false;

  TrackerHit *newHit = new TrackerHit();

  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(
      aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetMomentum(aStep->GetTrack()->GetMomentum());

  newHit->SetTrackName(aStep->GetTrack()->GetDefinition()->GetParticleName());

  // G4cout<<"
  // vol="<<aStep->GetPreStepPoint()->GetLogicalVolume()->GetName()<<G4endl;

  const G4String currentPhysicalName =
      aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  newHit->SetCurrentVolume(currentPhysicalName);

  fHitsCollection->insert(newHit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent *) {

  // G4LogicalVolumeStore::GetInstance()->GetVolume("");
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  G4int nofHits = fHitsCollection->entries();
  if (nofHits < 1) {
    return;
  }

  G4String pid_prim = "mu-"; // FIXME pid of gun particle

  // G4cout << "------------------  number of hits=" << nofHits << G4endl;

  // COMMENT: since the first layer is the inner tracker-> particle should hit
  // it first
  G4String namepart0 = (*fHitsCollection)[0]->GetTrackName();
  G4double px_hit0 = (*fHitsCollection)[0]->GetMomentum().x() / CLHEP::GeV;
  G4double py_hit0 = (*fHitsCollection)[0]->GetMomentum().y() / CLHEP::GeV;
  G4double pz_hit0 = (*fHitsCollection)[0]->GetMomentum().z() / CLHEP::GeV;
  G4double p_gen0 =
      sqrt(px_hit0 * px_hit0 + py_hit0 * py_hit0 + pz_hit0 * pz_hit0);
  G4double eta_gen0 = 0.5 * log((p_gen0 + pz_hit0) / (p_gen0 - pz_hit0));

  // check whether the particle is really a primary one
  if (namepart0 != pid_prim) {
    return;
  }
  analysisManager->FillH2(0, eta_gen0, p_gen0);
  // G4cout << " pgen=" << p_gen0 << "  etagen=" << eta_gen0
  //        << "   name_part=" << namepart0 << G4endl;
  bool is_in_L1 = false;
  int nbars_1 = 2880; // FIXME

  for (G4int i = 0; i < nofHits; i++) {

    if ((*fHitsCollection)[i]->GetTrackName() !=
        pid_prim) // only primary particle
      continue;
    for (G4int ibar = 0; ibar < nbars_1; ++ibar) {
      G4String name_l1 = std::string("scintillator_pv_") + std::to_string(ibar);

      if ((*fHitsCollection)[i]->GetCurrentVolume() == name_l1) {
        is_in_L1 = true;
        //        G4cout << " is in layer 1 name=" << name_l1 << " eta=" <<
        //        eta_gen0
        //               << " p=" << p_gen0 << G4endl;
      }
    }
  }
  int nbars_2 = 3200; // FIXME
  bool is_in_L2 = false;

  for (G4int i = 0; i < nofHits; i++) {
    if ((*fHitsCollection)[i]->GetTrackName() !=
        pid_prim) // only primary particle
      continue;
    for (G4int ibar = 0; ibar < nbars_2; ++ibar) {
      G4String name_l1 =
          std::string("scintillator2_pv_") + std::to_string(ibar);

      if ((*fHitsCollection)[i]->GetCurrentVolume() == name_l1) {
        is_in_L2 = true;
        //        G4cout << " is in layer 1 name=" << name_l1 << " eta=" <<
        //        eta_gen0
        //               << " p=" << p_gen0 << G4endl;
      }
    }
  }

  if (is_in_L1 && is_in_L2) {
    analysisManager->FillH2(1, eta_gen0, p_gen0);
  }

  if (verboseLevel > 1) {
    G4cout << "   ------ trackID=" << (*fHitsCollection)[0]->GetTrackID()
           << G4endl;

    G4cout << G4endl << "-------->Hits Collection: in this event they are "
           << nofHits << " hits in the tracker chambers: " << G4endl;
    for (G4int i = 0; i < nofHits; i++)
      (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
