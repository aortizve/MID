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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
          G4cout << "Using " << analysisManager->GetType() << G4endl;
  G4RunManager::GetRunManager()->SetPrintProgress(1000);     
          analysisManager->SetVerboseLevel(1);
          analysisManager->SetFirstHistoId(0); // default is 0
analysisManager->CreateH1("hHits","number of hits", 1000,0,1000);


analysisManager->CreateH2("hPVsEtaGen","before", 40,-2.,2.,60,0,6);
analysisManager->CreateH2("hPVsEtaRec","after", 40,-2.,2.,60,0,6);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{

// Get analysis manager
          G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

          // Open an output file
          //
          G4String fileName = "mid.root";
          analysisManager->OpenFile(fileName);
 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{

 // print histogram statistics
          //
          G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
          if ( analysisManager->GetH1(0) ) {
G4cout << G4endl << " mean hits = "<<analysisManager->GetH1(0)->mean()
               << G4endl;

}

          // save histograms & ntuple
          //
          analysisManager->Write();
          analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
