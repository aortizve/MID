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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class for MID
/// \author: Antonio Ortiz, antonio.ortiz@nucleares.unam.mx

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger
    *DetectorConstruction::fMagFieldMessenger = 0;

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(), fNbOfChambers(0), fLogicAbsorber(NULL),
      fLogicChamber_1(NULL), fAbsorberMaterial(NULL), fChamberMaterial(NULL),
      fStepLimit(NULL), fCheckOverlaps(true) {
  fMessenger = new DetectorMessenger(this);

  fNbOfChambers = 2;
  fLogicChamber_1 = new G4LogicalVolume *[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete[] fLogicChamber_1;
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() {
  // Material definition

  G4NistManager *nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Iron defined using NIST Manager
  fAbsorberMaterial = nistManager->FindOrBuildMaterial("G4_Fe");

  // Xenon gas defined using NIST Manager
  fChamberMaterial =
      nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes() {
  G4Material *air = G4Material::GetMaterial("G4_AIR");

  // Sizes of the principal geometrical components (solids)

  // absorber
  G4double absorberLength = 1000.0; // full length of Target
  G4int n_sec_z_1 = 10;
  G4int n_sec_z_2 = 10;
  G4double width = 5.0;
  G4double gap = 0.2;
  G4double reduction = 35.;
  G4double absorberThickness = 70.;
  G4double R_in_abs_o = 210.;

  G4double absorberInner_r = (R_in_abs_o - reduction);
  G4double absorberOuter_r = absorberInner_r + absorberThickness;

  const int nAbs_sec = 10;
  G4double abs_sec[nAbs_sec + 1] = {-500., -400., -300., -200., -100., 0.,
                                    100.,  200.,  300.,  400.,  500.};
  G4double abs_thickness[nAbs_sec] = {33., 40., 51., 63., 70.,
                                      70., 63., 51., 40., 33.};

  // MID
  G4double trackerLength_test = 1.2 * absorberLength;
  G4double worldLength = 1.25 * absorberLength;
  // G4double chamberWidth = 1000.0 * cm; // width of the chambers

  G4double midLength_1 = 1.0 * absorberLength;
  G4double midLength_2 = 1.05 * absorberLength;
  G4double R_ch1_o = 285.0; // radius first layer of chambers
  G4double R_ch2_o = 295.0; // radius second layer of chambers
  G4double R_ch1 = (R_ch1_o - reduction);
  G4double R_ch2 = (R_ch2_o - reduction);

  G4double chamberThickness = 1.;
  G4double sector_width = 105.0 * (R_ch1 / R_ch1_o);
  G4int nStrips = round((sector_width / (width + gap)));
  G4int n_sectors = 16; // it was 16
  G4int n_sectors_z = n_sec_z_1;

  G4double gap_z = 1.5;
  G4double length = midLength_1 / (n_sectors_z); // layer 1, length of the bars
  length -= gap_z;

  G4cout << " --- constructing MuonID layer 1, number of bars per z sector = "
         << nStrips * n_sectors << "    width=" << width
         << "  thickness=" << chamberThickness << "  length=" << length
         << G4endl;
  G4double omega;
  G4double phi_sect = 2.0 * M_PI / (1.0 * n_sectors);
  absorberLength = 0.5 * absorberLength; // Half length of the Target
  trackerLength_test = 0.5 * trackerLength_test;
  midLength_1 = 0.5 * midLength_1;
  midLength_2 = 0.5 * midLength_2;
  // G4double trackerSize = R_ch1-1.; // Radius of outer chamber

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / mm
         << " mm" << G4endl;

  G4Box *worldS = new G4Box("world", // its name
                            (worldLength / 2) * cm, (worldLength / 2) * cm,
                            (worldLength / 2) * cm);       // its size
  G4LogicalVolume *worldLV = new G4LogicalVolume(worldS,   // its solid
                                                 air,      // its material
                                                 "World"); // its name

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  G4VPhysicalVolume *worldPV =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        worldLV,         // its logical volume
                        "World",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

  // Absorber

  G4ThreeVector positionAbsorber = G4ThreeVector(0, 0, 0);

  G4Tubs *absorberS =
      new G4Tubs("absorber", absorberInner_r * cm, absorberOuter_r * cm,
                 absorberLength * cm, 0. * deg, 360. * deg);
  fLogicAbsorber = new G4LogicalVolume(absorberS, air, "Target", 0, 0, 0);
  new G4PVPlacement(0,                // no rotation
                    positionAbsorber, // at (x,y,z)
                    fLogicAbsorber,   // its logical volume
                    "Target",         // its name
                    worldLV,          // its mother volume
                    false,            // no boolean operations
                    0,                // copy number
                    fCheckOverlaps);  // checking overlaps

  G4VisAttributes *boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));

  G4cout << "There are " << nAbs_sec << " sectors in the absorber " << G4endl;
  G4LogicalVolume *fLogicAbsSec[nAbs_sec];
  for (G4int copyNo = 0; copyNo < nAbs_sec; copyNo++) {

    G4Tubs *absSectorS = new G4Tubs(
        "Abssec_solid",
        (absorberInner_r + (absorberThickness - abs_thickness[copyNo]) / 2.) *
            cm,
        (absorberInner_r + (absorberThickness + abs_thickness[copyNo]) / 2.) *
            cm,
        0.5 * 100. * cm, 0. * deg, 360. * deg);

    fLogicAbsSec[copyNo] = new G4LogicalVolume(absSectorS, fAbsorberMaterial,
                                               "Abssec_LV", 0, 0, 0);

    G4VisAttributes *absSVisAtt =
        new G4VisAttributes(G4Colour(.95, .95, .95, 0.3));
    fLogicAbsSec[copyNo]->SetVisAttributes(absSVisAtt);

    new G4PVPlacement(
        0, // no rotation
        G4ThreeVector(0 * cm, 0 * cm,
                      ((abs_sec[copyNo + 1] + abs_sec[copyNo]) / 2.) *
                          cm), // at (x,y,z)
        fLogicAbsSec[copyNo],  // its logical volume
        "Abssec_PV",           // its name
        fLogicAbsorber,        // its mother  volume
        false,                 // no boolean operations
        copyNo,                // copy number
        fCheckOverlaps);       // checking overlaps
  }

  G4cout << "Absorber is " << 2 * absorberLength << " cm of "
         << fAbsorberMaterial->GetName() << G4endl;

  // MID layer 1

  G4ThreeVector positionMid = G4ThreeVector(0, 0, 0);

  G4Tubs *Mid_1 =
      new G4Tubs("mid_1", (R_ch1)*cm, (R_ch1 + chamberThickness) * cm,
                 midLength_1 * cm, 0. * deg, 360. * deg);
  G4LogicalVolume *midLV_1 = new G4LogicalVolume(Mid_1, air, "Mid_1", 0, 0, 0);
  new G4PVPlacement(0,               // no rotation
                    positionMid,     // at (x,y,z)
                    midLV_1,         // its logical volume
                    "Mid_1",         // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps

  // Visualization attributes
  worldLV->SetVisAttributes(boxVisAtt);

  G4VisAttributes *midVisAtt =
      new G4VisAttributes(false, G4Colour(.0, .0, .0, 0.1));
  midLV_1->SetVisAttributes(midVisAtt);

  // Tracker segments
  G4double x_pos, y_pos, z_pos, phi;
  G4int ibar_global = 0;

  G4int ibar_global_t = 0;
  for (G4int isect_z = 0; isect_z < n_sectors_z; ++isect_z) {

    if (isect_z < n_sec_z_1 / 2) {
      z_pos = (0.5 + isect_z) * (gap_z + length);
      // z_pos = (0.5 + isect_z) * (gap_z + length) - (length / 2.0);
    } else {
      // z_pos = -(0.5 + (isect_z - 5)) * (gap_z + length) + (length / 2.0);
      z_pos = -(0.5 + (isect_z - (n_sec_z_1 / 2))) * (gap_z + length);
    }
    // sectors in phi
    for (G4int isect_r = 0; isect_r < n_sectors; ++isect_r) {
      G4int counter = 1;
      phi = isect_r * phi_sect;
      for (G4int ibar = 0; ibar < nStrips; ++ibar) {

        if (ibar <= G4int((sector_width / 2.0) / (width + gap))) {
          omega = std::atan(ibar * (gap + width) / R_ch1);
          x_pos = std::sqrt(ibar * (width + gap) * ibar * (width + gap) +
                            R_ch1 * R_ch1) *
                  std::cos(phi - omega);
          y_pos = std::sqrt(ibar * (width + gap) * ibar * (width + gap) +
                            R_ch1 * R_ch1) *
                  std::sin(phi - omega);
        } else {
          omega = std::atan(counter * (gap + width) / R_ch1);
          x_pos = std::sqrt(counter * (width + gap) * counter * (width + gap) +
                            R_ch1 * R_ch1) *
                  std::cos(phi + omega);
          y_pos = std::sqrt(counter * (width + gap) * counter * (width + gap) +
                            R_ch1 * R_ch1) *
                  std::sin(phi + omega);
          counter++;
        }

        G4Box *scintillator_s =
            new G4Box("aBoxSolid", 0.5 * chamberThickness * cm,
                      0.5 * width * cm, 0.5 * length * cm);

        G4LogicalVolume *scintillator_lv = new G4LogicalVolume(
            scintillator_s, fChamberMaterial, "Chamber_LV_1", 0, 0, 0);

        G4RotationMatrix *rm = new G4RotationMatrix;
        rm->rotateZ(-phi * rad);

        new G4PVPlacement(rm, G4ThreeVector(x_pos * cm, y_pos * cm, z_pos * cm),
                          scintillator_lv,
                          std::string("scintillator_pv_") +
                              std::to_string(ibar_global),
                          midLV_1, false, ibar_global_t, 0);

        G4VisAttributes *visAttributesScint =
            new G4VisAttributes(G4Colour(0.9, 0., 0.9, 0.05));
        scintillator_lv->SetVisAttributes(visAttributesScint);

        ibar_global++;
        ibar_global_t++;
      }
    }
  }

  // MID layer 2
  G4double sector_width_2 = 110.0 * (R_ch2 / R_ch2_o); // it was 110.
  G4int nStrips_ch2 = round((2. * midLength_2 / (width + gap)));
  G4int ibar_global2 = 0;
  G4double chamber2_length = (2. * midLength_2) / (n_sec_z_2);
  chamber2_length -= gap_z;
  G4int n_bars_chamber2 = round(chamber2_length / (width + gap));
  G4cout << " --- constructing MuonID layer 2, number of bars = "
         << nStrips_ch2 * n_sectors << "    width=" << width
         << "  thickness=" << chamberThickness << "  length=" << sector_width_2
         << G4endl;

  G4Tubs *Mid_2 =
      new G4Tubs("mid_2", (R_ch2)*cm, (R_ch2 + chamberThickness) * cm,
                 midLength_2 * cm, 0. * deg, 360. * deg);
  G4LogicalVolume *midLV_2 = new G4LogicalVolume(Mid_2, air, "Mid_2", 0, 0, 0);
  new G4PVPlacement(0,               // no rotation
                    positionMid,     // at (x,y,z)
                    midLV_2,         // its logical volume
                    "Mid_2",         // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps

  G4VisAttributes *mid2VisAtt =
      new G4VisAttributes(false, G4Colour(.0, .0, .0, 0.1));
  midLV_2->SetVisAttributes(mid2VisAtt);

  for (G4int isect_r = 0; isect_r < n_sectors; ++isect_r) {
    G4int ibar_neg_eta = 0;

    G4int counter2 = 0;
    phi = isect_r * phi_sect;
    x_pos = R_ch2 * std::cos(phi * rad);
    y_pos = R_ch2 * std::sin(phi * rad);

    G4RotationMatrix *rm = new G4RotationMatrix;
    rm->rotateZ(-phi * rad);
    int internal_bar_counter = 0;
    int j_sec_z_neg_eta = 0;
    for (G4int j_sec_z = 0; j_sec_z < n_sec_z_2; ++j_sec_z) {

      for (G4int ibar = internal_bar_counter;
           ibar < internal_bar_counter + n_bars_chamber2; ++ibar) {
        if (j_sec_z < int(n_sec_z_2 / 2.)) {
          z_pos = 0.5 * (2 * ibar + 1) * width + 0.5 * gap_z + j_sec_z * gap_z +
                  (ibar - j_sec_z) * gap;

        } else {
          z_pos = 0.5 * (2 * ibar_neg_eta + 1) * width + 0.5 * gap_z +
                  j_sec_z_neg_eta * gap_z +
                  (ibar_neg_eta - j_sec_z_neg_eta) * gap;
          z_pos *= -1.;
          ibar_neg_eta++;
        }
        counter2++;
        G4Box *scintillator_s =
            new G4Box("aBoxSolid", 0.5 * chamberThickness * cm,
                      0.5 * sector_width_2 * cm, 0.5 * width * cm);
        G4LogicalVolume *scintillator_lv = new G4LogicalVolume(
            scintillator_s, fChamberMaterial, "Chamber_LV_2", 0, 0, 0);

        new G4PVPlacement(rm, G4ThreeVector(x_pos * cm, y_pos * cm, z_pos * cm),
                          scintillator_lv,
                          std::string("scintillator2_pv_") +
                              std::to_string(ibar_global2),
                          midLV_2, false, ibar_global_t, 0);

        ibar_global2++;
        ibar_global_t++;

        G4VisAttributes *visAttributesScint =
            new G4VisAttributes(G4Colour(0.9, 0., 0.9, 0.05));
        scintillator_lv->SetVisAttributes(visAttributesScint);
      } // end loop bars
      internal_bar_counter = counter2;
      if (j_sec_z > int(n_sec_z_2 / 2.)) {
        j_sec_z_neg_eta++;
      }
    }
  }

  G4cout << "------- There are " << ibar_global << "  in the first layer"
         << "    width=" << width << "  thickness=" << chamberThickness
         << "  length=" << length << G4endl;

  G4cout << "------- There are " << ibar_global2 << "  in the second layer"
         << "    width=" << width << "  thickness=" << chamberThickness
         << "  length=" << sector_width_2 << G4endl;

  G4VisAttributes *absVisAtt = new G4VisAttributes(G4Colour(0., 0., 0., 0.));
  fLogicAbsorber->SetVisAttributes(absVisAtt);

  // small sensitive disk to check the true particle information before absorber
  G4ThreeVector positionTracker = G4ThreeVector(0, 0, 0);

  G4Tubs *trackerTest =
      new G4Tubs("trackerTest", (absorberInner_r - 0.1 - 10.) * cm,
                 (absorberInner_r - 10.) * cm, trackerLength_test * cm,
                 0. * deg, 360. * deg);
  G4LogicalVolume *trackerTestLV =
      new G4LogicalVolume(trackerTest, air, "Tracker", 0, 0, 0);
  new G4PVPlacement(0,               // no rotation
                    positionTracker, // at (x,y,z)
                    trackerTestLV,   // its logical volume
                    "TrackerTest",   // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,        // copy number
                    fCheckOverlaps); // checking overlaps

  G4Tubs *chamberStest =
      new G4Tubs("Chambertest_solid", (absorberInner_r - 0.1 - 10.) * cm,
                 (absorberInner_r - 10.) * cm, trackerLength_test * cm,
                 0. * deg, 360. * deg);

  G4LogicalVolume *fLogicChamberTest = new G4LogicalVolume(
      chamberStest, fChamberMaterial, "ChamberTest_LV", 0, 0, 0);

  G4VisAttributes *chamberVisAtt =
      new G4VisAttributes(false, G4Colour(.0, .0, 0., 0.01));
  fLogicChamberTest->SetVisAttributes(chamberVisAtt);
  trackerTestLV->SetVisAttributes(chamberVisAtt);

  new G4PVPlacement(0,                      // no rotation
                    G4ThreeVector(0, 0, 0), // at (x,y,z)
                    fLogicChamberTest,      // its logical volume
                    "ChamberTest_PV",       // its name
                    trackerTestLV,          // its mother  volume
                    false,                  // no boolean operations
                    0,                      // copy number
                    fCheckOverlaps);        // checking overlaps

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField() {
  // Sensitive detectors
  // test chamber
  G4String trackerChamberTestSDname = "TrackerChamberTestSD";
  TrackerSD *aTrackerTestSD =
      new TrackerSD(trackerChamberTestSDname, "TrackerTestHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerTestSD);
  // of "Chamber_LV".
  SetSensitiveDetector("ChamberTest_LV", aTrackerTestSD, true);

  // layer 1
  SetSensitiveDetector("Chamber_LV_1", aTrackerTestSD, true);
  // layer 2
  SetSensitiveDetector("Chamber_LV_2", aTrackerTestSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  // FIXME: Current implementation assumes a uniform magnetic field everywhere
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialName) {
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial = nistManager->FindOrBuildMaterial(materialName);

  if (fAbsorberMaterial != pttoMaterial) {
    if (pttoMaterial) {
      fAbsorberMaterial = pttoMaterial;
      if (fLogicAbsorber)
        fLogicAbsorber->SetMaterial(fAbsorberMaterial);
      G4cout << G4endl << "----> The absorber is made of " << materialName
             << G4endl;
    } else {
      G4cout << G4endl
             << "-->  WARNING from SetTargetMaterial : " << materialName
             << " not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetChamberMaterial(G4String materialName) {
  G4NistManager *nistManager = G4NistManager::Instance();

  G4Material *pttoMaterial = nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
    if (pttoMaterial) {
      fChamberMaterial = pttoMaterial;
      for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++) {
        if (fLogicChamber_1[copyNo])
          fLogicChamber_1[copyNo]->SetMaterial(fChamberMaterial);
      }
      G4cout << G4endl << "----> The chambers are made of " << materialName
             << G4endl;
    } else {
      G4cout << G4endl
             << "-->  WARNING from SetChamberMaterial : " << materialName
             << " not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStep(G4double maxStep) {
  if ((fStepLimit) && (maxStep > 0.))
    fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps) {
  fCheckOverlaps = checkOverlaps;
}
