# MID
MC simulations for MID
  -----------

 This example shows a simple implementation of the MID geometry in Geant4, an iron absorber (passive volume) is considered, this is followed by two layers of muon chambers. This geometry is intended for very basic studies like: acceptance, impact of magnetic field in the detector efficiency, ...

 1- GEOMETRY DEFINITION

     The geometry is constructed in the DetectorConstruction class.
     The setup consists of a passive absorber (Fe) followed by two concentric
     muon chambers (plastic scintillator). An inner chamber is considered
     to have the kinematic information of the muon before entering in the absorber

     In addition, a global, uniform,  magnetic field can be
     applied using G4GlobalMagFieldMessenger, instantiated in
     DetectorConstruction::ConstructSDandField with a non zero field value,
     or via interactive commands.
     For example:

     /globalField/setValue 0 0 2. tesla
     
2- PHYSICS LIST

    The particle's type and the physic processes which will be available
    in this example are set in the FTFP_BERT physics list. This physics list
    requires data files for electromagnetic and hadronic processes.
