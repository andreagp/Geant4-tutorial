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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 40*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Shape 1
  //  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 0*cm, -12*cm);
        
  // Conical section shape
  /*  
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Cons* solidShape1 =    
    new G4Cons("Shape1", 
    shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
    shape1_phimin, shape1_phimax);
  
  */
  //changed shape1 to G4Tubs cylinder
  
  G4double R1 = 0. *cm;
  G4double R2 = 3. *cm;
  G4double hz = 8. *cm;
  G4double stdeg = 0. *deg;
  G4double spanningangle = 360. *deg;
  
  G4Tubs* solidShape1 =
    new G4Tubs("tracker1", R1, R2, hz, stdeg, spanningangle);
   
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //     
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, 0*cm, 10*cm);

  //create detector wall
  G4double det_x = 12*cm;
  G4double det_y = 12*cm;
  G4double det_z = 2*cm;

  G4Box* detwall =
    new G4Box("detwall", //its name
	      0.5*det_x, 0.5*det_y, 0.5*det_z); //its size
  /*
  // Trapezoid shape       
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;      
  G4Trd* solidShape2 =    
    new G4Trd("Shape2",                      //its name
              0.5*shape2_dxa, 0.5*shape2_dxb, 
              0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
  */           
  G4LogicalVolume* logicdetwall =                         
    new G4LogicalVolume(detwall,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicdetwall,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //add collimator blades on detector wall
  // first collimator blade

  G4ThreeVector pos3 = G4ThreeVector(0, 0*cm, 7*cm);
  G4RotationMatrix *rotm1  = new G4RotationMatrix();
  rotm1->rotateZ(14.*deg);
  rotm1->rotateX(40.*deg);

  
  G4double det_xc1 = 0.01*cm;
  G4double det_yc1 = 12*cm;
  G4double det_zc1 = 4*cm;

  G4Box* detcol1 =
    new G4Box("col1", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol1 =                         
    new G4LogicalVolume(detcol1,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(rotm1,                       //no rotation
                    pos3,                    //at position
                    logiccol1,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // second collimator blade
  G4ThreeVector pos4 = G4ThreeVector(6*cm, 0*cm, 7*cm);
  
  G4Box* detcol2 =
    new G4Box("col2", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol2 =                         
    new G4LogicalVolume(detcol2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logiccol2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    // 3rd  collimator blade
  G4ThreeVector pos5 = G4ThreeVector(-6*cm, 0*cm, 7*cm);
  
  G4Box* detcol3 =
    new G4Box("col3", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol3 =                         
    new G4LogicalVolume(detcol3,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos5,                    //at position
                    logiccol3,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    // 4th collimator blade
  G4ThreeVector pos6 = G4ThreeVector(3*cm, 0*cm, 7*cm);
  
  G4Box* detcol4 =
    new G4Box("col4", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol4 =                         
    new G4LogicalVolume(detcol4,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos6,                    //at position
                    logiccol4,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

      // 5th collimator blade
  G4ThreeVector pos7 = G4ThreeVector(-3*cm, 0*cm, 7*cm);
  
  G4Box* detcol5 =
    new G4Box("col5", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol5 =                         
    new G4LogicalVolume(detcol5,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos7,                    //at position
                    logiccol5,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

      // 6th collimator blade
  G4ThreeVector pos8 = G4ThreeVector(1.5*cm, 0*cm, 7*cm);
  
  G4Box* detcol6 =
    new G4Box("col6", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol6 =                         
    new G4LogicalVolume(detcol6,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos8,                    //at position
                    logiccol6,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

      // 7th collimator blade
  G4ThreeVector pos9 = G4ThreeVector(-1.5*cm, 0*cm, 7*cm);
  
  G4Box* detcol7 =
    new G4Box("col7", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol7 =                         
    new G4LogicalVolume(detcol7,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos9,                    //at position
                    logiccol7,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

      // 8th collimator blade
  G4ThreeVector pos10 = G4ThreeVector(4.5*cm, 0*cm, 7*cm);
  
  G4Box* detcol8 =
    new G4Box("col8", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol8 =                         
    new G4LogicalVolume(detcol8,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos10,                    //at position
                    logiccol8,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

      // 9th collimator blade
  G4ThreeVector pos11 = G4ThreeVector(-4.5*cm, 0*cm, 7*cm);
  
  G4Box* detcol9 =
    new G4Box("col9", //its name
	      0.5*det_xc1, 0.5*det_yc1, 0.5*det_zc1); //its size
  
   G4LogicalVolume* logiccol9 =                         
    new G4LogicalVolume(detcol9,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos11,                    //at position
                    logiccol9,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //sphere placement
  G4Sphere* sphere1=
    new G4Sphere("sphere", 2*cm, 4*cm, 0*deg, 180*deg, 0*deg, 180*deg);

  G4LogicalVolume* logicsphere1 =                         
    new G4LogicalVolume(sphere1,         //its solid
                        shape2_mat,          //its material
                        "Sphere");           //its name

    G4ThreeVector spherepos = G4ThreeVector(0*cm, 0*cm, 0*cm);
    
    new G4PVPlacement(0,                       //no rotation
                    spherepos,                    //at position
                    logicsphere1,             //its logical volume
                    "Sphere",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicdetwall;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
