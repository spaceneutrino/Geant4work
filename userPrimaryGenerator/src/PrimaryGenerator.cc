#include "PrimaryGenerator.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//{  
PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{

  G4ThreeVector position = G4ThreeVector(0.0,0.0,0.0);
  
  G4double steptime = 76*ns;
  G4double zerotime = G4UniformRand()*steptime;
//  G4double lentime = 10*mks;
  
  for (G4int icollnum = 0; icollnum < 2; icollnum++){
  G4double ptime = steptime*icollnum;
  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, ptime);
  
  G4ParticleDefinition* particleDefinition
           = G4ParticleTable::GetParticleTable()->FindParticle("proton");
           
//  G4int iptclnum = RandPoisson 0.3
  G4int iptclnum = 5;
           
  for (G4int iptcl = 0; iptcl < iptclnum; iptcl++){
  G4PrimaryParticle* p = new G4PrimaryParticle(particleDefinition);
 
  
  G4double tetta = G4UniformRand()*2*3.14159;
  G4double phi = G4UniformRand()*2*3.14159;
 // G4double pt = G4UniformRand() + 0.1;
 // G4double px = pt*cos(phi);
  G4double px = std::sin(tetta)*std::cos(phi);
  G4double py = std::sin(tetta)*std::sin(phi);
  G4double pz = std::cos(tetta);
  G4ThreeVector direction = G4ThreeVector(px,py,pz);
  
  p->SetMomentumDirection(direction);
  p->SetKineticEnergy(1*GeV);
  
  vertex->SetPrimary(p);
 }
 event->AddPrimaryVertex(vertex);  
 }
  
  }
