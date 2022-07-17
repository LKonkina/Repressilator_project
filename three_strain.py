import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

from CellModeller.Signalling.GridDiffusion import GridDiffusion #add
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add

max_cells = 10000

#Specify parameter for solving diffusion dynamics #Add
grid_dim = (80, 80, 8) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-160, -160, -16) # where to place the diffusion space onto simulation space

n_signals = 3
n_species = 7


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False)
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [10.0, 10.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig)
    
    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)
    
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(-6.0,0,0))  #Add
    sim.addCell(cellType=1, pos=(0.0,0,0)) #Add
    sim.addCell(cellType=2, pos=(6.0,0,0)) #Add
    
    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    sim.addRenderer(sigrend) #Add
    
    sim.pickleSteps = 2

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    # Specify initial concentration of chemical species
    cell.species[:] = [0.0]*n_species
    # Specify initial concentration of signaling molecules
    cell.signals[:] = [0.0]*n_signals

cl_prefix = \
    '''
        const float Dc14 = 1.0f; ### to be determined
        const float Dc4 = 1.0f; ### to be determined
        const float Dscb2 = 1.0f; ### to be determined

        float c4_in = species[0];
        float c4 = signals[0];
        const float c4_degr_rate = 1.f; #### to be determined
        
        float scb2_in = species[1];
        float scb2 = signals[1];
        const float scb2_degr_rate = 1.f; #### to be determined
        
        float c14_in = species[2];
        float c14 = signals[2];
        const float c14_degr_rate = 1.f; #### to be determined

        float laci = species[3];
        const float laci_degr_rate = 1.f; #### to be determined

        float mcherry = species[4];
        const float mcherry_degr_rate = 1.f; #### to be determined

        float mven = species[5];
        float cfp = species[6]
        const float mven_degr_rate = 1.f; #### to be determined
        const float cfp_degr_rate = 1.f; #### to be determined

        const float k_plac = 10.f;  #### to be determined
        const float n_plac = 2.f;   #### to be determined
        const float a0_plac = 10000.f;   #### to be determined
        const float a1_plac = 0.1f;   #### to be determined
        
        const float k_prhl = 10.f;  #### to be determined
        const float n_prhl = 2.f;   #### to be determined
        const float a0_prhl = 10000.f;   #### to be determined
        const float a1_prhl = 0.1f;   #### to be determined

        const float k_pscb = 10.f;  #### to be determined
        const float n_pscb = 2.f;   #### to be determined
        const float a0_pscb = 10000.f;   #### to be determined
        const float a1_pscb = 0.1f;   #### to be determined
        
        const float k_pcin = 10.f;  #### to be determined
        const float n_pcin = 2.f;   #### to be determined
        const float a0_pcin = 10000.f;   #### to be determined
        const float a1_pcin = 0.1f;   #### to be determined        

        '''


# Dc14 = diffusion rate of C14 through the cell membrane
# Dc4 = diffusion rate of C4 through the cell membrane

def specRateCL(): # Add if/else, new species
    global cl_prefix
    return cl_prefix + '''
        if (cellType==0){
        rates[0] = Dc4*(c4-c4_in)*area/gridVolume - c4_degr_rate*c4;
        float prhl_r = (c4/k_prhl)**n_prhl;
        float prhl_expr_rate = (a0_prhl+a1_prhl*prhl_r)/(1+prhl_r);
        rates[3] = prhl_expr_rate - (laci_degr_rate + growthRate)*laci;

        float plac_r = (laci/k_plac)**n_plac;
        float plac_expr_rate = (a0_plac+a1_plac*plac_r)/(1+plac_r);
        rates[4] = plac_expr_rate - (mcherry_degr_rate + growthRate)*mcherry;
        rates[1] = plac_expr_rate - (scb2_degr_rate + growthRate)*scb2 + Dscb2*(scb2-scb2_in)*area/gridVolume;
        rates[2] = Dc14*(c14-c14_in)*area/gridVolume - c14_degr_rate*c14;
        

        } elif (cellType==1){
        rates[1] = Dscb2*(scb2-scb2_in)*area/gridVolume - scb2_degr_rate*scb2;
        float pscb_r = (scb2/k_pscb)**n_pscb;
        float pscb_expr_rate = (a0_pscb+a1_pscb*pscb_r)/(1+pscb_r);
        rates[3] = pscb_expr_rate - (laci_degr_rate + growthRate)*laci;

        float plac_r = (laci/k_plac)**n_plac;
        float plac_expr_rate = (a0_plac+a1_plac*plac_r)/(1+plac_r);
        rates[5] = plac_expr_rate - (mven_degr_rate + growthRate)*mven;
        rates[2] = plac_expr_rate - (c14_degr_rate + growthRate)*c14 + Dc14*(c14-c14_in)*area/gridVolume;
        rates[0] = Dc4*(c4-c4_in)*area/gridVolume - c4_degr_rate*c4;
        

        } else {
        rates[2] = Dc14*(c14-c14_in)*area/gridVolume - c14_degr_rate*c14;
        float pcin_r = (c14/k_pcin)**n_pcin;
        float pcin_expr_rate = (a0_pcin+a1_pcin*pcin_r)/(1+pcin_r);
        rates[3] = pcin_expr_rate - (laci_degr_rate + growthRate)*laci;

        float plac_r = (laci/k_plac)**n_plac;
        float plac_expr_rate = (a0_plac+a1_plac*plac_r)/(1+plac_r);
        rates[0] = plac_expr_rate - (c4_degr_rate + growthRate)*c4 + Dc4*(c4-c4_in)*area/gridVolume;
        rates[6] = plac_expr_rate - (cfp_degr_rate + growthRate)*cfp;
        rates[1] = Dscb2*(scb2-scb2_in)*area/gridVolume - scb2_degr_rate*scb2;
        

        }
        '''


def sigRateCL(): #Add
    global cl_prefix
    return cl_prefix + '''
        rates[0] = -Dc4*(c4-c4_in)*area/gridVolume - c4_degr_rate*c4;
        rates[1] = -Dscb2*(scb2-scb2_in)*area/gridVolume - scb2_degr_rate*scb2;
        rates[2] = -Dc14*(c14-c14_in)*area/gridVolume - c14_degr_rate*c14;
        
        '''


def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        cell.color = [cell.species[4]/500, cell.species[5]/500, cell.species[6]/500]
        if cell.volume > cell.targetVol:
            a = 1
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)
