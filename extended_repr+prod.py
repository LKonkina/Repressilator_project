import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

#Import Euler integrator for solving ODE system of chemical species inside the cells
from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator

max_cells = 100000


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False)
    integ = CLEulerIntegrator(sim, 9, max_cells)
    # or integ = CLEulerIntegrator(sim, 12, max_cells) if I code synthases

    # use this file for reg too
    regul = ModuleRegulator(sim)
    # Only biophys and regulation
    sim.init(biophys, regul, None, integ)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0)) 

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    sim.pickleSteps = 20



def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.0 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    # cell.growthRate = 0.6
    cell.growthRate = 1.0
    # Specify initial concentration of chemical species
    cell.species[:] = [0,0,0,0,0,0,0,0,0,0,0,0]


def specRateCL():
    return '''
    float laci = species[0];
    float mcherry = species[1];
    float scb2 = species[2];
    float tetr = species[3];
    float mven = species[4];
    float c14 = species[5];
    float ci = species[6];
    float cfp = species[7];
    float c4 = species[8];
    float scbb = species[9];
    float cini = species[10];
    float rhli = species[11];

    const float k_plac = 10.f;  #### to be determined
    const float n_plac = 2.f;   #### to be determined
    const float a0_plac = 10000.f;   #### to be determined
    const float a1_plac = 0.1f;   #### to be determined
    const float mcherry_degr_rate = 1.f; #### to be determined
    const float scbb_degr_rate = 1.f; #### to be determined
    const float scb2_k2 = 1.f; #### to be determined
    const float scb2_degr_rate = 1.f; #### to be determined

    const float k_pscb = 10.f;  #### to be determined
    const float n_pscb = 2.f;   #### to be determined
    const float a0_pscb = 0.1f;   #### to be determined
    const float a1_pscb = 10000.f;   #### to be determined
    const float tetr_degr_rate = 1.f; #### to be determined

    const float k_ptetr = 10.f;  #### to be determined
    const float n_ptetr = 2.f;   #### to be determined
    const float a0_ptetr = 10000.f;   #### to be determined
    const float a1_ptetr = 0.1f;   #### to be determined
    const float mven_degr_rate = 1.f; #### to be determined
    const float cini_degr_rate = 1.f; #### to be determined
    const float c14_k2 = 1.f; #### to be determined
    const float c14_degr_rate = 1.f; #### to be determined

    const float k_pcin = 10.f;  #### to be determined
    const float n_pcin = 2.f;   #### to be determined
    const float a0_pcin = 10000.f;   #### to be determined
    const float a1_pcin = 0.1f;   #### to be determined
    const float ci_degr_rate = 1.f; #### to be determined

    const float k_pci = 10.f;  #### to be determined
    const float n_pci = 2.f;   #### to be determined
    const float a0_pci = 10000.f;   #### to be determined
    const float a1_pci = 0.1f;   #### to be determined
    const float cfp_degr_rate = 1.f; #### to be determined
    const float rhli_degr_rate = 1.f; #### to be determined
    const float c4_k2 = 1.f; #### to be determined
    const float c4_degr_rate = 1.f; #### to be determined
    
    const float k_prhl = 10.f;  #### to be determined
    const float n_prhl = 2.f;   #### to be determined
    const float a0_prhl = 10000.f;   #### to be determined
    const float a1_prhl = 0.1f;   #### to be determined
    const float laci_degr_rate = 1.f; #### to be determined

    float plac_r = (laci/k_plac)**n_plac;
    float plac_expr_rate = (a0_plac+a1_plac*plac_r)/(1+plac_r);
    rates[1] = plac_expr_rate - (mcherry_degr_rate + growthRate)*mcherry;
    rates[9] = plac_expr_rate - (scbb_degr_rate + growthRate)*scbb;
    rates[2] = scbb * scb2_k2 - scb2 * scb2_degr_rate;

    float pscb_r = (scb2/k_pscb)**n_pscb;
    float pscb_expr_rate = (a0_pscb+a1_pscb*pscb_r)/(1+pscb_r);
    rates[3] = pscb_expr_rate - (tetr_degr_rate + growthRate)*tetr;

    float ptetr_r = (tetr/k_ptetr)**n_ptetr;
    float ptetr_expr_rate = (a0_ptetr+a1_ptetr*ptetr_r)/(1+ptetr_r);
    rates[4] = ptetr_expr_rate - (mven_degr_rate + growthRate)*mven;
    rates[10] = ptetr_expr_rate - (cini_degr_rate + growthRate)*cini;
    rates[5] = cini * c14_k2 - c14 * c14_degr_rate;

    float pcin_r = (c14/k_pcin)**n_pcin;
    float pcin_expr_rate = (a0_pcin+a1_pcin*pcin_r)/(1+pcin_r);
    rates[6] = pcin_expr_rate - (ci_degr_rate + growthRate)*ci;

    float pci_r = (ci/k_pci)**n_pci;
    float pci_expr_rate = (a0_pci+a1_pci*pci_r)/(1+pci_r);
    rates[7] = pci_expr_rate - (cfp_degr_rate + growthRate)*cfp;
    rates[11] = pci_expr_rate - (rhli_degr_rate + growthRate)*rhli;
    rates[8] = rhli * c4_k2 - c4 * c4_degr_rate;

    float prhl_r = (c4/k_prhl)**n_prhl;
    float prhl_expr_rate = (a0_prhl+a1_prhl*prhl_r)/(1+prhl_r);
    rates[0] = prhl_expr_rate - (laci_degr_rate + growthRate)*laci;

    '''

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        cell.color = [cell.species[1]/500, cell.species[4]/500, cell.species[7]/500]
        if cell.volume > cell.targetVol:
            a = 1
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.0 + random.uniform(0.0,0.5)
    d2.targetVol = 3.0 + random.uniform(0.0,0.5)
