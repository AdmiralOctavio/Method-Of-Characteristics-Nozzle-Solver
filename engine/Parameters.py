from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend
from rich import print
import cea
import numpy as np

# --- Default Constants & UI-Linked Parameters ---
P_combustion = 3.4 * 10**6                      
Oxidiser_Fuel_Ratio = 5.13                      
Refinement = 100                                
Thrust = 600                                    
M_exit = 2.23                                   
Shorten_Percentage = 0.75                       
Nozzle_Efficiency = 0.985                       
Combustion_Efficiency = 0.85  
Graph2d = False
Graph3d = True                  

# Aesthetic / Export Options
Graph3d_Fancy = False
Stl = False
Dxf = True
Temperature = False
Write = True
Material = "Copper"
Materials = {
    "Copper": "#db8d5c",
    "Steel": "#525252",
    "Inconel": "#958b87",
    "Titanium": "#B4B1A7",
    "Dodger Blue": "#1E90FF",
}

# Physical dimensions
L_combustion = 83.02                            
Contraction_ratio = 16                          
Chamber_Slope = 45                              
R1 = 10                                         
R2 = 50                                         

# --- Global Containers for Solver Access ---
T_combustion = 0.0
Cstr = 0.0
g = 1.2
Rs = 0.0
ISP_cea = 0.0

def update_engine_data(pc_pa, of_ratio):
    """Calculates chemistry based on UI or manual inputs."""
    global T_combustion, Cstr, g, Rs, ISP_cea
    
    # Fuel Definitions
    ethanol100 = newFuelBlend(["Ethanol", "H2O"], [100, 0])
    fuel_weights = np.array([1, 0, 0.0])
    oxidant_weights = np.array([0.0, 0.0, 1.0])
    T_reactant = np.array([298.15, 298.15, 298.15])
    reac_names = [b"C2H5OH(L)", b"H2O", b"N2O"]
    
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)
    pc_bar = pc_pa / 1e5
    
    # Solve CEA
    solver.solve(solution, weights, pc_bar, 33.55, subar=[16], supar=[5.5], hc=hc, iac=True)

    T_combustion = solution.T[0]
    Cstr = solution.c_star[0]
    ISP_cea = solution.Isp[2] / 9.80665 # Using exit plane Isp

    # Get Gamma and Gas Constant using the rocketcea object for precision
    obj = CEA_Obj(oxName="N2O", fuelName=ethanol100, temperature_units="degK", pressure_units="bar")
    config = obj.get_IvacCstrTc_ChmMwGam(Pc=pc_bar, MR=of_ratio, eps=5)
    
    g = float(config[4])
    Mw = config[3] / 1000  # kg/mol
    Rs = 8.314 / Mw

# Initial run to populate globals
update_engine_data(P_combustion, Oxidiser_Fuel_Ratio)