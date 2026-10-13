from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend
from rich import print
import cea
import numpy as np
from engine.Fuels import build_propellant

# Default Constants & UI-Linked Parameters
P_combustion = 3.4 * 10**6
Oxidiser_Fuel_Ratio = 5.13
Refinement = 100
Thrust = 600
M_exit = 2.23
Shorten_Percentage = 0.75
Nozzle_Efficiency = 0.985
Combustion_Efficiency = 0.875
Graph2d = False
Graph3d = True
Ambient_P = 101325
Selected_Oxidiser = "n2o"  # Default
Selected_Fuel = "ethanol"  # Default
Fuel_Concentration = 100.0    # wt% fuel in a fuel-water blend; ignored for unblended fuels

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
L_combustion = 83.02          # mm; used directly only when Lstar_Sizing is False.
                              # When Lstar_Sizing is True this is overwritten with
                              # the length computed from L* (for output/printouts).
Contraction_ratio = 16
Chamber_Slope = 45
R1 = 10
R2 = 50

# Chamber length sizing
Lstar_Sizing = True           # True: chamber length derived from L*; False: manual L_combustion
Lstar_Value = 1.0             # metres; chosen on the UI slider between the propellant's
                              # short (min) and elongated (max) L*. Used directly when set.
Lstar_Option = "nominal"      # legacy fallback: "short" | "nominal" | "elongated"

# Global Containers for Solver Access
T_combustion = 0.0
Cstr = 0.0
g = 1.2
Rs = 0.0
ISP_cea = 0.0
g_exit = 1.2


def get_propellant():
    """PropellantConfig for the selected oxidiser / fuel / fuel concentration."""
    return build_propellant(Selected_Oxidiser, Selected_Fuel, float(Fuel_Concentration))


def of_sweep(pc_pa, n=41):
    """Ideal CEA Isp (s) across an O/F sweep for the selected propellant.

    The sweep spans 0.3-1.4x the stoichiometric O/F, which brackets the Isp
    peak for every oxidiser / fuel pair. Returns (of_values, isp, of_stoich).
    """
    prop = get_propellant()
    of_stoich = prop.make_cea_obj().getMRforER(ERphi=1.0)
    of_values = np.linspace(0.3 * of_stoich, 1.4 * of_stoich, n)

    reac, prod = prop.make_mixture()
    solver = cea.RocketSolver(prod, reactants=reac)
    pc_bar = pc_pa / 1e5
    pressure_ratio = pc_pa / max(Ambient_P, 1.0)

    isp = np.full(n, np.nan)
    for i, of_ratio in enumerate(of_values):
        solution = cea.RocketSolution(solver)
        solver.solve(
            solution, prop.get_weights(of_ratio), pc_bar, pressure_ratio,
            subar=[16], supar=[5.5], hc=prop.get_enthalpy_coeff(of_ratio), iac=True,
        )
        isp[i] = solution.Isp[2] / 9.80665
    return of_values, isp, of_stoich


def update_engine_data(pc_pa, of_ratio):
    """Calculates chemistry based on UI or manual inputs."""
    global T_combustion, Cstr, g, Rs, ISP_cea, g_exit

    # Fuel Definitions
    prop = get_propellant()
    weights = prop.get_weights(of_ratio)
    reac, prod = prop.make_mixture()
    hc = prop.get_enthalpy_coeff(of_ratio)

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)
    pc_bar = pc_pa / 1e5

    # Solve CEA
    solver.solve(
        solution, weights, pc_bar, pc_pa/Ambient_P, subar=[16], supar=[5.5], hc=hc, iac=True
    )

    T_combustion = solution.T[0]
    Cstr = solution.c_star[0]
    ISP_cea = solution.Isp[2] / 9.80665  # Using exit plane Isp

    # Get Gamma and Gas Constant using the rocketcea object for precision
    obj = prop.make_cea_obj()

    config2 = obj.get_exit_MolWt_gamma(Pc=pc_bar, MR=of_ratio, eps=5.77)
    g_exit = config2[1]
    print(f"Exit gamma: {g_exit}")

    # For CFD Stuff

    config3 = obj.get_Throat_MolWt_gamma(Pc=pc_bar, MR=of_ratio, eps=5.77)
    g = config3[1]
    Mw = config3[0] / 1000
    print(f"{g} = THROPAT GAMMA MM TYUMMY")
    Rs = 8.314 / Mw

    # More CFD stuff

    config_chamber = obj.get_Chamber_MolWt_gamma(Pc=pc_bar, MR=of_ratio, eps=5.77)
    MolWt_chamber = config_chamber[0]
    gamma_chamber = config_chamber[1]
    Cp_ch, mu_ch, k_ch, Pr_ch = obj.get_Chamber_Transport(Pc=pc_bar, MR=of_ratio, eps=5.77)
    Cp_ch_SI = Cp_ch * 4184      # cal/(g·K) -> J/(kg·K), same conversion as throat
    mu_ch_SI = mu_ch * 1e-4      # millipoise -> Pa·s
    print(f"{MolWt_chamber} = Chamber Molecular Weight")
    print(f"{gamma_chamber} = Chamber Gamma")
    print(f"{Cp_ch_SI} = Cp in SI chamber")
    print(f"{mu_ch_SI} = Chamber mu")
    




# Initial run to populate globals
update_engine_data(P_combustion, Oxidiser_Fuel_Ratio)
