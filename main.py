import Output
import IsentropicTools as IT
import Parameters as P
import numpy as np
import solver_iteration as SI
from scipy.optimize import fsolve


def run():

    def engine_residuals(variables):
        mdot_guess, mach_guess = variables
        A_t = (mdot_guess * P.Cstr) / P.P_combustion
        r_t = np.sqrt(A_t / np.pi) * 1000 
        Cf, thrust_new, mdot_out, P_exit = SI.main(mdot_guess, r_t, mach_guess)

        thrust_error = thrust_new - P.Thrust
        pressure_error = P_exit - 101325

        return [thrust_error, pressure_error]

    CF_in = IT.estimate_CF(5.5, P.g, P.P_combustion)    # Thrust Coefficient
    mdot_initial = (P.P_combustion * (P.Thrust / (CF_in * P.P_combustion))) / P.Cstr

    initial_guess = [mdot_initial, P.M_exit]

    solution = fsolve(engine_residuals, x0=initial_guess)

    mdot_final, mach_final = solution

    A_t_final = (mdot_final * P.Cstr) / P.P_combustion
    radius_throat_final = np.sqrt(A_t_final / np.pi) * 1000  

    # Initial parameters:
   
    '''Cstr = P.Cstr                                       # C*
    A_t = P.Thrust / (CF_in * P.P_combustion)           # Throat Area
    mdot = P.P_combustion * A_t / Cstr                  # Mass flow-rate
    r_t = np.sqrt(A_t / np.pi) * 1000                   # Throat Radius
    mach_initial = P.M_exit

    Cf, thrust_old, mdot, P_exit = SI.main(mdot, r_t, mach_initial)   # New parameters for above inputs

    A_t = mdot * Cstr / (P.P_combustion)                # New throat area
    r_t = np.sqrt(A_t / np.pi) * 1000                   # New throat radius

    difference = 1

    while difference > 0.001:
        output = SI.main(mdot, r_t, mach_initial)
        thrust_new = output[1]
        difference = abs((thrust_new - thrust_old) / thrust_old)
        mdot = output[2]
        A_t = mdot * Cstr / (P.P_combustion)
        r_t = np.sqrt( A_t/ np.pi) * 1000
        thrust_old = thrust_new

    radius_throat_final = np.sqrt( (output[1] / (output[0] * P.P_combustion)) / (np.pi) ) * 1000
    mdot_final = output[2]'''

    

    Output.outputTable(radius_throat_final, mdot_final, mach_final)


if __name__ == "__main__":
    run()
