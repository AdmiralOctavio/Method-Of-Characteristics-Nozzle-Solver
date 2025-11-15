from rocketcea.cea_obj import CEA_Obj

P_combustion = 3.3 * 10 ** 6 #Pascal

obj = CEA_Obj(oxName = "N2O", fuelName = "Ethanol")
s = obj.get_full_cea_output(Pc = 33, MR = 2, eps = 5, pc_units = 'bar')
configuration = obj.get_IvacCstrTc_ChmMwGam(Pc = 33, MR = 2, eps = 5)

L_combustion = 100 #mm
D_combustion = 50 #mm

M_exit = 2.32
g = float(configuration[4])
T_combustion = float(configuration[2])

print(configuration)

R = 8.314
Mw = 0.0265 #kg/mol
Rs = R / Mw
mdot = 0.25 #kg/s

L = 5 #theoretical throat radius in mm

Shorten_Percentage = 0.9