from dataclasses import dataclass, field
from functools import lru_cache
from typing import Optional
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
from rocketcea.blends import newFuelBlend
import cea


"""
Fuel definitions
Feel free to add your own oxidisers and fuels depending on what you need to simulate
Parameters.py utilises the new and old cea libraries, so both inputs must be given

"""


@dataclass
class PropellantConfig:
    # For new CEA
    name: str
    reac_names: list
    fuel_weights: np.ndarray
    oxidant_weights: np.ndarray
    T_reactant: np.ndarray

    # For old rocketcea
    ox_name: str = ""
    fuel_name: str = ""

    # Characteristic chamber length L* in metres: (short, nominal, elongated).
    # Drives combustion-chamber volume sizing in CombustionChamber.py.
    Lstar_range: tuple = (1.0, 1.14, 1.27)

    # J/kg added to each reactant's thermo.lib enthalpy (e.g. a latent heat
    # where the library only has the gas phase). None: no correction.
    enthalpy_offsets: Optional[np.ndarray] = None

    def get_lstar(self, option: str = "nominal") -> float:
        """Return characteristic length L* (metres) for the chosen sizing option."""
        idx = {"short": 0, "nominal": 1, "elongated": 2}.get(str(option).lower(), 1)
        return self.Lstar_range[idx]

    def make_mixture(self):
        # Return (reac, prod, weights_fn) for new cea
        reac = cea.Mixture(self.reac_names)
        prod = cea.Mixture(self.reac_names, products_from_reactants=True)
        return reac, prod

    def get_weights(self, of_ratio: float) -> np.ndarray:
        # Convert O/F ratio to mass fractions using this propellant's weight vectors.
        reac, _ = self.make_mixture()
        return reac.of_ratio_to_weights(
            self.oxidant_weights, self.fuel_weights, of_ratio
        )

    def get_enthalpy_coeff(self, of_ratio: float) -> float:
        """Return hc (enthalpy / R) for this propellant at the given O/F."""
        reac, _ = self.make_mixture()
        weights = self.get_weights(of_ratio)
        h = reac.calc_property(cea.ENTHALPY, weights, self.T_reactant)
        if self.enthalpy_offsets is not None:
            # weights are relative, not normalised mass fractions
            h += float(np.dot(weights, self.enthalpy_offsets) / np.sum(weights))
        return h / cea.R

    def make_cea_obj(self) -> CEA_Obj:
        """Return a rocketcea CEA_Obj configured for this propellant."""
        return CEA_Obj(
            oxName=self.ox_name,
            fuelName=self.fuel_name,
            temperature_units="degK",
            pressure_units="bar",
        )


@dataclass(frozen=True)
class Oxidiser:
    name: str
    reac_names: tuple       # new CEA species making up the oxidiser
    weights: tuple          # mass fractions of those species
    T_reactant: tuple
    ox_name: str            # rocketcea name
    # Default L* band for this oxidiser; overridden per fuel in _LSTAR_OVERRIDES.
    Lstar_range: tuple = (1.0, 1.14, 1.27)


@dataclass(frozen=True)
class Fuel:
    name: str
    reac_name: bytes        # new CEA species
    T_reactant: float
    fuel_name: str          # rocketcea name
    # True: fuel can be diluted with water, concentration set on the UI slider.
    water_blend: bool = False
    # J/kg added to the new-CEA species enthalpy (see IPA below).
    enthalpy_offset: float = 0.0


OXIDISERS: dict[str, Oxidiser] = {
    "n2o": Oxidiser(
        name="N2O",
        reac_names=(b"N2O",),
        weights=(1.0,),
        T_reactant=(298.15,),
        ox_name="N2O",
        # ESTIMATE: N2O not tabulated in Sutton; N2O vaporisation tends to
        # need a longer chamber. Tune against test data if available.
        Lstar_range=(0.5, 1.0, 1.8),
    ),
    "lox": Oxidiser(
        name="LOX",
        reac_names=(b"O2(L)",),
        weights=(1.0,),
        T_reactant=(90.15,),
        ox_name="LOX",
        # Sutton & Biblarz: O2/LH2 L* = 0.76-1.02 m. Also used for O2/alcohol
        # (ESTIMATE: not directly tabulated; Sutton's O2/ammonia band is the same).
        Lstar_range=(0.76, 0.89, 1.02),
    ),
    "peroxide98": Oxidiser(
        name="Peroxide (98%)",
        reac_names=(b"H2O2", b"H2O"),
        weights=(0.98, 0.02),
        T_reactant=(298.15, 298.15),
        ox_name="H2O2",
        # Sutton & Biblarz: H2O2/RP-1 L* = 1.52-1.78 m (long; slow decomposition).
        Lstar_range=(1.52, 1.65, 1.78),
    ),
}

# thermo.lib only carries gas-phase 2-propanol, and rocketcea's stock
# "Isopropanol" card also uses the gas heat of formation. Shift both to the
# liquid: dHf(l) = -318.1 kJ/mol vs dHf(g) = -272.7 kJ/mol.
_IPA_H_VAP = 755.5e3  # J/kg
add_new_fuel(
    "Isopropanol(L)",
    """
    fuel C3H8O-2propanol(L) C 3 H 8 O 1   wt%=100.
    h,cal=-76030.     t(k)=298.15   rho=0.786
    """,
)

FUELS: dict[str, Fuel] = {
    "ethanol": Fuel(
        name="Ethanol",
        reac_name=b"C2H5OH(L)",
        T_reactant=298.15,
        fuel_name="Ethanol",
        water_blend=True,
    ),
    "ipa": Fuel(
        name="IPA",
        reac_name=b"C3H8O,2propanol",
        T_reactant=298.15,
        fuel_name="Isopropanol(L)",
        water_blend=True,
        enthalpy_offset=-_IPA_H_VAP,
    ),
    "rp1": Fuel(
        name="RP-1",
        reac_name=b"RP-1",
        T_reactant=298.15,
        fuel_name="RP-1",
    ),
    "lh2": Fuel(
        name="LH2",
        reac_name=b"H2(L)",
        T_reactant=20.15,
        fuel_name="LH2",
    ),
}

# (oxidiser, fuel) pairs whose L* differs from the oxidiser's default band.
_LSTAR_OVERRIDES = {
    # Sutton & Biblarz, Rocket Propulsion Elements: O2/RP-1 L* = 1.02-1.27 m.
    ("lox", "rp1"): (1.02, 1.14, 1.27),
}


@lru_cache(maxsize=None)
def build_propellant(ox_key: str, fuel_key: str, fuel_conc: float = 100.0) -> PropellantConfig:
    """Assemble a PropellantConfig from an oxidiser, a fuel and, for
    water-blendable fuels, the fuel concentration in wt% (balance water)."""
    ox = OXIDISERS[ox_key]
    fuel = FUELS[fuel_key]
    if not fuel.water_blend:
        fuel_conc = 100.0

    fuel_species = [fuel.reac_name]
    fuel_fracs = [fuel_conc / 100]
    fuel_T = [fuel.T_reactant]
    fuel_offsets = [fuel.enthalpy_offset]
    fuel_name = fuel.fuel_name
    name = f"{ox.name} / {fuel.name}"
    if fuel.water_blend:
        name += f" ({fuel_conc:g}%)"
    if fuel_conc < 100:
        fuel_species.append(b"H2O(L)")
        fuel_fracs.append(1 - fuel_conc / 100)
        fuel_T.append(298.15)
        fuel_offsets.append(0.0)
        fuel_name = newFuelBlend([fuel.fuel_name, "H2O"], [fuel_conc, 100 - fuel_conc])

    n_fuel, n_ox = len(fuel_species), len(ox.reac_names)
    return PropellantConfig(
        name=name,
        reac_names=fuel_species + list(ox.reac_names),
        fuel_weights=np.array(fuel_fracs + [0.0] * n_ox),
        oxidant_weights=np.array([0.0] * n_fuel + list(ox.weights)),
        T_reactant=np.array(fuel_T + list(ox.T_reactant)),
        enthalpy_offsets=np.array(fuel_offsets + [0.0] * n_ox),
        ox_name=ox.ox_name,
        fuel_name=fuel_name,
        Lstar_range=_LSTAR_OVERRIDES.get((ox_key, fuel_key), ox.Lstar_range),
    )
