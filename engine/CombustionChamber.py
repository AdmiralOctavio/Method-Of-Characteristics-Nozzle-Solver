import numpy as np
from scipy.integrate import trapezoid
import engine.Parameters as P
from engine.Fuels import PROPELLANTS


def _convergent_volume(R1, R2, Rt, alpha, R_inlet, L_straight, L_total):
    """Volume (mm^3) of the convergent section, inlet plane (R_inlet) -> throat (Rt).

    Integrates pi * y(x)^2 dx over the full convergent wall profile
    (inlet arc R2 -> conical line -> throat arc R1), which is fully determined
    by the contraction geometry and independent of the chamber length.
    """
    theta_inlet = np.linspace(0, alpha, 200)
    x_in = -L_total + R2 * np.sin(theta_inlet)
    y_in = R_inlet - R2 * (1 - np.cos(theta_inlet))

    x_t1, y_t1 = x_in[-1], y_in[-1]
    x_t2 = x_t1 + L_straight * np.cos(alpha)
    y_t2 = y_t1 - L_straight * np.sin(alpha)

    theta_throat = np.linspace(alpha, 0, 200)
    x_th = -R1 * np.sin(theta_throat)
    y_th = Rt + R1 * (1 - np.cos(theta_throat))

    x = np.concatenate([x_in, [x_t2], x_th])
    y = np.concatenate([y_in, [y_t2], y_th])
    return trapezoid(np.pi * y ** 2, x)


def _length_from_lstar(R1, R2, Rt, alpha, R_inlet, L_straight, L_total):
    """Total chamber length (mm) needed to hit the propellant's target L*.

    L* defines a chamber *volume* (V_c = L* * A_throat). The convergent volume
    is fixed by geometry, so the remaining volume sets the cylindrical barrel
    length. Total length = convergent axial length + barrel length.
    """
    Lstar_m = getattr(P, "Lstar_Value", None)
    if Lstar_m is None:  # legacy fallback to the discrete short/nominal/elongated option
        prop = PROPELLANTS[P.Selected_Propellant]
        Lstar_m = prop.get_lstar(getattr(P, "Lstar_Option", "nominal"))
    Lstar_mm = Lstar_m * 1000.0  # m -> mm

    A_throat = np.pi * Rt ** 2
    V_target = Lstar_mm * A_throat
    V_conv = _convergent_volume(R1, R2, Rt, alpha, R_inlet, L_straight, L_total)

    L_barrel = (V_target - V_conv) / (np.pi * R_inlet ** 2)
    if L_barrel < 0:
        print(
            f"Warning: L* ({Lstar_mm/1000:.3f} m) is too small for this contraction "
            f"geometry — the convergent section alone exceeds the target volume. "
            f"Clamping barrel length to 0. Reduce the contraction ratio or pick a "
            f"longer L* option."
        )
        L_barrel = 0.0
    return L_total + L_barrel


def CombustionChamber(wall_x, wall_y, R1, R2, L_throat):
    alpha = np.deg2rad(P.Chamber_Slope)
    Rt = L_throat
    Contraction_Ratio = P.Contraction_ratio
    R_inlet = Rt * np.sqrt(Contraction_Ratio)

    dy_straight = (R_inlet - Rt) - (R1 + R2) * (1 - np.cos(alpha))
    L_straight = dy_straight / np.sin(alpha)

    L_total = (R1 + R2) * np.sin(alpha) + L_straight * np.cos(alpha)

    if getattr(P, "Lstar_Sizing", False):
        Length = _length_from_lstar(R1, R2, Rt, alpha, R_inlet, L_straight, L_total)
        P.L_combustion = Length  # expose computed length for output / printouts
    else:
        Length = P.L_combustion

    theta_inlet = np.linspace(0, alpha, 50)
    x_arc_inlet = -L_total + R2 * np.sin(theta_inlet)
    y_arc_inlet = R_inlet - R2 * (1 - np.cos(theta_inlet))

    x_line_cont = -Length
    y_line_cont = R_inlet

    x_t1, y_t1 = x_arc_inlet[-1], y_arc_inlet[-1]
    x_t2 = x_t1 + L_straight * np.cos(alpha)
    y_t2 = y_t1 - L_straight * np.sin(alpha)
    x_line = np.array([x_t1, x_t2])
    y_line = np.array([y_t1, y_t2])

    theta_throat = np.linspace(alpha, 0, 50)
    x_arc_throat = -R1 * np.sin(theta_throat)
    y_arc_throat = Rt + R1 * (1 - np.cos(theta_throat))

    wall_x = np.append(x_arc_throat, wall_x)
    wall_x = np.append(x_line, wall_x)
    wall_x = np.append(x_arc_inlet, wall_x)
    wall_x = np.append(x_line_cont, wall_x)

    wall_y = np.append(y_arc_throat, wall_y)
    wall_y = np.append(y_line, wall_y)
    wall_y = np.append(y_arc_inlet, wall_y)
    wall_y = np.append(y_line_cont, wall_y)

    return wall_x, wall_y
