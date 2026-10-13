from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.rule import Rule
import engine.Parameters as Param
import numpy as np
import csv
import utils.stlgenerator as stlgenerator
import os
import engine.solver_lcars as solver_lcars
import engine.solver_crt as solver_crt
import engine.solver_standard as solver_standard


def m(s):
    return Text.from_markup(s)


def divider(color):
    return [Rule(style=color)]


folder_name = "Output_Files"
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

Graph2d = Param.Graph2d
Graph3d = Param.Graph3d
Graph3d_Fancy = Param.Graph3d_Fancy
Write = Param.Write
Stl = Param.Stl
Dxf = Param.Dxf

'''
Output file. Just outputs a table with the important design parameters, along with activating the write functions if required.

'''

def write_params_csv(data, filepath):
    """Write a params.csv holding every inputted and outputted parameter.

    Each row is (category, parameter, value, unit) so the file is
    self-documenting and easy to diff between runs.
    """
    if Param.Lstar_Sizing:
        chamber_len_row = ("Input", "Chamber L*", Param.Lstar_Value, "m")
    else:
        chamber_len_row = ("Input", "Chamber Length", Param.L_combustion, "mm")

    rows = [
        # ── Inputs ────────────────────────────────────────────────────────────
        ("Input", "Propellant", Param.get_propellant().name, "-"),
        ("Input", "Chamber Pressure", Param.P_combustion, "Pa"),
        ("Input", "Ambient Pressure", Param.Ambient_P, "Pa"),
        ("Input", "O/F Ratio", Param.Oxidiser_Fuel_Ratio, "-"),
        ("Input", "Target Thrust", Param.Thrust, "N"),
        ("Input", "Combustion Efficiency", Param.Combustion_Efficiency, "-"),
        ("Input", "Refinement Factor", Param.Refinement, "-"),
        ("Input", "Nozzle Length Fraction", Param.Shorten_Percentage, "-"),
        ("Input", "L* Sizing", Param.Lstar_Sizing, "bool"),
        chamber_len_row,
        ("Input", "Contraction Ratio", Param.Contraction_ratio, "-"),
        ("Input", "Chamber Slope", Param.Chamber_Slope, "deg"),
        ("Input", "Radius 1", Param.R1, "mm"),
        ("Input", "Radius 2", Param.R2, "mm"),
        # ── Outputs ───────────────────────────────────────────────────────────
        ("Output", "Nozzle Length", data["wall_x"][-1], "mm"),
        ("Output", "Total Length", data["total_length"], "mm"),
        ("Output", "Throat Radius", data["y_min"], "mm"),
        ("Output", "Throat Diameter", 2 * data["y_min"], "mm"),
        ("Output", "Chamber Diameter",
         2 * data["y_min"] * np.sqrt(Param.Contraction_ratio), "mm"),
        ("Output", "Exit Radius", data["exit_radius"], "mm"),
        ("Output", "Exit Angle", data["Exit_Angle"], "deg"),
        ("Output", "Combustion Temperature", data["T_comb"], "K"),
        ("Output", "Gamma", data["g"], "-"),
        ("Output", "Optimal Exit Mach", data["M_opt"], "-"),
        ("Output", "Optimal Expansion Ratio", data["AR_opt"], "-"),
        ("Output", "Theoretical Expansion Ratio", data["AR_theoretical"], "-"),
        ("Output", "True Expansion Ratio", data["AR_true"], "-"),
        ("Output", "Design Exit Mach", data["M_exit_design"], "-"),
        ("Output", "Predicted Exit Mach", data["M_exit_predicted"], "-"),
        ("Output", "Total Thrust", data["Thrust_total"], "N"),
        ("Output", "Momentum Thrust", data["Thrust_momentum"], "N"),
        ("Output", "Pressure Thrust", data["Thrust_pressure"], "N"),
        ("Output", "Exit Pressure", data["P_exit"], "Pa"),
        ("Output", "Exit Velocity", data["Ve"], "m/s"),
        ("Output", "Specific Impulse (Design)", data["Isp_design"], "s"),
        ("Output", "Specific Impulse (CEA)", data["Isp_cea"], "s"),
        ("Output", "Mass Flow Rate", data["mdot"], "kg/s"),
    ]

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Category", "Parameter", "Value", "Unit"])
        writer.writerows(rows)


def outputTable(R, mdot, mach, theme):
    console = Console()

    if theme == 'lcars':

        (
            wall_x,
            Exit_Angle,
            y_min,
            T_combustion,
            g,
            P_combustion,
            M_optimal,
            AreaRatio,
            wall_y,
            M_exit,
            M_exit_characteristic,
            Thrust,
            A_exit,
            P_exit,
            Ve,
            ISP_cea,
            y_calc,
            L,
            wall_y_mirrored,
            mdot,
            fig,
        ) = solver_lcars.main(mdot, R, mach)
    elif theme == 'crt': 
        (
            wall_x,
            Exit_Angle,
            y_min,
            T_combustion,
            g,
            P_combustion,
            M_optimal,
            AreaRatio,
            wall_y,
            M_exit,
            M_exit_characteristic,
            Thrust,
            A_exit,
            P_exit,
            Ve,
            ISP_cea,
            y_calc,
            L,
            wall_y_mirrored,
            mdot,
            fig,
        ) = solver_crt.main(mdot, R, mach)

    elif theme == 'standard':
        (
            wall_x,
            Exit_Angle,
            y_min,
            T_combustion,
            g,
            P_combustion,
            M_optimal,
            AreaRatio,
            wall_y,
            M_exit,
            M_exit_characteristic,
            Thrust,
            A_exit,
            P_exit,
            Ve,
            ISP_cea,
            y_calc,
            L,
            wall_y_mirrored,
            mdot,
            fig,
        ) = solver_standard.main(mdot, R, mach)

        

    table = Table(
        show_header=False,
        box=None,
        pad_edge=True,
        expand=False,
        padding=(0, 2),
        collapse_padding=False,
    )
    table.add_column(max_width=50, justify="left")  # label 1
    table.add_column(max_width=50, justify="right")  # value 1
    table.add_column(min_width=5)  # SPACER COLUMN
    table.add_column(max_width=50, justify="left")  # label 2
    table.add_column(max_width=50, justify="right")  # value 2

    table.add_row(
        "[cyan3]Nozzle length:",
        f"[cyan3]{wall_x[-1]:.2f} mm",
        " ",
        "[cyan3]Total length:",
        f"[cyan3]{wall_x[-1] - wall_x[0]:.2f} mm",
    )

    exit_style = "red" if Exit_Angle > 6 else "dark_turquoise"
    table.add_row(
        m("[dark_turquoise]Exit Angle:[/dark_turquoise]"),
        m(f"[{exit_style}]{Exit_Angle:.2f} deg"),
        " ",
        m("[dark_turquoise]Exit radius:"),
        f"[dark_turquoise]{wall_y[-1]:.2f} mm",
    )

    table.add_row(
        "[cyan3]True throat radius:",
        f"[cyan3]{y_min:.2f} mm",
        " ",
        "[cyan3]True throat diameter:",
        f"[cyan3]{2 * y_min:.2f} mm",
    )

    table.add_row(
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
    )

    table.add_row(
        "[orange_red1]Combustion temperature:",
        f"[orange_red1]{T_combustion:.1f} K",
        " ",
        "[orange_red1]Gamma:",
        f"[orange_red1]{g:.3f}",
    )

    table.add_row(
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
    )

    table.add_row(
        "[light_sky_blue3]Optimal pressure ratio:",
        f"[light_sky_blue3]{P_combustion / 101325:.2f}",
        "",
        "",
    )
    table.add_row(
        "[sky_blue2]Optimal exit Mach:", f"[sky_blue2]{M_optimal:.2f}", "", ""
    )
    table.add_row(
        "[light_sky_blue3]Optimal expansion ratio:",
        f"[light_sky_blue3]{AreaRatio:.2f}",
        "",
        "",
    )

    table.add_row(
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
    )

    table.add_row(
        "[green_yellow]Theoretical expansion ratio:",
        f"[green_yellow]{(wall_y[-1] ** 2 / L**2):.2f}",
        " ",
        "[green_yellow]True expansion ratio:",
        f"[green_yellow]{(y_calc**2 / y_min**2):.2f}",
    )

    table.add_row(
        "[light_green]Design exit Mach:",
        f"[light_green]{M_exit:.2f}",
        " ",
        "[light_green]Predicted exit Mach:",
        f"[light_green]{M_exit_characteristic:.2f}",
    )

    table.add_row(
        "[green_yellow]Predicted thrust:", f"[green_yellow]{Thrust:.0f} N", "", ""
    )

    table.add_row(
        "[light_green]> Thrust (massflow):",
        f"[light_green]{Thrust - (P_exit - Param.Ambient_P) * A_exit:.2f} N",
        " ",
        "[light_green]> Thrust (pressure):",
        f"[light_green]{(P_exit - Param.Ambient_P) * A_exit:.2f} N",
    )

    if P_exit < 0.25 * Param.Ambient_P:
        warn = Text("FLOW SEPARATION WILL OCCUR", style="bold red")
    elif P_exit < 0.4 * Param.Ambient_P:
        warn = Text("Flow separation may occur", style="bold orange1")
    elif P_exit < 0.5 * Param.Ambient_P:
        warn = Text("Nearing exit instability region", style="yellow")
    else:
        warn = Text("Flow stable", style="green")

    table.add_row(
        "[green_yellow]Exit pressure:", f"[green_yellow]{P_exit:.0f} Pa", "", ""
    )
    table.add_row(warn, "", "", "")

    table.add_row(
        "[light_green]Specific impulse:",
        f"[light_green]{Thrust / (mdot * 9.80665):.2f} s",
        " ",
        "[light_green]CEA ISp:",
        f"[light_green]{ISP_cea:.2f} s",
    )
    table.add_row("[green_yellow]Mass flow-rate:", f"[green_yellow]{mdot:.3f} kg/s")

    panel = Panel(
        table,
        title="[bold red]Output Nozzle Design Specifications[/bold red]",
        border_style="cyan",
        padding=(1, 2),
        expand=True,
    )

    console.print(panel)

    if Write:
        wall_z = np.zeros_like(wall_x)
        combined = np.stack((wall_x, wall_y), axis=1)
        filename = f"Contour.csv"
        filepath = os.path.join(folder_name, filename)
        np.savetxt(
            filepath,
            combined,
            delimiter=",",
            comments="",
        )

    if Stl:
        stlgenerator.create_stl(wall_x, wall_y, M_exit)

    if Dxf:
        stlgenerator.create_dxf(wall_x / 1000, wall_y_mirrored / 1000, M_exit)

    result = {
        "wall_x": wall_x,
        "wall_y": wall_y,
        "Exit_Angle": Exit_Angle,
        "y_min": y_min,
        "T_comb": T_combustion,
        "g": g,
        "P_comb": P_combustion,
        "M_opt": M_optimal,
        "AR_opt": AreaRatio,
        "M_exit_design": M_exit,
        "M_exit_predicted": M_exit_characteristic,
        "Thrust_total": Thrust,
        "Thrust_momentum": mdot * Ve,
        "Thrust_pressure": (P_exit - Param.Ambient_P) * A_exit,
        "P_exit": P_exit,
        "Ve": Ve,
        "Isp_cea": ISP_cea,
        "Isp_design": Ve / 9.80665,
        "mdot": mdot,
        "AR_theoretical": (wall_y[-1] ** 2 / L**2),
        "AR_true": (y_calc**2 / y_min**2),
        "fig": fig,
        "total_length": wall_x[-1] - wall_x[0],
        "exit_radius": wall_y[-1],
    }

    if Write:
        params_path = os.path.join(folder_name, "params.csv")
        write_params_csv(result, params_path)

    return result
