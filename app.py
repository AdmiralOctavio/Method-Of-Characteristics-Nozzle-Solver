import streamlit as st
import engine.Parameters as P
import engine.solver as solver
import main
import matplotlib.pyplot as plt
import utils.Output as Output

st.set_page_config(page_title="MoC Nozzle Designer", layout="wide")

st.title("Method of Characteristics Nozzle Solver")

with st.sidebar:
    st.header("1. Performance Goals")
    pc_mPa = st.number_input("Chamber Pressure (Bar)", 10, 200, 34)
    of_val = st.number_input("O/F Ratio", value=5.13)
    thrust_val = st.number_input("Target Thrust (N)", value=600)
    P.Refinement = st.number_input("Refinement Factor", value = 100)
    P.Shorten_Percentage = st.slider("Shorten Percentage", 0.5, 1.0, 0.75, 0.05)

    st.header("2. Chamber Design")
    P.L_combustion = st.number_input("Total chamber length", value = 83.02)
    P.Contraction_ratio = st.number_input("Contraction Ratio", value = 16)
    P.Chamber_Slope = st.number_input("Convergent Chamber Slope", value = 45)
    P.R1 = st.number_input("Radius 1", value = 10)
    P.R2 = st.number_input("Radius 2", value = 50)
    
    st.header("3. Plotting & Files")
    P.Graph2d = st.toggle("Show 2D Characteristic Grid", value=True)
    P.Graph3d = st.toggle("Show 3D Characteristic Grid", value=False)
    P.Stl = st.toggle("Generate STL", value=False)
    P.Dxf = st.toggle("Generate DXF", value=True)

if st.sidebar.button("Run Design Optimization", use_container_width=True):
    P.P_combustion = pc_mPa * 1e5
    P.Oxidiser_Fuel_Ratio = of_val
    P.Thrust = thrust_val
    
    P.update_engine_data(P.P_combustion, P.Oxidiser_Fuel_Ratio)
    
    with st.spinner("Converging on geometry..."):
        # 1. Run optimization logic
        opt_results = main.run(gui_mode=True)
        # 2. Get the dictionary of all parameters we just mapped in Output.py
        data = Output.outputTable(opt_results['rt'], opt_results['mdot'], opt_results['mach'])

    # --- TOP LEVEL DASHBOARD ---
    st.header("Final Nozzle Design Specifications")
    
    # Panel 1: Dimensions
    with st.container(border=True):
        st.subheader("Geometry & Dimensions")
        c1, c2, c3 = st.columns(3)
        c1.metric("Nozzle Length", f"{data['wall_x'][-1]:.2f} mm")
        c1.metric("Total Length", f"{data['total_length']:.2f} mm")
        
        c2.metric("Throat Radius", f"{data['y_min']:.2f} mm")
        c2.metric("Throat Diameter", f"{data['y_min']*2:.2f} mm")
        
        c3.metric("Exit Radius", f"{data['exit_radius']:.2f} mm")
        # Recreating your red/turquoise logic for exit angle
        angle_color = "normal" if data['Exit_Angle'] <= 6 else "inverse"
        st.metric("Exit Angle", f"{data['Exit_Angle']:.2f}°", delta_color=angle_color)

    # Panel 2: Performance & Chemistry
    with st.container(border=True):
        st.subheader("Performance & Combustion")
        p1, p2, p3 = st.columns(3)
        p1.metric("Chamber Temperature:", f"{data['T_comb']:.1f} K")
        p1.metric("Gamma:", f"{data['g']:.3f}")
        
        p2.metric("Optimal Pressure Ratio:",  f"{data['P_comb']/101325:.2f}")
        p2.metric("Optimal Exit Mach:",  f" {data['M_opt']:.2f}")
        p2.metric("Optimal Expansion Ratio:",  f" {data['AR_opt']:.2f}")
        
        p3.metric("Design Exit Mach:",  f" {data['M_exit_design']:.2f}")
        p3.metric("Predicted Exit Mach:",  f" {data['M_exit_predicted']:.2f}")

    # Panel 3: Thrust & Stability (Red/Green Logic)
    with st.container(border=True):
        st.subheader("Thrust & Flow Stability")
        t1, t2 = st.columns(2)
        with t1:
            st.metric("Total Predicted Thrust", f"{data['Thrust_total']:.1f} N")
            st.write(f"**Momentum Component:** {data['Thrust_momentum']:.2f} N")
            st.write(f"**Pressure Component:** {data['Thrust_pressure']:.2f} N")
        
        with t2:
            st.metric("Exit Pressure", f"{data['P_exit']/1000:.1f} kPa")
            # Stability status
            if data['P_exit'] < 0.4 * 101325:
                st.error("⚠️ FLOW SEPARATION RISK HIGH")
            else:
                st.success("✅ Flow Stable")

    # Panel 4: Efficiency
    with st.container(border=True):
        st.subheader("Efficiency & Mass Flow")
        e1, e2 = st.columns(2)
        e1.metric("Specific Impulse (Design)", f"{data['Isp_design']:.1f} s")
        e1.write(f"**CEA Ideal Isp:** {data['Isp_cea']:.1f} s")
        e2.metric("Mass Flow Rate", f"{data['mdot']:.3f} kg/s")

    # Visualization
    if data['fig']:
        st.divider()
        st.subheader("Method of Characteristics Grid")
        st.pyplot(data['fig'])