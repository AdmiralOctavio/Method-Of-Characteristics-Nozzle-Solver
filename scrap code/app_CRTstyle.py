
import streamlit as st
import engine.Parameters as P
import engine.solver as solver
import main
import matplotlib.pyplot as plt
import utils.Output as Output

st.set_page_config(page_title="MoC Nozzle Designer", layout="wide")

# ── CRT / Retro Terminal CSS ──────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Share+Tech+Mono&display=swap');

/* ── Root palette ── */
:root {
    --crt-bg:        #0a0500;
    --crt-surface:   #0f0800;
    --crt-border:    #7a3000;
    --crt-amber:     #ff6a00;
    --crt-glow:      #ffaa44;
    --crt-dim:       #7a3000;
    --crt-red:       #ff2200;
    --crt-green:     #88ff00;
    --crt-blue:      #00ccff;
}

/* ── Global background & font ── */
html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"] {
    background-color: var(--crt-bg) !important;
    color: var(--crt-amber) !important;
    font-family: 'Share Tech Mono', monospace !important;
}

/* ── Main content area ── */
[data-testid="stMain"], .main .block-container {
    background-color: var(--crt-bg) !important;
    padding-top: 1.5rem;
}

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background-color: var(--crt-surface) !important;
    border-right: 1px solid var(--crt-border) !important;
}
[data-testid="stSidebar"] * {
    color: var(--crt-amber) !important;
    font-family: 'Share Tech Mono', monospace !important;
}
[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2,
[data-testid="stSidebar"] h3 {
    color: var(--crt-glow) !important;
    border-bottom: 1px solid var(--crt-border);
    padding-bottom: 4px;
    letter-spacing: 0.1em;
}

/* ── Page title ── */
h1 {
    color: var(--crt-glow) !important;
    font-family: 'Share Tech Mono', monospace !important;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    border-bottom: 1px solid var(--crt-border);
    padding-bottom: 0.4rem;
}
h2, h3 {
    color: var(--crt-amber) !important;
    font-family: 'Share Tech Mono', monospace !important;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}

/* ── Containers / panels ── */
[data-testid="stVerticalBlockBorderWrapper"] > div {
    background-color: var(--crt-surface) !important;
    border: 1px solid var(--crt-border) !important;
    border-radius: 0 !important;
    padding: 1rem !important;
}

/* ── Metric widgets ── */
[data-testid="stMetric"] {
    background-color: var(--crt-bg) !important;
    border: 1px solid var(--crt-border) !important;
    border-radius: 0 !important;
    padding: 0.5rem 0.75rem !important;
}
[data-testid="stMetricLabel"] {
    color: var(--crt-dim) !important;
    font-size: 0.7rem !important;
    text-transform: uppercase;
    letter-spacing: 0.1em;
}
[data-testid="stMetricValue"] {
    color: var(--crt-glow) !important;
    font-size: 1.4rem !important;
}
[data-testid="stMetricDelta"] { color: var(--crt-green) !important; }

/* ── Inputs & number inputs ── */
input, textarea, [data-baseweb="input"] input,
[data-baseweb="textarea"] textarea {
    background-color: var(--crt-bg) !important;
    color: var(--crt-amber) !important;
    border: 1px solid var(--crt-border) !important;
    border-radius: 0 !important;
    font-family: 'Share Tech Mono', monospace !important;
    caret-color: var(--crt-amber);
}
input:focus {
    border-color: var(--crt-amber) !important;
    box-shadow: 0 0 6px var(--crt-amber) !important;
    outline: none !important;
}

/* ── Slider ── */
[data-testid="stSlider"] [data-baseweb="slider"] div[role="slider"] {
    background-color: var(--crt-amber) !important;
    border-color: var(--crt-amber) !important;
}
[data-testid="stSlider"] [data-baseweb="slider"] [data-testid="stThumbValue"] {
    color: var(--crt-amber) !important;
}

/* ── Toggle ── */
[data-testid="stToggle"] span {
    background-color: var(--crt-border) !important;
}
[data-testid="stToggle"] input:checked + span {
    background-color: var(--crt-amber) !important;
}

/* ── Buttons ── */
[data-testid="stButton"] button,
button[kind="primary"] {
    background-color: var(--crt-bg) !important;
    color: var(--crt-amber) !important;
    border: 1px solid var(--crt-amber) !important;
    border-radius: 0 !important;
    font-family: 'Share Tech Mono', monospace !important;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    transition: box-shadow 0.15s ease, background-color 0.15s ease;
}
[data-testid="stButton"] button:hover {
    background-color: var(--crt-amber) !important;
    color: var(--crt-bg) !important;
    box-shadow: 0 0 12px var(--crt-amber);
}

/* ── Spinner ── */
[data-testid="stSpinner"] {
    color: var(--crt-amber) !important;
}

/* ── Divider ── */
hr {
    border-color: var(--crt-border) !important;
}

/* ── Alerts: success / error ── */
[data-testid="stAlert"][data-baseweb="notification"] {
    border-radius: 0 !important;
    font-family: 'Share Tech Mono', monospace !important;
}

/* ── Pyplot container ── */
[data-testid="stImage"] img,
[data-testid="stPyplotRootElement"] {
    border: 1px solid var(--crt-border);
}

/* ── Scrollbar ── */
::-webkit-scrollbar { width: 6px; background: var(--crt-bg); }
::-webkit-scrollbar-thumb { background: var(--crt-border); border-radius: 0; }
::-webkit-scrollbar-thumb:hover { background: var(--crt-amber); }

/* ── Hide the sidebar collapse button (shows as raw icon text) ── */
[data-testid="stSidebarCollapsedControl"],
[data-testid="collapsedControl"],
button[kind="headerNoPadding"],
[data-testid="stSidebarNavCollapseIcon"] {
    display: none !important;
}

/* ── Top header / toolbar bar ── */
[data-testid="stHeader"],
header[data-testid="stHeader"] {
    background-color: var(--crt-bg) !important;
    border-bottom: 1px solid var(--crt-border) !important;
}
/* Deploy button and toolbar icons */
[data-testid="stToolbar"] *,
[data-testid="stDecoration"],
[data-testid="stStatusWidget"] {
    color: var(--crt-amber) !important;
    fill: var(--crt-amber) !important;
}
/* The "Deploy" button specifically */
[data-testid="stAppDeployButton"] button {
    background: transparent !important;
    color: var(--crt-amber) !important;
    border: 1px solid var(--crt-border) !important;
    border-radius: 0 !important;
    font-family: 'Share Tech Mono', monospace !important;
}
[data-testid="stAppDeployButton"] button:hover {
    border-color: var(--crt-amber) !important;
    box-shadow: 0 0 8px var(--crt-amber) !important;
}

/* ── Number input +/- buttons ── */
button[data-testid="stNumberInputStepDown"],
button[data-testid="stNumberInputStepUp"] {
    background-color: var(--crt-bg) !important;
    color: var(--crt-amber) !important;
    border: 1px solid var(--crt-border) !important;
    border-radius: 0 !important;
}
button[data-testid="stNumberInputStepDown"]:hover,
button[data-testid="stNumberInputStepUp"]:hover {
    background-color: var(--crt-amber) !important;
    color: var(--crt-bg) !important;
    box-shadow: 0 0 8px var(--crt-amber) !important;
}
/* BaseWeb internal step buttons fallback */
[data-baseweb="input"] + div button,
[data-baseweb="spinner"] button {
    background-color: var(--crt-bg) !important;
    color: var(--crt-amber) !important;
    border-color: var(--crt-border) !important;
    border-radius: 0 !important;
}
/* Force SVG icons inside buttons to amber */
[data-testid="stNumberInputStepDown"] svg,
[data-testid="stNumberInputStepUp"] svg,
[data-baseweb="spinner"] button svg {
    fill: var(--crt-amber) !important;
    color: var(--crt-amber) !important;
}

/* ── CRT scanline overlay (subtle) ── */
[data-testid="stAppViewContainer"]::after {
    content: "";
    pointer-events: none;
    position: fixed;
    inset: 0;
    background: repeating-linear-gradient(
        to bottom,
        transparent 0px,
        transparent 3px,
        rgba(0,0,0,0.18) 3px,
        rgba(0,0,0,0.18) 4px
    );
    z-index: 9999;
}
</style>
""", unsafe_allow_html=True)

# ── Page header ──────────────────────────────────────────────────────────────
st.title("// MOC NOZZLE DESIGNER")
st.caption("METHOD OF CHARACTERISTICS  ·  AXISYMMETRIC NOZZLE SOLVER  ·  T-001V")

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("01 // PERFORMANCE")
    pc_mPa = st.number_input("Chamber Pressure (Bar)", 10, 200, 34)
    P.Ambient_P = st.number_input("Ambient Pressure (Pa)", 0, 101325, 101325)
    of_val = st.number_input("O/F Ratio", value=5.13)
    thrust_val = st.number_input("Target Thrust (N)", value=600)
    P.Refinement = st.number_input("Refinement Factor", value=100)
    P.Shorten_Percentage = st.slider("Total Nozzle Length", 50, 100, 75, 5) / 100

    st.header("02 // CHAMBER")
    P.L_combustion = st.number_input("Total Chamber Length", value=83.02)
    P.Contraction_ratio = st.number_input("Contraction Ratio", value=16)
    P.Chamber_Slope = st.number_input("Convergent Chamber Slope", value=45)
    P.R1 = st.number_input("Radius 1", value=10)
    P.R2 = st.number_input("Radius 2", value=50)

    st.header("03 // OUTPUT")
    P.Graph2d = st.toggle("Show 2D Characteristic Grid", value=True)
    P.Graph3d = st.toggle("Show 3D Characteristic Grid", value=False)
    P.Stl = st.toggle("Generate STL", value=False)
    P.Dxf = st.toggle("Generate DXF", value=True)

# ── Run button ────────────────────────────────────────────────────────────────
if st.sidebar.button("▶  RUN DESIGN OPTIMIZATION", use_container_width=True):
    P.P_combustion = pc_mPa * 1e5
    P.Oxidiser_Fuel_Ratio = of_val
    P.Thrust = thrust_val

    P.update_engine_data(P.P_combustion, P.Oxidiser_Fuel_Ratio)

    with st.spinner("CONVERGING ON GEOMETRY..."):
        opt_results = main.run(gui_mode=True)
        data = Output.outputTable(opt_results['rt'], opt_results['mdot'], opt_results['mach'])

    # ── Dashboard ─────────────────────────────────────────────────────────────
    st.header("FINAL NOZZLE DESIGN SPECIFICATIONS")

    # Panel 1: Geometry
    with st.container(border=True):
        st.subheader("GEOMETRY + DIMENSIONS")
        c1, c2, c3 = st.columns(3)
        c1.metric("Nozzle Length",  f"{data['wall_x'][-1]:.2f} mm")
        c1.metric("Total Length",   f"{data['total_length']:.2f} mm")
        c2.metric("Throat Radius",  f"{data['y_min']:.2f} mm")
        c2.metric("Throat Diameter",f"{data['y_min']*2:.2f} mm")
        c3.metric("Exit Radius",    f"{data['exit_radius']:.2f} mm")
        angle_color = "normal" if data['Exit_Angle'] <= 6 else "inverse"
        st.metric("Exit Angle", f"{data['Exit_Angle']:.2f}°", delta_color=angle_color)

    # Panel 2: Performance
    with st.container(border=True):
        st.subheader("PERFORMANCE + COMBUSTION")
        p1, p2, p3 = st.columns(3)
        p1.metric("Chamber Temperature", f"{data['T_comb']:.1f} K")
        p1.metric("Gamma",               f"{data['g']:.3f}")
        p2.metric("Optimal Pressure Ratio", f"{data['P_comb']/101325:.2f}")
        p2.metric("Optimal Exit Mach",      f"{data['M_opt']:.2f}")
        p2.metric("Optimal Expansion Ratio",f"{data['AR_opt']:.2f}")
        p3.metric("Predicted Exit Mach",    f"{data['M_exit_predicted']:.2f}")

    # Panel 3: Thrust
    with st.container(border=True):
        st.subheader("THRUST + FLOW STABILITY")
        t1, t2 = st.columns(2)
        with t1:
            st.metric("Predicted Thrust", f"{data['Thrust_total']:.1f} N")
        with t2:
            st.metric("Exit Pressure", f"{data['P_exit']/1000:.1f} kPa")
            if data['P_exit'] < 0.4 * P.Ambient_P:
                st.error("⚠  FLOW SEPARATION RISK HIGH")
            else:
                st.success("✅  FLOW STABLE")

    # Panel 4: Efficiency
    with st.container(border=True):
        st.subheader("EFFICIENCY + MASS FLOW")
        e1, e2 = st.columns(2)
        e1.metric("Specific Impulse (Design)", f"{data['Isp_design']:.1f} s")
        e1.write(f"**CEA Ideal Isp:** {data['Isp_cea']:.1f} s")
        e2.metric("Mass Flow Rate", f"{data['mdot']:.3f} kg/s")

    # ── Plot ──────────────────────────────────────────────────────────────────
    if data['fig']:
        st.divider()
        st.subheader("METHOD OF CHARACTERISTICS GRID")
        st.pyplot(data['fig'])




