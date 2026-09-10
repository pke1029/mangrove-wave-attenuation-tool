import streamlit as st
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Set the title and favicon that appear in the Browser's tab bar.
st.set_page_config(
    page_title='Mangrove Wave Attenuation Tool',
    page_icon=':ocean:', # This is an emoji shortcode. Could be a URL too.
)

"# :ocean: Mangrove Wave Attenuation Tool v0.1"

st.image("https://images.unsplash.com/photo-1589556183130-530470785fab?q=80&w=1170&auto=format&fit=crop&ixlib=rb-4.1.0&ixid=M3wxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8fA%3D%3D")

'''
Mangrove Wave Attenuation Tool is a simple interactive calculator that 
estimates wave attenuation through mangrove forests given wave input 
parameters. It helps visualize how mangroves reduce wave energy and 
support coastal resilience. The wave attenuation model is based on the 
work by Pang and Tay (https://arxiv.org/abs/2606.11653).
'''

"### Wave conditions"
col1, col2, col3 = st.columns(3)
water_depth = col1.number_input("Water Depth h [m]", value=0.5)
wave_period = col2.number_input("Wave Period T [s]", value=1.0)
wave_height = col3.number_input("Wave Height H [m]", value=0.1)
wave_angular_frequency = 2*np.pi/wave_period
wave_number = fsolve(lambda k: 9.81*k*np.tanh(k*water_depth) - wave_angular_frequency**2, 1.0)
wave_number = wave_number[0]
wave_length = 2*np.pi / wave_number
wave_length = col1.number_input("Wavelength L [m]", value=wave_length, disabled=True)
wave_steepness =  col2.number_input("Wave Steepness H/L [-]", value=wave_height/wave_length, disabled=True)
relative_water_depth = col3.number_input("Relative Water Depth h/L [-]", value=water_depth/wave_length, disabled=True)

"### Mangrove root properties"

"Below are suggested values obtained from literatures. Using site specific values (especially the root density which has the highest variability) yield better results."

# is_enabled = st.toggle("Enable edit", value=False)
df = pd.DataFrame({
    "Species/Genus":["Rhizophora [1,3]", "Sonneratia [1,2]", "Avicennia [4]"],
    "Root density [1/m²]": [150, 500, 250],
    "Average root height [m]": [0.6, 0.06, 0.07],
    "Average root diameter [m]":[0.025, 0.005, 0.006],
    # "Enabled": [True, False, False]
})
df = st.data_editor(df, disabled=False, hide_index=True, num_rows="dynamic")

species = df["Species/Genus"]
df_distribution = pd.DataFrame(columns=species)
z = np.linspace(0, 2, 101)
for row in df.itertuples(index=False):
    df_distribution[row[0]] = row[1] * row[3] * np.exp(-z/row[2]) 
df_distribution["Elevation [m]"] = z
st.line_chart(df_distribution, x="Elevation [m]", y_label="Frontal cover [m/m²]") 
st.caption("Plot shows the frontal area distribution of each species of mangrove in terms of elevation.", text_alignment="center")

show_source = st.toggle("Show source", value=True)
if show_source:
    st.markdown(
    """
    <small>
    [1] Horstman, Erik M., et al. "Wave attenuation in mangroves: A quantitative approach to field observations." Coastal engineering 94 (2014): 47-62. <br>
    [2] Liénard, Jean, et al. "Efficient three-dimensional reconstruction of aquatic vegetation geometry: Estimating morphological parameters influencing hydrodynamic drag." Estuarine, Coastal and Shelf Science 178 (2016): 77-85. <br>
    [3] Mori, Nobuhito, et al. "Parameterization of mangrove root structure of Rhizophora stylosa in coastal hydrodynamic model." Frontiers in Built Environment 7 (2022): 782219. <br>
    [4] Horstman, Erik M., et al. "Are flow-vegetation interactions well represented by mimics? A case study of mangrove pneumatophores." Advances in water resources 111 (2018): 360-371.
    </small>
    """, unsafe_allow_html=True
    )

"### Estimated wave attenuation"

# "Here we used the drag coefficient from Mendez and Losada (2004)."
# st.latex(r"C_D = 0.47\mathrm{e}^{-0.052K_C}, \qquad K_C=\frac{uT}{D}")
drag_coefficient = st.number_input("Drag Coefficient [-]", value=1.0)

# normalised KD
def KD(k, h, a0, b):
    kh = k*h
    bh = b*h
    abar = a0/bh * (1.0 - np.exp(-bh))
    f = lambda x: (np.exp(x)-1.0)/x # if np.abs(x)>0.001 else 0.5*x+1.0
    phi = a0*h/(8*np.sinh(kh)**3) * (f(3*kh-bh) + 3*f(kh-bh) + 3*f(-kh-bh) + f(-3*kh-bh))
    val = 2*k**2*np.tanh(kh)**2 / (kh + np.tanh(kh) - kh*np.tanh(kh)**2) * phi
    val = val / (3*np.pi) * drag_coefficient
    return val 

df_KD = pd.DataFrame({"Species/Genus":df["Species/Genus"]})
df_KD["Wave dacay coefficient [1/m²]"] = KD(wave_number, water_depth, df["Root density [1/m²]"]*df["Average root diameter [m]"], 1/df["Average root height [m]"])
df_KD.loc[len(df_KD)] = ["Combined", df_KD["Wave dacay coefficient [1/m²]"].sum()]
st.table(df_KD)

x_range = st.slider("Mangrove belt width [m]", 100, 1000, 500)

df_attenuation = pd.DataFrame(columns=species)
x = np.linspace(0, x_range, 101)
for row in df_KD.itertuples(index=False):
    df_attenuation[row[0]] = 1/(1 + row[1] * wave_height * x)
df_attenuation["Distance [m]"] = x
st.line_chart(df_attenuation, x="Distance [m]", y_label="Wave attenuation factor [-]", x_label="Distance along mangrove belt [m]") 

col1, col2, col3 = st.columns(3)
col1.metric("Incident wave height", f"{wave_height:.3f} m", border=True)
col2.metric("Transmitted wave height", f"{df_attenuation["Combined"].iloc[-1]*wave_height:.3f} m", border=True)
col3.metric("Combined height reduction", f"{1-df_attenuation["Combined"].iloc[-1]:.2f} %", border=True)

"### Sensitivity to sea-level rise"

def slr_sensitivity(k, h, b):
    kh = k*h
    bh = b*h
    m = 2.0 - bh/(np.exp(bh)-1.0) - kh*(3*np.cosh(kh)*np.sinh(kh)+kh)/(np.sinh(kh)**2 + kh*np.tanh(kh))
    T = lambda f, kh, bh: f(3*kh-bh) + 3*f(kh-bh) + 3*f(-kh-bh) + f(-3*kh-bh)
    m +=  T(lambda x: np.exp(x), kh, bh) / T(lambda x: (np.exp(x)-1.0)/x, kh, bh)
    return m

df["SLR sensitivity [-]"] = slr_sensitivity(wave_number, water_depth, 1/df["Average root height [m]"])
df["Efficiency change per cm SLR [%/cm]"] = df["SLR sensitivity [-]"] / water_depth
st.table(df[["Species/Genus", "SLR sensitivity [-]", "Efficiency change per cm SLR [%/cm]"]])
