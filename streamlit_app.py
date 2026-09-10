import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from scipy.optimize import fsolve

# Set the title and favicon that appear in the Browser's tab bar.
st.set_page_config(
    page_title='Mangrove Wave Attenuation Tool',
    page_icon=':ocean:', # This is an emoji shortcode. Could be a URL too.
)

"# :ocean: Mangrove Wave Attenuation Tool v0.2"

'''
Estimate how much a mangrove belt attenuates an incoming wave, from root structure and local water conditions. 
The wave attenuation model is based on the work by Pang and Tay (https://arxiv.org/abs/2606.11653).
'''

st.divider()

col1, col2 = st.columns(2)

col1.write("### :material/waves: Wave conditions")

water_depth = col1.slider("Water depth, h [m]", 0.1, 5.0, 0.5)
wave_height = col1.slider("Wave height, H [m]", 0.1, 2.0, 0.3)
wave_period = col1.slider("Wave period, T [s]", 1.0, 8.0, 3.0)
x_range = col1.slider("Mangrove belth width [m]", 10, 1000, 100, step=10)

wave_angular_frequency = 2*np.pi/wave_period
wave_number = fsolve(lambda k: 9.81*k*np.tanh(k*water_depth) - wave_angular_frequency**2, 1.0)
wave_number = wave_number[0]

col1.write("### :material/nest_eco_leaf: Vegetations")

col11, col12 = col1.columns(2)
if col11.button("Rhizophora [1,2]", width="stretch"):
    st.session_state.root_diameter = 25
    st.session_state.root_height = 0.6
    st.session_state.root_density = 150
if col12.button("Sonneratia [1,3]", width="stretch"):
    st.session_state.root_diameter = 5
    st.session_state.root_height = 0.06
    st.session_state.root_density = 500
if col11.button("Avicennia [4]", width="stretch"):
    st.session_state.root_diameter = 6
    st.session_state.root_height = 0.07
    st.session_state.root_density = 250

root_diameter = col1.slider("Average root diameter [mm]", 0, 100, 25, key="root_diameter") 
root_height = col1.slider("Average root height [m]", 0.01, 2.0, 0.6, key="root_height")
root_density = col1.slider("Root density [1/m²]", 0, 2000, 150, key="root_density", step=10)
drag_coefficient = col1.slider("Drag coefficient [-]", 0.0, 6.0, 1.0)

col1.markdown(
"""
<small>
[1] Horstman, Erik M., et al. "Wave attenuation in mangroves: A quantitative approach to field observations." Coastal engineering 94 (2014): 47-62. <br>
[2] Mori, Nobuhito, et al. "Parameterization of mangrove root structure of Rhizophora stylosa in coastal hydrodynamic model." Frontiers in Built Environment 7 (2022): 782219. <br>
[3] Liénard, Jean, et al. "Efficient three-dimensional reconstruction of aquatic vegetation geometry: Estimating morphological parameters influencing hydrodynamic drag." Estuarine, Coastal and Shelf Science 178 (2016): 77-85. <br>
[4] Horstman, Erik M., et al. "Are flow-vegetation interactions well represented by mimics? A case study of mangrove pneumatophores." Advances in water resources 111 (2018): 360-371.
</small>
""", unsafe_allow_html=True
)

def KD(k, h, a0, b):
    kh = k*h
    bh = b*h
    abar = a0/bh * (1.0 - np.exp(-bh))
    f = lambda x: (np.exp(x)-1.0)/x # if np.abs(x)>0.001 else 0.5*x+1.0
    phi = a0*h/(8*np.sinh(kh)**3) * (f(3*kh-bh) + 3*f(kh-bh) + 3*f(-kh-bh) + f(-3*kh-bh))
    val = 2*k**2*np.tanh(kh)**2 / (kh + np.tanh(kh) - kh*np.tanh(kh)**2) * phi
    val = val / (3*np.pi) * drag_coefficient
    return val 

wave_decay_coefficient = KD(wave_number, water_depth, 0.001*root_diameter*root_density, 1/root_height)
col2.metric(
    "Wave decay coefficient", 
    f"{wave_decay_coefficient:.3f} m⁻²", 
    border=True, 
    help="The wave height across the beltfollows the equation H/(1+KHx)"
)

x = np.linspace(0, x_range, 101)
y = wave_height/(1+wave_decay_coefficient*wave_height*x)
df = pd.DataFrame({"x":x, "y":y})
# col2.area_chart(df, x="x", y="y", x_label="Distance into belt [m]", y_label="", color="#008CFF37") 
fig = px.area(df, x="x", y="y")
fig.update_layout(xaxis_title="Distance into belt [m]")
fig.update_layout(yaxis_title="")
fig.update_layout(title="Wave height across the belt [m]")
col2.plotly_chart(fig)


col21, col22 = col2.columns(2)
col21.metric(
    "Transmitted wave", 
    f"{y[-1]:.3f} m", 
    delta=f"{(1-y[-1]/y[0])*100:.0f} %", 
    delta_arrow="down", 
    icon=":material/tsunami:"
)
col22.metric(
    "Transmitted energy", 
    f"{0.125*1000*9.81*(y[-1]**2):.0f} Jm⁻²", 
    delta=f"{0.125*1000*9.81*(y[0]**2-y[-1]**2):.0f} Jm⁻²", 
    delta_arrow="down", 
    icon=":material/flash_on:"
)


def slr_sensitivity(k, h, b):
    kh = k*h
    bh = b*h
    m = 2.0 - bh/(np.exp(bh)-1.0) - kh*(3*np.cosh(kh)*np.sinh(kh)+kh)/(np.sinh(kh)**2 + kh*np.tanh(kh))
    T = lambda f, kh, bh: f(3*kh-bh) + 3*f(kh-bh) + 3*f(-kh-bh) + f(-3*kh-bh)
    m +=  T(lambda x: np.exp(x), kh, bh) / T(lambda x: (np.exp(x)-1.0)/x, kh, bh)
    return m

m = slr_sensitivity(wave_number, water_depth, 1/root_height)
t = np.arange(0, 11, 1)
smin = np.array([0.28, 0.35, 0.44, 0.50, 0.56])
smed = np.array([0.41, 0.48, 0.58, 0.64, 0.72])
smax = np.array([0.60, 0.68, 0.80, 0.87, 0.97])
data = pd.DataFrame({
    "Scenarios":["SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"],
    "Lower":smin * m / water_depth,
    "Median":smed * m / water_depth,
    "Upper":smax * m / water_depth,
})
data["Difference"] = data["Upper"] - data["Lower"]
# data
fig = px.bar(
    data, 
    x="Scenarios", 
    y="Difference", 
    base="Lower",
    orientation='v',
    color="Upper",
    title="Efficiency subject to sea-level rise [%/year]",
    color_continuous_scale=px.colors.sequential.Sunset,
    opacity=0.9
)
fig.update_layout(coloraxis_showscale=False)
fig.update_layout(yaxis_title="")
col2.plotly_chart(fig)

if m > 0:
    col2.metric(
        "Sensitivity to sea-level rise [%/cm]", 
        f"{m/water_depth:.3f}", 
        border=True, 
        delta="Efficiency increase", 
        delta_color="green"
    )
else:
    col2.metric(
        "Sensitivity to sea-level rise [%/cm]", 
        f"{m/water_depth:.3f}", 
        border=True, 
        delta="Efficiency decrease", 
        delta_color="red", 
        delta_arrow="down"
    )