import streamlit as st
import pandas as pd
import numpy as np
from scipy.optimize import fsolve

# Set the title and favicon that appear in the Browser's tab bar.
st.set_page_config(
    page_title='Mangrove Wave Attenuation Tool',
    page_icon=':ocean:', # This is an emoji shortcode. Could be a URL too.
)

'''
# :ocean: Mangrove Wave Attenuation Tool

Mangrove Wave Attenuation Tool is a simple interactive calculator that 
estimates wave attenuation through mangrove forests using wave input 
parameters. It helps visualize how mangroves reduce wave energy and 
support coastal resilience. The wave attenuation model is based on the 
work by Pang and Tay (https://arxiv.org/abs/2606.11653).
'''

water_depth = st.number_input("Water Depth h [m]", value=1.0)
wave_period = st.number_input("Wave Period T [s]", value=1.0)
wave_height = st.number_input("Wave Height H [m]", value=1.0)
gravity = st.number_input("Gravity g [m/s²]", value=9.81, disabled=True)
wave_angular_frequency = 2*np.pi/wave_period
wave_number = fsolve(lambda k: gravity*k*np.tanh(k*water_depth) - wave_angular_frequency**2, 1.0)
wave_number = st.number_input("Wave Number k [1/m]", value=wave_number[0], disabled=True)
drag_coefficient = st.number_input("Drag Coefficient Cd [-]", value=1.0)
tree_density = st.number_input("Tree Density [trees/m²]", value=0.5)
mangrove_belt_width = st.number_input("Mangrove Belt Width [m]", value=100.0)
mangrove_species = st.menu_button("Mangrove Root Parameters Preset", options=["Rhizophora", "Sonneratia"])

tree_beta = 0.0
tree_a0 = 0.0
if mangrove_species == "Rhizophora":
    tree_beta = 1.7
    tree_a0 = 7.0
elif mangrove_species == "Sonneratia":
    tree_beta = 6.0
    tree_a0 = 1.0
tree_beta = st.number_input("Root Shape Parameters [1/m]", value=tree_beta)
tree_a0 = st.number_input("Root Number * Root Diameter [m]", value=tree_a0)

# normalised KD
def KD(k, h, a0, b):
    kh = k*h
    bh = b*h
    abar = a0/bh * (1.0 - np.exp(-bh))
    f = lambda x: (np.exp(x)-1.0)/x if np.abs(x)>0.001 else 0.5*x+1.0
    phi = a0*h/(8*np.sinh(kh)**3) * (f(3*kh-bh) + 3*f(kh-bh) + 3*f(-kh-bh) + f(-3*kh-bh))
    val = 2*k**2*np.tanh(kh)**2 / (kh + np.tanh(kh) - kh*np.tanh(kh)**2) * phi
    return val

wave_decay_coefficient = st.number_input("Wave decay coefficient [1/m²]", value=KD(wave_number, water_depth, tree_a0, tree_beta+0.0001))

x = np.linspace(0, mangrove_belt_width, 100)
y = wave_height/(1+wave_decay_coefficient*wave_height*x)
df = pd.DataFrame({'Distance [m]': x, 'Wave Height [m]': y})

st.line_chart(df, x="Distance [m]", y="Wave Height [m]")
