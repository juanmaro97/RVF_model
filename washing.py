import matplotlib.pyplot as plt
import math
from sklearn import metrics
from scipy.integrate import odeint
import numpy as np
from RVDF import RVDF, slurry_cake

""" ------ CÁLCULO DE ETAPA DE LAVADO ------ """
### Bioseparaciones (2011), Tejeda, et. al
### Operaciones unitarias en ingeniería química (2007), McCabe, et al
### Solid/Liquid separations Principles of Industrial Filtration (2005) Wakeman, Tarleton

def water_wash(filtro,slurry,e,wsh_Q,thickness,Saturation):
    wsh_angle = filtro.wsh_angle    # ángulo de lavado [°]
    w = filtro.w                    # velocidad de rotación [rpm]
    rd = filtro.rd                  # radio del tambor [m]
    L = filtro.L                    # Longitud del tambor [m]
    epsilon = slurry.epsilon        # Porosidad de la torta []

    water_visc = 0.8 # 0.8 cP Viscosidad del agua de lavado a 30 °C
    wsh_time = RVDF.angle_to_time(wsh_angle,w)          # tiempo de lavado [s]
    wsh_A = RVDF.drum_filter_area(rd,L,wsh_angle)       # área de lavado [m2]

    wsh_ratio = (wsh_Q*wsh_time)/(epsilon*wsh_A*thickness*Saturation)    # Tasa de lavado
    r = (1-e)**wsh_ratio                                                 # Retención de soluto en la torta
    Vf_wsh = wsh_Q*wsh_time - epsilon*wsh_A*thickness*(1-Saturation)     # Cantidad de volumen filtrado en lavado

    print("###### WASHING #######")
    print("Eficiencia del lavado: ","{:.4f}".format(e))
    print("Retención de soluto en cachaza: ","{:.4f}".format(r))
    print("Tiempo de lavado [s]: ","{:.4f}".format(wsh_time))
    print("Cantidad de lavado empleado [m3]: ","{:.6f}".format(wsh_Q*wsh_time))
    print("Wash ratio: ","{:.4f}".format(wsh_ratio)+"\n")

    return Vf_wsh, r
