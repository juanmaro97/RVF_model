import matplotlib.pyplot as plt
from sklearn import metrics
from scipy.integrate import odeint
import numpy as np

""" ------ CÁLCULO DE LA POL EN CACHAZA ------ """
#### The balance between capacity and performance
#### of rotary mud filters (1997), Wright, Steggles, Steindl

def pol_calc(msr,filter_area,wwr):
    """
    msr: mud solid rate - tasa de sólidos de los lodos (kg/s, ton/h)
    wwr: wash water rate - tasa de agua de lavado (kg/s, ton/h)
    filter_area: área de filtración del filtro (m2)
    MSL: mud solids loading (t/h/m2)
    WWL: wash water loading (t/h/m2)

    """
    msr = (msr*3600)/1000   # kg/s -> ton/h
    wwr = (wwr*3600)/1000   # kg/s -> ton/h
    MSL = (100*msr)/(filter_area)
    WWL = (100*wwr)/(filter_area)
    # print("msr: ",msr)
    # print("wwr: ",wwr)
    # print("WWL: ",WWL)
    # print("MSL: ",MSL)
    # print("area: ",filter_area)
    pol_ms = 1.265*MSL**(2.5)+30.78*WWL**(-0.25)-12.43

    return pol_ms
