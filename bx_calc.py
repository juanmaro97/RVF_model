import matplotlib.pyplot as plt
from sklearn import metrics
from scipy.integrate import odeint
import numpy as np

from dewatering import dewatering_process

""" ------ CÁLCULO DEL BRIX DEL JUGO FILTRADO TOTAL EN EL CICLO ------ """

def calc_dens(Bx):
    dens = 1000*(0.956-0.005*Bx)     # Peña (2009)
    return dens

def calc_Bx(Bx_0,Vf,SS_agua,V_agua,r_torta):
    '''
    Bx_0: Brix inicial del lodo
    Vf: Volumen filtrado en la etapa de formación de la torta
    SS_agua: sólido solubles en el agua de lavado
    V_agua: volumen de agua empleada en el lavado. Se estima que la densidad es 1 y se iguala a kg.
    r_torta: retención de soluto en la torta después del proceso de lavado
    '''
    dens_filtrado = 1000*(0.956-0.005*Bx_0)
    # print("BX form dens: ",dens_filtrado)
    m_filtrado = Vf*dens_filtrado   
    SS_filtrado = (Bx_0/100)*m_filtrado

    Bx_jugo_final = ((SS_filtrado+SS_agua)/(m_filtrado+V_agua))*(1-r_torta)*100

    return Bx_jugo_final