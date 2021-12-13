import matplotlib.pyplot as plt
import math
from numpy.lib.financial import irr
from sklearn import metrics
from scipy.integrate import odeint
import numpy as np

""" ------ SIMULACIÓN DINÁMICA: FILTRACIÓN EN UN TAMBOR ROTATORIO DE VACÍO ------ """
### Bioseparaciones (2011), Tejeda, et. al
### Operaciones unitarias en ingeniería química (2007), McCabe, et al
### Solid/Liquid separations Principles of Industrial Filtration (2005) Wakeman, Tarleton

class RVDF:
    def __init__(self,P_diff,rd,L,Af,tf, filtration_angle,wsh_angle,dew_angle1,dew_angle2,w,Rm,nivel):
        '''
        PARÁMETROS DEL FILTRO DE TAMBOR ROTATORIO AL VACÍO

        P_diff: Presión del vacío en el tambor [Pa]
        rd: Radio del tambor [m]
        L: Largo del tambor [m]
        Af: Área de filtración [m2]
        tf: Tiempo de filtración
        filtration_angle: Ángulo de filtración del tambor [°]
        wsh_angle: Ángulo de lavado del tambor [°]
        dew_angle: Ángulo de secado del tambor [°]
        w: Velocidad de rotación [rpm]
        Rm: Resistencia del medio filtrante
        nivel: Nivel del tanque de lodos []
        '''
        self.P_diff = P_diff    # Caída de presión [Pa]
        self.rd = rd              # [m]
        self.L = L              # [m]
        self.filtration_angle = filtration_angle   # [°]
        self.wsh_angle = wsh_angle   # [°]
        self.dew_angle1 = dew_angle1   # [°]
        self.dew_angle2 = dew_angle2   # [°]
        self.w = w             # [rpm]
        self.Af = Af            # Área de filtración [m2]
        self.tf = tf            # Tiempo de filtración [s]
        self.Rm = Rm      # Resistencia del medio filtrante
        self.nivel = nivel      # nivel del tanque []

    def drum_filter_area(rd,L,filtration_angle):
        '''
        CÁLCULO DEL ÁREA DE FILTRACIÓN DEL TAMBOR
        rd: Radio del tambor [m]
        L: Largo del tambor [m]
        filtration_angle: Ángulo de filtración del tambor [°]
        '''
        A = 2*math.pi*rd*L*(filtration_angle/360)
        
        return A

    def angle_to_time(phi,omega):
        '''
        CÁLCULO DEL ÁNGULO DE FILTRACIÓN A TIEMPO DE FILTRACIÓN
        phi: Ángulo [°]
        omega: Velocidad de rotación del tambor [rpm]
        '''
        phi = np.deg2rad(phi)
        omega = omega/60 # revolutions per min to revolutions per second
        tf = phi/(2*math.pi*omega)   # [s]

        return tf

    def real_filtration_angle(filtration_angle, dew1_angle, nivel):
        theta_perdido =  np.rad2deg(np.arcsin(1-nivel))
        filtration_angle = filtration_angle - theta_perdido
        dew1_angle = dew1_angle + theta_perdido

        return filtration_angle, dew1_angle

class slurry_cake:

    """ PARÁMETROS DE LOS LODOS Y LA TORTA """

    def __init__(self,alpha,u,c,k,epsilon,solid_dens,filtrate_dens,surface_tension,incompressible,s):
        self.u = u     # Viscosidad de la suspensión [Pa*s]
        self.c = c        # Sólidos secos por unidad de volumen filtrado [kg/m3]
        self.k = k      # Permeabilidad de la torta[m2]
        self.s = s         # índice de compresibilidad []
        self.epsilon = epsilon      # Porosidad de la torta []
        self.solid_dens = solid_dens         # Densidad de sólidos [kg/m3]
        self.filtrate_dens = filtrate_dens        # Densidad del jugo filtrado[kg/m3]
        self.surface_tension = surface_tension      # Tensión superficial del filtrado [N/m]
        self.incompressible = incompressible        # True / False
        self.alpha = alpha                        # Resistencia específica de la torta [m/kg]
                  

    def calc_alpha(incompressible,k,solid_dens,epsilon,alpha_prima,P_diff,s):

        """ CÁLCULO DE LA RESISTENCIA DINÁMICA DE LA TORTA """

        if incompressible:
            s = 0
            alpha = 1/(k*solid_dens*(1-epsilon)) # [m/kg]
        else:
            alpha = alpha_prima*(P_diff**s)    # [m/kg]

        return alpha