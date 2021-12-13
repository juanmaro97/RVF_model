import matplotlib.pyplot as plt
from numpy.core.defchararray import array
from numpy.lib.function_base import append
from sklearn import metrics
from scipy.integrate import odeint
import numpy as np
from RVDF import RVDF, slurry_cake
import filtration
from pol_calc import pol_calc
import washing 
import dewatering
import bx_calc

""" ------ SCRIPT DE ARRANQUE DE SIMULACIÓN ------ """

####### ENTRADAS PARÁMETROS TAMBOR #########
P_diff = 20000                                          # Caída de presión [Pa]
rd = 1.5                                                # [m]
L = 4.3                                                 # [m]
filtration_angle = 75                                   # [°]
wsh_angle = 15                                          # [°]
dew_angle1 = 15                                         # [°]
dew_angle2 = 35                                         # [°]
w = 1.2                                                 # [rpm]
tf = RVDF.angle_to_time(filtration_angle,w)             # Tiempo de filtración [s]
Af = RVDF.drum_filter_area(rd,L,filtration_angle)       # Área de filtración [m2]
Rm = 0                                                  # Resistencia del medio filtrante
nivel_lodos = 0.95                                      # Nivel del tanque de lodos []

filtration_angle, dew_angle1 = RVDF.real_filtration_angle(filtration_angle,dew_angle1,nivel_lodos)     # estimar cambio en el ánglo de filtración por nivel de lodos

####### ENTRADAS PARÁMETROS TORTA Y LODOS ########
u = 400*10**-3                                 # Viscosidad de la suspensión [Pa*s]
c = 50                                         # Sólidos secos por unidad de volumen filtrado [kg/m3]
k = 2*10**-13                                  # Permeabilidad [m2]
s = 0                                          # índice de compresibilidad []
Bx_0 = 15                                      # Grados brix de lodo de entrada [kg/kg]
epsilon = 0.5                                  # Porosidad [espacio vacío/ espacio total]
solid_dens = 500                               # Densidad de sólidos[kg/m3]
filtrate_dens = bx_calc.calc_dens(Bx_0)        # Densidad del jugo filtrado[kg/m3]
surface_tension = 0.07                         # Tensión superficial del filtrado [N/m]
incompresible = True                           # DETERMINAR SI LA TORTA ES INCOMPRESIBLE O NO
alpha = slurry_cake.calc_alpha(incompresible,k,solid_dens,epsilon,0,P_diff,s)
wsh_Q = 0.0015                                 # tasa de lavado con agua [m3/s]

####### CREACIÓN DEL FILTRO Y LOS LODOS ###########
RVDF1 = RVDF(P_diff,rd,L,Af,tf,filtration_angle,wsh_angle,dew_angle1,dew_angle2,w,Rm,nivel_lodos)
LODOS_TORTA = slurry_cake(alpha,u,c,k,epsilon,solid_dens,filtrate_dens,surface_tension,incompresible,s)

""" --------- SIMULACIÓN ---------- """
ts = 0.1                  # tiempo de muestreo
t = np.arange(0, tf, ts)  # arreglo temporal

""" --------- ETAPA DE FILTRACIÓN Y FORMACIÓN DE LA TORTA ----------- """
v, Vf, q, l, Q_mean, W_cake = filtration.volume(RVDF1,LODOS_TORTA,ts,tf)                        

""" --------- ETAPA DE SECADO 1 ----------- """
S, M, irreduc_S, Vf_dew1, Vf_dew1_arr = dewatering.dewatering_process(RVDF1, LODOS_TORTA, thickness = float(l[-1]), 
                            dew_A = RVDF.drum_filter_area(rd,L,dew_angle1), 
                            dew_t = RVDF.angle_to_time(dew_angle1,w), zone = 1, r = 0, ts = ts)

""" --------- ETAPA DE LAVADO ----------- """
Vf_wsh, r = washing.water_wash(RVDF1, LODOS_TORTA, e = 0.8, wsh_Q = wsh_Q, thickness = float(l[-1]), Saturation=S)

""" --------- ETAPA DE SECADO 2 ----------- """
S, M, irreduc_S, Vf_dew2, Vf_dew2_arr = dewatering.dewatering_process(RVDF1, LODOS_TORTA, thickness = float(l[-1]), 
                            dew_A = RVDF.drum_filter_area(rd,L,dew_angle2), 
                            dew_t = RVDF.angle_to_time(dew_angle2,w), zone = 2, r = r, ts = ts)

""" --------- CÁLCULO DE BRIX DE SALIDA ----------- """                            
Bx_jugo_total = bx_calc.calc_Bx(Bx_0 = Bx_0, Vf = Vf, SS_agua = 0, V_agua = Vf_wsh, r_torta = r)

""" --------- CÁLCULO POL EN CACHAZA ------------- """
Pol_cachaza = pol_calc(msr = W_cake, filter_area = Af, wwr = wsh_Q)


print("####### RESUMEN ########")
print("Duración del ciclo de rotación [s]: ", 60/w)
print("Nivel de lodos [%]: ", nivel_lodos*100)
print("Ángulo real de filtración [°]: ", "{:.4f}".format(filtration_angle))
print("Volumen de filtrado en un ciclo [m3]: " "{:.4f}".format((Vf + Vf_dew1 + Vf_wsh + Vf_dew2)))
print("Tasa de filtración [m3/s]: " "{:.6f}".format((Vf + Vf_dew1 + Vf_wsh + Vf_dew2)/(60/w)))
print("Tasa de filtración [m3/h]: " "{:.6f}".format(3600*(Vf + Vf_dew1 + Vf_wsh + Vf_dew2)/(60/w)))
print("Brix de jugo filtrado en el ciclo [°]: ""{:.4f}".format(Bx_jugo_total))
print("Pol en cachaza [%]: ""{:.4f}".format(Pol_cachaza)+"\n")


""" ----- PLOT TOTAL PROCESS ---- """

### volumen arrays ###

Vf_dew1_arr = Vf_dew1_arr + v[-1]  # suma el último valor de formación con la filtración realizada en dewatering 1

arrays = append(v,Vf_dew1_arr)   # concatena vector de formación de torta con dewatering 1

watering_wash_vector = wsh_Q*(np.arange(0, RVDF.angle_to_time(wsh_angle,w), ts))+Vf_dew1_arr[-1]  # suma el último valor de dewatering 1 a la filtración en lavado

arrays = append(arrays,watering_wash_vector)    # concatena vectores hasta lavado

Vf_dew2_arr = Vf_dew2_arr + arrays[-1]     # suma último valor de lavado a la filtración de dewatering 2

arrays = append(arrays,Vf_dew2_arr)   # concatena vectores hasta dewatering 2

times_array = np.arange(0, ((arrays.shape[0]-ts)*ts), ts)   # crea arreglo temporal con el tiempo total de filtración

#### plot ###

plt.suptitle('Drum total cicle filtration')
plt.plot(times_array,arrays, 'r',label="Cane juice filtration [m3]")
plt.ylabel('Volume [m3]')
plt.xlabel('Time [s]')
plt.grid()
plt.legend()
plt.show()