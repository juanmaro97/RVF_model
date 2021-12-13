import matplotlib.pyplot as plt
import math
from sklearn import metrics
from scipy.integrate import odeint
import numpy as np

""" ------ CÁLCULO DE LA ETAPA DE SECADO O DESHIDRATACIÓN DE LA TORTA ------ """
######### CAP 5: Solids/Liquids separations Principles of Industrial Filtration (2005) Wakeman, Tarleton ########

def dewatering_process(filtro,slurry,thickness,dew_A,dew_t,zone,r,ts):
    
    P_diff = filtro.P_diff                      # Presión diferencial del filtro
    epsilon = slurry.epsilon                    # Porosidad 
    alpha = slurry.alpha                        # Resistencia dinámica de la torta
    solid_dens = slurry.solid_dens              # Densidad de los sólidos
    filtrate_dens = slurry.filtrate_dens        # Densidad del filtrado
    surface_tension = slurry.surface_tension    # Tensión superficial
    viscosity = slurry.u                        # Viscosidad de lodos
    k = slurry.k                                # Permeabilidad de la torta
    

    x = 13.4*((1-epsilon)/(alpha*solid_dens*epsilon**3))**0.5        # Tamaño de partículas promedio [m]                   
    Ncap = ((epsilon**3)*x*x*(filtrate_dens*9.81*thickness+P_diff))/(((1-epsilon)**2)*thickness*surface_tension)     # Número capilar []
    irreduc_sat = 0.155*(1+0.031*Ncap**-0.49)   # Saturación irreducible []
    pb = (4.6*(1-epsilon)*surface_tension)/(epsilon*x)    # Presión normalizada 
    p_dimesionless = P_diff/pb                  # Presión adimensional [N/m2]
    dew_time = np.arange(0, dew_t, ts)        # Tiempo de secado [s]
    dew_t_theta = (k*pb*dew_time)/(viscosity*epsilon*(1-irreduc_sat)*thickness**2)      # Tiempo adimensional []
    
    if dew_t_theta[-1]*p_dimesionless <= 1.915:
        SR = 1/(1+1.08*(dew_t_theta*p_dimesionless)**0.88)    # Reducción de saturación []
    else:
        SR = 1/(1+1.46*(dew_t_theta*p_dimesionless)**0.48)    # Reducción de saturación []
    S = SR*(1-irreduc_sat)+irreduc_sat     # Saturación []
    M = (S*epsilon*filtrate_dens*100)/((1-epsilon)*solid_dens)
    M = (M/(100+M))*100     # Humedad [%]

    V_filtrate_dew = (1-float(S[-1]))*epsilon*dew_A*thickness    # Cantidad de filtrado en etapa de secado [m3]
    V_filtate_dew_arr = (1-S)*epsilon*dew_A*thickness

    plt.suptitle('Dewatering process '+str(zone))
    plt.plot(dew_time,S*100, label="Saturation")
    # plt.plot(dew_time,V_filtate_dew_time, label="Saturation")
    plt.plot(dew_time,M, label="Moisture")
    plt.ylabel('Percentage [%]')
    plt.xlabel('Time [s]')
    plt.grid()
    plt.legend()
    plt.show()

    if zone == 1:
        print("###### DEWATERING 1 ########")
        print("Filtrado en dewatering 1 [m3]: ","{:.6f}".format(V_filtrate_dew)+"\n")
        
    if zone == 2:
        print("###### DEWATERING 2 ########")
        print("Filtrado en dewatering 2 [m3]: ","{:.6f}".format(V_filtrate_dew))
        print("Humedad en cachaza [%]: ","{:.4f}".format(float(M[-1])))
        print("Saturación irreducible de la torta []: ","{:.4f}".format(irreduc_sat))
        print("Saturación de final []: ","{:.4f}".format(float(S[-1])))
        print("Retención final de soluto en cachaza []: ","{:.4f}".format(float(S[-1])*r)+"\n")

    return float(S[-1]), float(M[-1]), irreduc_sat,  V_filtrate_dew, V_filtate_dew_arr # retornar valores finales de saturación y humedad