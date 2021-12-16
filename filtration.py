import matplotlib.pyplot as plt
import math
from sklearn import metrics
from scipy.integrate import odeint

import numpy as np

""" ------ CÁLCULO DE LA ETAPA DE FILTRACIÓN Y FORMACIÓN DE LA TORTA ------ """
### Bioseparaciones (2011), Tejeda, et. al
### Operaciones unitarias en ingeniería química (2007), McCabe, et al
### Solid/Liquid separations Principles of Industrial Filtration (2005) Wakeman, Tarleton

def f(V,t,P_diff,A,u,alpha,c,Rm):
    """this is the rhs of the ODE to integrate, i.e. dV/dt=f(V,t)"""
    '''
    V: Volumen filtrado [m3]
    t: Tiempo [s]
    P: Diferencial de presión [Pa]
    A: Área de filtración [m2]
    u: Viscosidad [Pa*s]
    alpha: Resistencia específica de la torta [m/kg]
    c: Sólidos secos por unidad de volumen filtrado [kg/m3]
    Rm: Resistencia del medio filtrante [1/m]
    '''

    return (P_diff*A)/(u*((alpha*c*(V/A))+Rm))

def calc_Q(rd,L,phi,omega,P_diff,u,alpha,c,s):
    '''
    Q: Flujo volumétrico [m3/t]
    rd: Radio del tambor [m]
    L: Largo del tambor [m]
    phi: Ángulo de filtración del tambor [°]
    omega: Velocidad de rotación del tambor [rpm]
    P_diff: Diff de presión del tambor [Pa]
    u: Viscosidad [Pa*s]
    alpha: Resistencia específica de la torta [m/kg]
    c: Sólidos secos por unidad de volumen filtrado [kg/m3]
    '''
    omega = omega/60 # revolutions per min to revolutions per second
    phi = np.deg2rad(phi)

    Q = rd*L*((4*math.pi*phi*omega*(P_diff**(1-s)))/(u*alpha*c))**0.5

    return Q

def volume(filtro,slurry,ts,tf):

    P_diff = filtro.P_diff
    A = filtro.Af
    Rm = filtro.Rm
    rd = filtro.rd
    L = filtro.L
    w = filtro.w
    filtration_angle = filtro.filtration_angle
    u = slurry.u
    alpha = slurry.alpha
    c = slurry.c
    s = slurry.s
    k = slurry.k
    
    v0 = 0.0001
    t = np.arange(0, tf, ts)  # values of t for
                          # which we require
                          # the solution y(t)

    v = odeint(f, v0, t, args = (P_diff,A,u,alpha,c,Rm))  # actual computation of v(t)
    
    Vf = float(v[-1])   # Volumen filtrado [m3] (valor final de v)

    q = (A*A*P_diff)/(u*alpha*c*v)     # q instantánea [m3/s]

    # l = (alpha*c*v*k)/A

    l = (k*A*P_diff)/(u*q)      # cálculo de el espesor de la torta

    Q_mean = calc_Q(rd,L,filtration_angle,w,P_diff,u,alpha,c,s)    # velocidad de filtración media [m3/s]

    W_cake = c*calc_Q(rd,L,filtration_angle,w,P_diff,u,alpha,c,s)   # velocidad de formación de la torta [kg/s]

    #### PLOTTING #### 

    print("\n------ SIMULACIÓN DINÁMICA: FILTRACIÓN EN UN TAMBOR ROTATORIO DE VACÍO ------\n")
    print("###### FILTRATION ########")
    print("Tiempo de filtración [s]: ","{:.2f}".format(tf))
    print("Q media [m3/s]: ","{:.6f}".format(Q_mean))
    print("Espesor torta [mm] es: ","{:.6f}".format(float(l[-1]*1000)))
    print("V filtrado [m3]: ","{:.6f}".format(metrics.auc(t,q)))
    print("Producción de sólidos [kg/s]: ","{:.6f}".format(W_cake)+"\n")

    fig, axarr = plt.subplots(3, 1)
    fig.canvas.set_window_title('RVF') 
    fig.suptitle('Filtración en tambor rotatorio al vacío')

    axarr[0].plot(t, v, 'b', label='V [m3]')
    axarr[0].legend(loc='best')
    axarr[0].grid()
    axarr[0].set(ylabel = 'V [m3]')
    axarr[1].set_title('AUC Q: ' +str("{:.5f}".format(metrics.auc(t,q))+ " [m3]"))
    axarr[1].plot(t, q, 'r', label='Q [m3/s]')
    axarr[1].legend(loc='best')
    axarr[1].grid()
    axarr[1].set(ylabel = 'Q [m3/s]')
    axarr[2].plot(t, l, 'g', label='L [m]')
    axarr[2].legend(loc='best')
    axarr[2].grid()
    axarr[2].set(xlabel = 'Time [s]', ylabel = 'L [m]')

    plt.show()

    return v, Vf, q, l, Q_mean, W_cake

