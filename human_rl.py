import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import random 
import config
import os
from scipy.interpolate import CubicSpline
from mpl_toolkits.axes_grid1 import make_axes_locatable


def write_data(timestep , X, S, E):
    with open('plot.txt','a') as plotting_file:
        plotting_file.write(f"{timestep},{X},{S},{E}\n")
        plotting_file.close()

if __name__ == '__main__':

    if os.path.exists('plot.txt'):
        with open("plot.txt",'w') as f:
            f.truncate(0)
            f.close()
    else:
        with open('myfile.txt', 'w') as fp:
            pass
            fp.close()

    weibull_values = pd.read_csv('pdfcsv.csv')
    xvalues = weibull_values['x']
    y_values = weibull_values['pdf']

    # Create a cubic spline interpolation model
    cs = CubicSpline(xvalues, y_values)


    # Initial conditions
    X0 = config.X0  # CDW g/L 
    S0 = config.S0 # mol/L
    E0 = config.E0 # U/L 

    # Process conditions
    T = config.INIT_TEMPRATURE #'C 
    T_opt = config.OPTIMUM_TEMPERATURE #'C
    A = config.INIT_AGITATION # RPM

    # model parameters
    Ks = config.KS    # mol/L
    Yxs =  config.YXS  # CDW g/mol
    MuE_opt = config.MUE_OPT    # U/CDW g
    mu_max = config.MU_MAX   # /h
    del_t = config.DEL_T # hours ie. 36 seconds
    t_end = config.T_END
    total_sim_steps = int(t_end/del_t)
    tvec = [del_t + i * del_t for i in range(total_sim_steps)]
    ns = len(tvec)
    kl = config.KL # mol/L
    cell_death_timer = config.CELL_DEATH_TIMER
    cell_death_time = config.CELL_DEATH_TIME

    # Substrate addition calculations
    tank_capacity = config.TANK_CAPACITY # L
    substrate_in_tank_liters = config.SUBSTRATE_IN_TANK_LITERS # L
    max_substrate_limit_liters = config.MAX_SUBSTRATE_LIMIT_LITERS # L
    substrate_transfer_amount_liters = config.SUBSTRATE_TRANSFER_AMOUNT_LITERS # L 
    media_transfer_gap =  config.MEDIA_TRANSFER_GAP # Hours this is after 10 steps ie. 6 minutes
    media_transfer_step = int(media_transfer_gap/del_t)
    substrate_concentration = S0
    sub_in_tank_moles = S0 * substrate_in_tank_liters # mol

    # external media tank configurations
    ext_tank_substrate_conc = config.EXT_TANK_SUBSTRATE_CONC # mol/L
    substrate_transfer_moles = substrate_transfer_amount_liters * ext_tank_substrate_conc # mol 

    # X S E delE delX, muE
    D = np.zeros((ns+1, 6))
    D[0][0] = X0
    D[0][1] = substrate_concentration
    D[0][2] = E0


    for i in range(ns):

        if D[i][1] < 0.001 or D[i][1] > 0.045:
            MuX = 0
        else:
            MuX = 0.2

        # new cells that are generated
        dXdt = MuX * D[i][0]

        # SUbstrate consumption
        dSdt = -(1/Yxs)* dXdt 

        # Find change in cells 
        delX = dXdt * del_t
        # Find change in substrate
        delS = dSdt * del_t

        # Update cells
        D[i+1][0] = D[i][0] + delX

        # # Update substrate
        if substrate_in_tank_liters < max_substrate_limit_liters:
            if i%(media_transfer_step) == 0 and i != 0:
                print("------------- ENTER SUBSTRATE FLOW VOLUME : (0 - 100 ml) ------------")
                substrate_action = int(input())/1000
                if True:
                    # get substrate concentraion for current timestep
                    substrate_concentration = D[i][1]
                    sub_in_tank_moles = substrate_concentration * substrate_in_tank_liters
                    # Add substrate and calculate new concentraion
                    substrate_in_tank_liters = substrate_in_tank_liters + substrate_action # Liters of media
                    sub_in_tank_moles = sub_in_tank_moles + ext_tank_substrate_conc * substrate_action  # grams
                    substrate_concentration = sub_in_tank_moles/substrate_in_tank_liters # substrate grams/liter
                    D[i][1] = substrate_concentration

        # Check if substrate is less than or close to 0
        if D[i][1] + delS < 0.000001:
            D[i+1][1] = 0
        else:
            D[i+1][1] = D[i][1] + delS

        # Cells start dying if no substrate for more than 2 hours
        if D[i+1][1] == 0:
            cell_death_timer += 1
        else:
            cell_death_timer = 0
        
        # if cell_death_timer == 2 hours then cells start dying
        if int(cell_death_time/del_t) <= int(cell_death_timer):
            D[i+1][0] = D[i+1][0] - (D[i+1][0]*0.10)

        # Enzyme determination 
        sub_cell_ratio = D[i][1]/D[i][0]
        sub_cell_ratio = sub_cell_ratio * 1e6

        if sub_cell_ratio > 11000:
            weibull = 0
        else:
            weibull = cs(sub_cell_ratio)
        MuE = MuE_opt * weibull

        D[i][5] = MuE

        if dXdt == 0:
            MuE = 0

        # new enzyme from fresh cells
        delE = MuE * D[i][0] * del_t

        # Change in enzyme    
        D[i+1][3] = delE
        
        # Update enzyme variable
        D[i+1][2] = D[i][2] + delE

        # terminate if tank capacity is full and cells start dying
        if substrate_in_tank_liters >= max_substrate_limit_liters: #and D[i][1] <= 0:
            break

        write_data(i, D[i][0], D[i][1], D[i][2])
    
    X = D[:,0][0:i]
    S = D[:,1][0:i]
    E = D[:,2][0:i]
    MUE = D[:,5][0:i]
    S_C_R = S/X
    S_C_R = S_C_R * 1e6

    max_enzyme = round(max(E),3)
    sub_max = round(max(S),3)
    max_cells = round(max(X),3)

    '''
    ================================================================================================================================================
                                                        Plotting      
    ================================================================================================================================================
    '''
    
    print(" -----------------------------------------------------------------------")
    print(" ----------------------------  Results  --------------------------------")
    print(" MAX ENZYME ACTIVITY : ", max(D[:,2]))

    fig, ax1 = plt.subplots()
    fig.suptitle(f"Max Cells: {max_cells} CDWg/L, Max Substrate: {sub_max} mol/L, Max Enzyme Activity: {max_enzyme} U/L" )
    #fig.suptitle(f"Simulation: Fed Batch Mode")
    ax1.plot(tvec[0:i] ,E ,color="red", label= "Enzyme Activity U/L")

    ax2 = ax1.twinx()
    ax2.plot(tvec[0:i], S, color="orange", label="Substrate mol/L")

    ax3 = ax1.twinx()
    ax3.plot(tvec[0:i] , X ,color="blue", label="Cells CDW g/L")
    #ax3.spines['right'].set_position(('outward',60))
    ax3.spines['right'].set_position(('axes',1.15))

    ax1.set_ylabel("Enzyme Activity U/L", color="red")
    ax1.set_xlabel("Time (hours)")
    ax2.set_ylabel("Substrate mol/L ", color="orange")
    ax3.set_ylabel("Cells CDW g/L", color="blue")

    ax1.tick_params(axis='y',colors="red")
    ax2.tick_params(axis='y',colors="orange")
    ax3.tick_params(axis='y',colors="blue")

    ax2.spines['right'].set_color("orange")
    ax3.spines['right'].set_color("blue")
    ax3.spines['left'].set_color("red")
    #fig.legend()
    plt.show()

    plt.plot(tvec[:len(S_C_R)], S_C_R, label="Substrate to cell ratio (mol/g)")
    print("Initial substrate to cell ratio: ",S_C_R[0])
    plt.xlabel("Time H")
    plt.ylabel("Substrate to cell ratio (mol/g)")
    plt.legend()

    






