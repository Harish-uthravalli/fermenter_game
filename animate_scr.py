import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import config


style.use('fivethirtyeight')

fig, ax1 =  plt.subplots()

def animate(i):
    graph_data = open('plot.txt','r').read()
    lines = graph_data.split('\n')
    T = []
    SCR = []
    for line in lines:
        if len(line) > 1:
            t, x , s ,_ = line.split(',')
            ratio = float(s)/float(x)
            SCR.append(ratio * 1e6)
            T.append(float(t)/100)

    ax1.clear()
    
    ax1.plot(T , SCR ,color="red", label= "Temperature")
    ax1.set_ylabel("Substrate to Cell ratio", color="red")
    ax1.set_xlabel("Time (hours)")
    ax1.tick_params(axis='y',colors="red")
    ax1.axhline(y=config.OPT_SUB_CELL_RATIO, color='blue', linestyle='--', label="Threshold")
    
    

ani = animation.FuncAnimation(fig, animate, interval=1)
plt.show()