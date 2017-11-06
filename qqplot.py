import numpy as np, matplotlib.pyplot as plt
from matplotlib import rcParams

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 12
rcParams['legend.numpoints'] = 1

def read_q_q_file(file_name):
    f = open(file_name, 'r')
    lines = f.readlines()[1:]
    N_data = len(lines)

    # Initialise vectors
    quantiles = np.zeros(N_data)
    values = np.zeros(N_data)
    mn_values = np.zeros(N_data)
    # Load in data
    for i in range (0, N_data):
       line = lines[i].strip().split(" ")
       quantiles[i]=float(line[0])
       values[i]=float(line[1])
       mn_values[i]=float(line[2])
    f.close()
    return quantiles,values,mn_values

def make_q_q_plots(snv1,values1,mn_values1,snv2,values2,mn_values2):

   threshold = 0.005
   flag = 0
   min_length = 200
   range1 = np.ptp(values1)
   #print "Relief range: ", range1
   for i in range(0,len(snv1)):
        if (snv1[i] <= 0):
            frac_diff = abs((values1[i] - mn_values1[i]))/range1
            if (frac_diff < threshold):
                if (flag == 0):
                    flag = 1
                    count = 0
                    for j in range(1,min_length+1):
                        next_frac = abs((values1[i+j] - mn_values1[i+j]))/range1
                        if (next_frac < threshold):
                            count = count+1
                    if (count == min_length):
                        relief_thresh = snv1[i]
                        print "Relief threshold: ", values1[i]
                    else:
                        flag = 0

   flag = 0
   range2 = np.ptp(values2)
   print "Slope range: ", range2
   for i in range(0,len(snv2)):
       if (snv2[i] <= 0):
           frac_diff = abs((values2[i] - mn_values2[i]))/range2
           if (frac_diff < threshold):
                if (flag == 0):
                    flag = 1
                    count = 0
                    for j in range(1,min_length):
                        next_frac = abs((values2[i+j] - mn_values2[i+j]))/range2
                        if (next_frac < threshold):
                            count = count+1
                    if (count == min_length-1):
                        slope_thresh = snv2[i]
                        print "Slope threshold: ", values2[i]
                    else:
                        flag = 0
   print slope_thresh, relief_thresh

   plt.figure(1, facecolor='White',figsize=[10,5])
   ax1 = plt.subplot(1,2,1)
   ax1.plot(snv1,values1,linewidth=2,color="blue",label="Real data")
   ax1.plot(snv1,mn_values1,"--",linewidth=2,color="red",label="Normal distribution")
   ax1.axvline(x=relief_thresh,linestyle='--',linewidth=1.5,color='black')
   xmin,xmax = ax1.get_xlim()
   ax1.axvspan(xmin, relief_thresh, alpha = 0.2, color='blue')
   ax1.legend(loc = 2)
   ax1.set_xlabel("Standard Normal Variate", fontsize=rcParams['font.size']+2)
   ax1.set_ylabel("Channel relief ($m$)", fontsize=rcParams['font.size']+2)
   ax1.set_xlim(xmin,xmax)
   ax1.grid(True)


   ax2 = plt.subplot(1,2,2)
   ax2.plot(snv2,values2,linewidth=2,color="blue",label="Real data")
   ax2.plot(snv2,mn_values2,"--",linewidth=2,color="red",label="Normal distribution")
   ax2.axvline(x=slope_thresh,linestyle='--',linewidth=1.5,color='black')
   xmin2,xmax2 = ax2.get_xlim()
   ax2.axvspan(xmin2, slope_thresh, alpha = 0.2, color='blue')
   #ax2.legend(loc = 2)
   ax2.set_xlabel("Standard Normal Variate", fontsize=rcParams['font.size']+2)
   ax2.set_ylabel("Gradient", fontsize=rcParams['font.size']+2)
   ax2.set_xlim(xmin2,xmax2)
   ax2.grid(True)
   plt.tight_layout()

if __name__ == "__main__":
    DataDirectory="/media/fionaclubb/terrace_lidar/DEMs_for_analysis/Miss_test"
    DEM_name = 'Miss_test'
    relief_file=DataDirectory+DEM_name+"_qq_relief.txt"
    slope_file=DataDirectory+DEM_name+"_qq_slope.txt"
    OutputName = DataDirectory+DEM_name+"_qq_plots"
    dot = "."
    OutputFormat = "png"
    x,y1,y2=read_q_q_file(relief_file)
    slope_x,slope_y1,slope_y2 = read_q_q_file(slope_file)
    #plt.show()
    make_q_q_plots(x,y1,y2,slope_x,slope_y1,slope_y2)
    #plt.show()
    plt.savefig((OutputName+dot+OutputFormat), format=OutputFormat)
    plt.clf()
