import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

def get_results(filepath):
    datafile =open(filepath,"r")
    var_names = []
    values = []
    var_data =datafile.read().split("VARIABLES=")[1].split("\n\n")
    var_names =var_data[0].replace('"','').strip().split(" ")
    values = np.genfromtxt(filepath,skip_header=24)
    return [var_names,values]

#basepath = "E:\\sim"
basepath ="\\\\sshfs\\s2133342@tfe2.ctw.utwente.nl\\thesis"

#filepath="C:\\Users\\Borger\\Documents\\uni\\Thesis\\sim\\"
#filepath="E:\\sim\\temp\\1eq\\"
filepath=basepath+"\\eppler_sim\\Re2e5_Ma06\\"
#filepath= basepath + "\\NACA_sim\\"


matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['lines.markersize'] = 8

#group constant residuals
#group constant values
index_residuals=[1, 2, 3, 4]

j=0
prev_path = ""
for path, folders, files in os.walk(filepath):
    # List contain of folder
    for folder_name in folders:
        if folder_name == "results":
            if os.path.exists(path+"\\"+folder_name+"\\result.monitor.pval.dat"):
                [var_names,data] = get_results(path+"\\"+folder_name+"\\result.monitor.pval.dat")
            
            elif os.path.exists(path+"\\"+folder_name+"\\result.monitor.tmp.dat"):
                [var_names,data] = get_results(path+"\\"+folder_name+"\\result.monitor.tmp.dat")

            else:
                continue
        else:
            continue
        path_folders = path.split("\\")
        alpha = path_folders[-1].split("_")[-1]
        if prev_path != path_folders[-2]:
            if prev_path != "":
                axs1.legend(loc = "upper left", bbox_to_anchor=[1.2,1.1])
                fig1.savefig(basepath+'\\postprocessing\\convergence\\'+prev_path+".png", bbox_inches='tight')
            prev_path = path_folders[-2]
            fig1 = plt.figure()
            axs1=fig1.add_subplot(111)

            axs1.set_xlabel("Iter")
            axs1.set_ylabel("R_rho")
            axs1.set_yscale("log")
            axs1.set_title('Density residuals '+ path_folders[-2])
        
        No_of_points = list(range(0, len(data), 1000))
        fig2 = plt.figure()
        axs2=fig2.add_subplot(111)

        axs2.set_xlabel("Iter")
        axs2.set_ylabel("Residuals")
        axs2.set_yscale("log")
        axs2.set_title(path_folders[-2]+" "+path_folders[-1])
        #plot data
        axs1.plot(data[No_of_points,0],data[No_of_points,1], label=path_folders[-1])
        for i in index_residuals:
            axs2.plot(data[No_of_points,0],data[No_of_points,i], label=var_names[i])
        axs2.legend(loc='upper right')

        
            

axs1.legend(loc = "upper left", bbox_to_anchor=[1.2,1.1])
fig1.savefig(basepath+'\\postprocessing\\convergence\\'+prev_path+".png", bbox_inches='tight')