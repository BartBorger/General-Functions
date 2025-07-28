import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import re
import pyvista as pv

def plot_cp_cf(multiblock):           
    fig1= plt.figure()
    axs1=fig1.add_subplot(111)

    axs1.set_xlabel("X")
    axs1.set_ylabel("Cp")
    axs1.set_title(path.split("\\")[-2]+" " +path.split("\\")[-1])
    axs1r = axs1.twinx()
    axs1r.set_ylabel("Cf")

    # Loop through blocks
    for i, block in enumerate(multiblock):
        if block is not None:
            
            # Extract Cp and Cf from point data (or cell data as fallback)
            point_data = block.point_data
            cp = point_data.get("cp")  if "cp" in point_data else None
            cf = point_data.get("cf")  if "cf" in point_data else None
            # Extract X-coordinates
            x_coords = block.points[:, 0]  
            label_plot = [None,None]
            if multiblock.get_block_name(i) == "lower" or multiblock.get_block_name(i) == "lower_surface":
                x_coords[0] = 0.0
                label_plot = ["Cp","Cf"]
                axis_max =max(cp)
                axis_min =min(cp)
            N_points = int(len(x_coords)/2)
            cp = cp[0:N_points]
            cf = cf[0:N_points]
            x_coords = x_coords[0:N_points]
            axs1.plot(x_coords, cp, label=label_plot[0], color='black')
            axs1r.plot(x_coords, cf, label=label_plot[1], color='red')
    [exp_lower,exp_upper] = get_exp_results(path_folders)
    
    if len(exp_lower)>0:
        axs1.scatter(exp_lower[:,0], exp_lower[:,1], label="Exp", color='black', marker='s')
    if len(exp_upper)>0:
        axs1.scatter(exp_upper[:, 0], exp_upper[:, 1], color='black', marker='s')
 
    
    lines_1, labels_1 = axs1.get_legend_handles_labels()
    lines_2, labels_2 = axs1r.get_legend_handles_labels()
    axs1.set_ylim(max(axis_max,max(cp))+0.2,min(axis_min,min(cp))-0.2)
    axs1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left', bbox_to_anchor=[1.2,1.1])
    save_folder = basepath+'\\postprocessing\\Cp_Cf\\'+path_folders[-2]+"\\"
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    fig1.savefig(save_folder+path_folders[-1]+".png", bbox_inches='tight')
def plot_comparison(multiblock,axs1,axs1r):         
    # Loop through blocks
    linestyles = ['-','-', '--', '--', '-.','-.', ':', ':']
    for i, block in enumerate(multiblock):
        if block is not None:
            
            # Extract Cp and Cf from point data (or cell data as fallback)
            point_data = block.point_data
            cp = point_data.get("cp")  if "cp" in point_data else None
            cf = point_data.get("cf")  if "cf" in point_data else None
            # Extract X-coordinates
            x_coords = block.points[:, 0]  
            label_plot = [None,None]
            if multiblock.get_block_name(i) == "lower" or multiblock.get_block_name(i) == "lower_surface":
                x_coords[0] = 0.0
                label_plot = ["Cp "+path_folders[-2],"Cf "+path_folders[-2]]
            N_points = int(len(x_coords)/2)
            cp = cp[0:N_points]
            cf = cf[0:N_points]
            x_coords = x_coords[0:N_points]
            axs1.plot(x_coords, cp, label=label_plot[0], color='black', linestyle=linestyles[len(axs1.lines) % len(linestyles)])
            axs1r.plot(x_coords, cf, label=label_plot[1], color='red',  linestyle=linestyles[len(axs1r.lines) % len(linestyles)])
    
def get_multiblock(folderpath):

    # Regex pattern to match files
    pattern = re.compile(r"result\.surface\.pval\.(\d+)\.vtm")
    # Get all filenames in the folder
    files = os.listdir(folderpath)

    # Filter and find the one with the highest number
    matched_files = [f for f in files if pattern.match(f)]
    if not matched_files:
        print("No matching files found in", folderpath)
        return None
    max_file = max(matched_files, key=lambda x: int(pattern.match(x).group(1)))
    
    # Load the .vtm file
    multiblock = pv.read(folderpath + "\\" + max_file)
    return multiblock
            
def get_exp_results(path_folders):
    if path_folders ==[]:
        return [[],[]]
    alpha = path_folders[-1].split("_")[-1]
    match path_folders[5]:
        case "eppler_sim":
            exp_path ="E:\\Sim\\eppler_sim\\Postprocessing\\Exp\\McGhee1989\\Re2e5\\McGhee1989_Fig17_Re2e5_alpha" + alpha
            # Load experimental data
            lower_data = []
            upper_data = []
            if os.path.exists(exp_path+".dat"):
                lower_data = np.genfromtxt(exp_path+".dat",skip_header=1)
            
            if os.path.exists(exp_path+"_upper.dat"):
                upper_data = np.genfromtxt(exp_path+"_upper.dat",skip_header=1)
            
            return [lower_data, upper_data]
        case "nlf_sim":
            exp_path ="E:\\Sim\\nlf_sim\\Postprocessing\\Exp\\nlf0416_Re4.0e6_Ma0.1_cp.dat"
            # Load experimental data
            lower_data = []
            upper_data = []
            datablocks = open(exp_path, 'r').read().split("zone t = ")
            for i in range(0,len(datablocks)):
                if datablocks[i].find("a "+  alpha) == -1:
                    continue

                lower_data= np.genfromtxt(datablocks[i].split("\n")[1:])[:,[0,3]]
                upper_data= np.genfromtxt(datablocks[i+1].split("\n")[1:])[:,[0,3]]
                break
            
            return [lower_data, upper_data]
        case "NACA_sim":
            exp_path ="E:\\Sim\\NACA_sim\\Postprocessing\\Exp\\cp\\cp_exp_inboard_nosweep.txt"
            exp_path2 ="E:\\Sim\\NACA_sim\\Postprocessing\\Exp\\cp\\cp_exp_outboard_nosweep.txt"
            if alpha.find(".") == -1:
                alpha = alpha + ".0"
            # Load experimental data
            lower_data = []
            upper_data = []
            datablocks = open(exp_path, 'r').read().split("zone t = ")
            datablocks2 = open(exp_path2, 'r').read().split("zone t = ")
            for i in range(0,len(datablocks)):
                if datablocks[i].find("alpha "+  alpha) == -1:
                    continue

                lower_data1= np.genfromtxt(datablocks[i].split("\n")[1:])
                lower_data2= np.genfromtxt(datablocks2[i].split("\n")[1:])
                lower_data = (lower_data1+lower_data2)/2
                break
            
            return [lower_data, upper_data]
    return [[],[]]
    

matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['lines.markersize'] = 8

#basepath = "E:\\sim\\"
basepath = "\\\\sshfs\\s2133342@tfe2.ctw.utwente.nl\\thesis\\"

#filepath=basepath+"eppler_sim\\Re2e5_Ma06\\"
#filepath=basepath+"nlf_sim\\1eq_prep_released\\"
filepath=basepath+"NACA_sim\\"
    
#compare between cases
cases = [filepath + name+"\\" for name in os.listdir(filepath) ]
angles = []
for case in cases:
    [angles.append(angle) for angle in os.listdir(case)]

angles=list(set(angles))

for angle in angles:
    if angle.find("eppler") == -1 and angle.find("NLF") == -1 and angle.find("naca") == -1:
        continue
    fig1 = plt.figure()
    axs1 = fig1.add_subplot(111)
    axs1.set_xlabel("X")
    axs1.set_ylabel("Cp")
    axs1.set_title("Cp comparison at angle " + angle)
    fig2 = plt.figure()
    axs2 = fig2.add_subplot(111)
    axs2.set_xlabel("X")
    axs2.set_ylabel("Cf")
    axs2.set_title("Cf comparison at angle " + angle)
    path_folders = []
    for case in cases:
        if not os.path.exists(case + angle + "\\results"):
            continue
        multiblock = get_multiblock(case + angle + "\\results")
        if multiblock is None:
            continue
        path = case + angle 
        path_folders =path.split("\\")
        plot_comparison(multiblock,axs1,axs2)
    
    [exp_lower,exp_upper] = get_exp_results(path_folders)
    
    if len(exp_lower)>0:
        axs1.scatter(exp_lower[:,0], exp_lower[:,1], label="Exp", color='black', marker='s')
    if len(exp_upper)>0:
        axs1.scatter(exp_upper[:, 0], exp_upper[:, 1], color='black', marker='s')
    axs1.set_ylim(axs1.get_ylim()[1],axs1.get_ylim()[0])
    axs1.legend(loc='upper left', bbox_to_anchor=[1.2,1.1])
    axs2.legend(loc='upper left', bbox_to_anchor=[1.2,1.1])

    save_folder = basepath + '\\postprocessing\\Cp_comparison\\'
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    if not len(axs1.lines) == 0:
        fig1.savefig(save_folder + angle + '.png', bbox_inches='tight')
        fig2.savefig(save_folder + "Cf " + angle + '.png', bbox_inches='tight')

j=0
prev_path = ""
for path, folders, files in os.walk(filepath):
    # List contain of folder
    for folder_name in folders:
        if folder_name == "results":
            path_folders =path.split("\\")
            if path_folders[-2] != prev_path:
                prev_path = path_folders[-2]
                print(prev_path)

            multiblock = get_multiblock(path + "\\" + folder_name)
            if multiblock is None:
                continue
            #plot each case separately
            plot_cp_cf(multiblock)
            
        

