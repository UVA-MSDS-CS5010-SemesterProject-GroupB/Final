#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#read in most current data set
df= pd.read_csv("RB3_RoiInfo.csv")


#find how many unique neurons (bodyID's) in the right brian
unique_body_ID = df.groupby('bodyId')
print("")
print("There are " + str(len(unique_body_ID)) + " completed and unique neurons on the right side of the hemibrain")


#find the average & the most connected neuron (bodyID) and the number of connections it has
most_connected_body_ID = '' 
num_connect_body_ID = 0
sum_connect_body_ID = 0
for name, group in unique_body_ID:
    sum_connect_body_ID += len(group)
    if len(group)>num_connect_body_ID:
        num_connect_body_ID = len(group)
        most_connected_body_ID = str(name)

ave_connect_body_ID = sum_connect_body_ID/len(unique_body_ID)    

print("The average number of connections per neurons is: " + str(round(ave_connect_body_ID)))
print("The most connected neuron (body ID) is: " + most_connected_body_ID + " with " + str(num_connect_body_ID) + " connections")
print("")




#find how many unique roi regions on the rightside of the Hemibrain
unique_rois = df.groupby('roiRegion')
print("There are " + str(len(unique_rois)) + " unique ROI regions in the right side of the hemibrain")


#find the average connections & the roi region that has the most neurons connected
most_connected_neuron = '' 
num_connect_roi = 0
sum_connect_roi = 0
roi_list=[]
for name, group in unique_rois:
    res=group.count()
    roi_list.append([name, res[1]])
    sum_connect_roi += len(group)
    if len(group)>num_connect_roi:
        num_connect_roi = len(group)
        most_connected_neuron = str(name)
        
ave_connect_roi = sum_connect_roi/len(unique_rois)   



print("The average number of connections for the ROIs is: " + str(round(ave_connect_roi)))
print("The most connected ROI is " + most_connected_neuron + " with " + str(num_connect_roi) + " connections")
print("")



# find MBON-type unique neurons (BodyIDs), how may are they?
df['type_trunc']=df.type.str[:4]
df_MBON = df[df['type_trunc'] == "MBON"] 
unique_body_ID_MBON = df_MBON.groupby('bodyId')
print("There are " + str(len(unique_body_ID_MBON)) + " unique MBONs on the right side of the hemibrain")



#find the average connection of them and compare with the entire group
sum_connect_body_ID_MBON = 0
for name, group in unique_body_ID_MBON:
    sum_connect_body_ID_MBON += len(group)
   
ave_connect_body_ID_MBON = sum_connect_body_ID_MBON/len(unique_body_ID_MBON)    

print("The average number of connections for MBONs is " + str(round(ave_connect_body_ID_MBON))+ ",\n"
      + "while the average number of connections for all right side neurons is " + str(round(ave_connect_body_ID)))
print("")




# average, min, and max number of presynapses and postsynapses for all neurons (bodyIDs) on the right side of the hemibrain
#presynaptic a
average_pre = np.mean(unique_body_ID['pre'].mean())
min_pre = np.min(unique_body_ID['pre'].mean())
max_pre = np.max(unique_body_ID['pre'].mean())
average_post = np.mean(unique_body_ID['post'].mean())
min_post = np.min(unique_body_ID['post'].mean())
max_post = np.max(unique_body_ID['post'].mean())
print("The average number of presynapses/T-Bars is " + str(round(average_pre)) + "\nmin: " + str(min_pre) + "\nmax: " + str(max_pre))
print("The average number of postsynapse/PSD is " + str(round(average_post))+ "\nmin: " + str(min_post) + "\nmax: " + str(max_post))
print("")




#scatter plot for pre vs post
plt.figure()
sns.regplot(unique_body_ID['pre'].mean(), unique_body_ID['post'].mean())  
plt.title('Presynapses / T-Bars versus Postsynapses / PSD')



