#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 01:16:36 2020

@author: jing
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

import Neuron_FetchData
import matplotlib as mpl


def ROI_heatmap():           

    roi_list = Neuron_FetchData.fetch_PrimaryROI()
    
    
    new_list, removed_list = Neuron_FetchData.make_unwanted_columns_list(roi_list, "(L)")
    
    
    ####### testing fetch connections################################
    c=Neuron_FetchData.connection_setup()
    roiConnectivity = c.fetch_roi_connectivity()
    roiConnectivity.to_csv("headmap_ini.csv")
    
    
    # remove ROI in the left brain in  "from_roi" column#################
    RightBrain = roiConnectivity[roiConnectivity['from_roi'].isin(new_list)]
    
    
    # remove ROI in the left brain in  "to_roi" column
    RightBrain = RightBrain[RightBrain['to_roi'].isin(new_list)]
    
    RightBrain["weight"] = np.log2(RightBrain["weight"])
    print(RightBrain['weight'].describe()) # get statistical info of 'weight'
    
    RightBrain.to_csv("headmap_df.csv")
    
    
    # create pivot table, row: from_roi, colomn: to_roi
    df = pd.pivot_table(RightBrain, values = 'weight', index='from_roi', columns = 'to_roi', fill_value=0).sort_index(axis=1).sort_index()

            
    ax=plt.axes()        
    sns.set(rc={'figure.figsize':(20,20)}, font_scale=3)
    
    res=sns.heatmap(df, cmap="Blues",square=True, ax=ax)
    
    # res.set_yticklabels(res.get_yticklabels(), fontsize = 18)
    # # res.get_xticklabels
    # # res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 22)
    
    ax.set_title("Right Brain Primary ROIs Connectivity")
    plt.show()
    plt.close()







mpl.rcParams['font.size'] = 9.0

plt.rc('figure', figsize=(20, 15))

# # read right brain data
rb = pd.read_csv("rightbrain3.csv", index_col=0, low_memory=False)



###################################################################
# neuron count in each ROI

# bar plot of neuron number in each primary ROI region
primaryRois =Neuron_FetchData.fetch_PrimaryROI()

df =pd.DataFrame()
temp_list=[]
for col in rb.columns:
    if col in primaryRois:
        res=rb.groupby(col)["bodyId"].count()
        temp_list.append([col, res[1]])

df = pd.DataFrame(temp_list, columns = ["Roi", "count"])
#print(df)

df.plot.bar(x='Roi', y='count', title='Numbers of neurons in ROIs')
plt.show()


##################################################################

# cell types , summary
print("============================")



rb.dropna(subset=['type'], inplace=True)
type_list = ['MBON','KCab', 'KCg','LC','ER', 'LC', 'LA', 'PAM', 'SLP', 'SMP','PVLP'] 

rb['new_type']=""

#for i in rf.iterrows()
for types in type_list:
#    if rb["type"].str.contains(types):
    rb.loc[rb["type"].str.contains(types), 'new_type'] = types

rb.loc[rb["new_type"]=="",'new_type'] = 'others'

cell_type=(rb.groupby(rb["new_type"]))
print(len(cell_type))
print(cell_type["new_type"].describe())

a = cell_type.size()
print(a)  # summary of cell types

# rb['bins'] = pd.cut(rb['size'],bins=[0,17,59,120])
a.plot.pie(title="Cell types of Right Brain" )
print("")
####################################################################
# box plot of neuron size in each primary ROI 
df =pd.DataFrame()
temp_list=[]

# create a new dataframe based on Roi
for col in rb.columns:
    if col in primaryRois:
        rslt_df = rb[rb[col]==True]
        if len(rslt_df) > 0:
            rslt_df['Roi']=col
           
            df = df.append(rslt_df[['Roi', 'size']])
# df = pd.DataFrame(temp_list)
#print(df)

sns.set(rc={'figure.figsize':(20,20)}, font_scale=0.9)
sns.catplot(x='size', y='Roi', kind='box',data=df, showfliers = False)
plt.show()
plt.close()
# ##################################################################

# pre, post synapses , summary
print("============================")



df =pd.DataFrame()
temp_list=[]

# create a new dataframe based on Roi
for col in rb.columns:
    if col in primaryRois:
        rslt_df = rb[rb[col]==True]
        if len(rslt_df) > 0:
            rslt_df['Roi']=col
           
            df = df.append(rslt_df[['Roi', 'pre', 'post']])
# df = pd.DataFrame(temp_list)
#print(df)

# fig = plt.figure()

# ax1 = fig.add_subplot(1, 2, 1)

# ax2 = fig.add_subplot(1, 2, 2)

sns.set(font_scale=1)
sns.catplot(x='pre', y='Roi', kind='box',data=df, showfliers = False)

sns.catplot(x='post', y='Roi', kind='box',data=df, showfliers = False)
plt.show()
plt.close()
###################################################################
# primary ROIs connectivity,
# plot heatmap showing the  connectivity
ROI_heatmap()

###################################################################
# pre synapses vs. post synapses, per neuron

print("============================")
ax1 = df.plot.scatter(x='pre', y='post',  c='DarkBlue')

