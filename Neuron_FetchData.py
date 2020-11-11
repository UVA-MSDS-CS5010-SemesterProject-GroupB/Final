#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 21:24:36 2020

@author: jing
"""

import pandas as pd

import neuprint
from neuprint import Client
from neuprint import fetch_primary_rois
from dotenv import load_dotenv
import json
#from neuprint import fetch_roi_completeness
import os
import ast



def connection_setup():
    load_dotenv()
    authtoken = os.environ.get("api-token")
    c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=authtoken)
    
    return c
    

def fetch_PrimaryROI():
    connection_setup()
    primaryRois =fetch_primary_rois()
    return primaryRois
    
def fetch_neurons():
    c = connection_setup()
        
    q = """\
        MATCH (n :Neuron)
        
        RETURN n 
        
    """
    results = c.fetch_custom(q)
    
    
    results.to_csv("first.csv")
    print("Retrieved " +  str(len(results)) + " results from NeuPrint")
    return(len(results))


def Dict2df():
    df=pd.read_csv('first.csv')
    datalist = []
    
    for index, row in df.iterrows():
        
        a = dict(ast.literal_eval(row[1]))
        datalist.append(a)
    
    df = pd.DataFrame(datalist)
    df.to_csv("all_neurons.csv")
    
    successful=True
    
    return successful

# function create a new coloum list by removing unwanted string
def make_unwanted_columns_list(old_column, remove_str):
    new_list=[]
    removed_list=[]
    for name in old_column:
        if not remove_str in name:
            new_list.append(name)
        else:
            removed_list.append(name)
    return new_list, removed_list



def left_or_right_brain(df):
    columns_name = list(df.columns)
    # create a left_brain column list, which contains "(L)"
    new_columns, left_brains = make_unwanted_columns_list(columns_name, "(L)")
    new_columns, right_brain = make_unwanted_columns_list(columns_name, "(R)")
    
    # create on new column "Right_Brain"
    df["Right_Brain"]=""
    
    for index, row in df.iterrows():
        if (True in row[left_brains].values):
            df.loc[index, "Right_Brain"] = False
        if (True in row[right_brain].values):
            df.loc[index, "Right_Brain"] = True
    
    return df
    
    
def Neuron_filter(rb):
    
    # remove Orphan, Leaves, untraced neurons    
    new_rb = rb[rb["statusLabel"].isin(["Traced", "Roughly traced"])] 
#    new_rb.to_csv("rightbrain3.csv")
   
    # the df contains Roi info, which is about 200 columns, so we decided to remove the ROI columns
    c = connection_setup()
    allRoi =neuprint.queries.fetch_all_rois()
    col_list=[]
    for col in new_rb.columns:
        if col in allRoi:
            col_list.append(col)
    
    new_rb_noRoi = new_rb.drop(columns = col_list)
    
#    new_rb.to_csv("rightbrain3_noRoi.csv")
    
    return new_rb, new_rb_noRoi

############################################################################
def checkKey(dict, key): 
    if key in dict: 
        return(dict[key]) 
    # If the key is not in the dictionary, the function returns an empty string 
    else: 
        return ''

def expandSOMALocation(dataFrame, fileName):
    # dataFrame.assign(somaLocation=dataFrame.somaLocation).explode('somaLocation')
    # expand df.tags into its own dataframe
    
    somaLocation_df = pd.DataFrame(data={})
    
    for index, row in dataFrame.iterrows():
        # somaLocation_df = dataFrame['somaLocation'].apply(pd.Series)
        
        if (type(row['somaLocation']) == 'list'):
            # somaLocation_df = row['somaLocation'].apply(pd.Series)
            somaLocation_df = somaLocation_df.rename(columns = lambda x : 'somaLocation_' + str(x))
        
        # print(somaLocation_df)

    pd.concat([dataFrame[:], somaLocation_df[:]], axis=1)
    dataFrame.to_csv(fileName)

def expandROIInfo(dataFrame, fileName):
    rowDFList = []
    
    # Loop over all the rows in customQueryExample
    for index, row in dataFrame.iterrows():
        # print(row['roiInfo'])
        
        row['roiInfo']=json.loads(row['roiInfo'])
        roiInfo_roi_list = list(row['roiInfo'])
        # print(roiInfo_roi_list)
        
        for i in roiInfo_roi_list:
            rowROI = row['roiInfo'][i]
            # print(rowROI)
            # print(i)
            # print(rowROI[index])
            
            somaLocation = row['somaLocation']
            somaX = ''
            somaY = ''
            somaZ = ''
            
            if (type(somaLocation) == 'list'):
                somaX = row['somaLocation'][0]
                somaY = row['somaLocation'][1]
                somaZ = row['somaLocation'][2]
                
    
            newDict = {};
            newDict['key'] = str(row['bodyId']) + str(i)
            newDict['bodyId'] = row['bodyId']
            newDict['instance'] = row['instance']
            newDict['type'] = row['type']
            newDict['pre'] = row['pre']
            newDict['post'] = row['post']
            newDict['status'] = row['status']
            newDict['size'] = row['size']
            newDict['somaLocationX'] = somaX
            newDict['somaLocationY'] = somaY
            newDict['somaLocationZ'] = somaZ
            newDict['somaRadius'] = row['somaRadius']
            newDict['roiRegion'] = i
            newDict['roiRegion - Pre'] = checkKey(rowROI,'pre') 
            newDict['roiRegion - Post'] = checkKey(rowROI,'post')  
            newDict['roiRegion - downstream'] = checkKey(rowROI,'downstream')  
            newDict['roiRegion - upstream'] = checkKey(rowROI,'upstream')  
            
            final = pd.DataFrame(newDict, index=[newDict['key']])
        
            rowDFList.append(final)

    roiExpanded = pd.concat(rowDFList, ignore_index=True)
    roiExpanded.to_csv(fileName)  





if __name__ == "__main__":   
# fetch neuron data from NeuPrint, 

    fetch_neurons()
    
    # convert neuron data to dataframe, dataframe save as a csv, "all_neurons.csv"
    Dict2df()
    
    ##############################################################################
    # clean up the data, only keep right brain neurons 
    
    # read csv file to a dataframe neurons
    neurons = pd.read_csv("all_neurons.csv", index_col=0, low_memory=False)
        
    
    print("============================")    
    print("Please be patient, data processing takes about 15 minutes!")   
    # run function to mark L and R brain
    neurons = left_or_right_brain(neurons)
        
    
    #### filtered right brain, remove left side neurons---------------------------
    print("============================")
    print("number of row before filtering: {}".format(len(neurons)))
    Rneurons = neurons[neurons["Right_Brain"]==True]
    print("number of row after filtering: {}".format(len(Rneurons)))
    
    
    #### futher clean up right brain dataframe to remove unwanted columns
    # remove left brains, labelled with(L)------------------------------------
    print("============================")
    columns_name = list(Rneurons.columns)
    print("orignal dataset has {} columns".format(len(columns_name)))
    new_columns, leftbrain_columns = make_unwanted_columns_list(columns_name, "(L)")
    print("after removing left brain, new dataset has {} columns".format(len(new_columns)))
    
    # remove useless Rois with in Left side of the brain
    Rneurons = Rneurons[new_columns]
    
    # save as a csv file
    Rneurons.to_csv("Right_Brain.csv")
    
    ##### filter out un-traced neurons----------------------------------------
    Rneurons, Rneuron_NoRoi=Neuron_filter(Rneurons)
    Rneurons.to_csv("rightbrain3.csv")
    Rneuron_NoRoi.to_csv("rightbrain3_noRoi.csv")
    ##############################################################################
    # generate RoiInfo file
    print("============================")
    print("generating RoiInfo file, it will take about 15 minutes, almost there... ")
    test =pd.read_csv("rightbrain3.csv",index_col=0, low_memory=False)
    expandROIInfo(test, "RB3_RoiInfo.csv")
    print("Done!")










