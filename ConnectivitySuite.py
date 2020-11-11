# -*- coding: utf-8 -*-
"""
CS 5010 Project:
Connectivity suite
aaw3ff
"""
import numpy as np
import pandas as pd
import neuprint as neuprint
from neuprint import Client
import networkx as nx


authtoken = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImFhdzNmZkB2aXJnaW5pYS5lZHUiLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGg2Lmdvb2dsZXVzZXJjb250ZW50LmNvbS8tb08yX1c1SU9TTW8vQUFBQUFBQUFBQUkvQUFBQUFBQUFBQUEvQU1adXVjbmhIejdMajRXV0FDSkFUbVB4ZWdIMXc3eTFlUS9zOTYtYy9waG90by5qcGc_c3o9NTA_c3o9NTAiLCJleHAiOjE3ODMxNDc2NTZ9.8sqoj4DI3c1dxS0E3xTLeXycWEm5fWeVNxDJy32_iv8'

c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=authtoken)

#Create a connector class to hold connectivity matrices and graph representations
#for roi-roi connections
class Connector:
    #This class has the following attributes: client, the neuprint client used
    #to query the DB, connections and weight matrices stored in pd dataframes, 
    #and directed aand undirected graphs to 
    def __init__(self,client, connections_matrix = None, weight_matrix = None):
        self.imports = 'client from neuprint, networkkx as nx, pandas as pd'
        #Not super resilient to input error as this is not user-facing, but 
        #politely remind us to have the appropriate modules imported
        try:
            self.client = client
            
            #Set placeholders for graph objects until they are computed
            self.undirected_graph = 'No current graphs initialized for this type'
            self.directed_graph = 'No Current graphs initialized for this type'
            
            #If using pre-computed Janelia Data, define connectivity matrices
            if connections_matrix == None:
                #Fetch pre computed values for connectivity from Janelia DB (PD DF)
                self.roi_connectivity = self.client.fetch_roi_connectivity()
                #Compute matrix of number of connections between each ROI in Janelia set
                self.connectivity_matrix = pd.pivot_table(self.roi_connectivity, \
                values = 'count', index='from_roi', columns = 'to_roi',\
                fill_value=0).sort_index(axis=1).sort_index()
                #Compute matrix of weights of connections between each ROI
                self.weight_matrix = pd.pivot_table(self.roi_connectivity, \
                values = 'weight', index='from_roi', columns = 'to_roi',\
                fill_value=0).sort_index(axis=1).sort_index()
            
            #If not using Janelia pre-computed data, accept user inputs
            else:
                self.weight_matrix = weight_matrix
                self. connectivity_matrix = connections_matrix
        except  ModuleNotFoundError:
            print('This object requires the following imports' + self.imports)
    
    #Str function just reports whether the class contains any graphs
    def __str__(self):
        intro = 'This connector object contains: \n'
        if type(self.undirected_graph) == str:
            undirected = 'No undirected graph object \n'
        else:
            undirected_nodes_count = str(len(list(self.undirected_graph)))
            undirected = 'An undirected graph with ' + undirected_nodes_count\
            + ' nodes. \n'
        if type(self.directed_graph) == str:
            directed = 'No directed graph object \n'
        else:
            directed_nodes_count = str(len(list(self.directed_graph)))
            directed = 'A directed graph with ' + directed_nodes_count\
            + ' nodes. \n'    
        return intro + undirected + directed
    
    #Generate an undirected graph, weight = false by default to generate graph based on
    #connections, weight != false for graph based on weights of connections
    def generate_undirected_graph(self, roi_list, weight = False):
        #Generate a non-directional graph object
        self.undirected_graph = nx.Graph()
        self.undirected_graph.add_nodes_from(roi_list)
        
        #If weight is False (default) generate graph based on connectivity
        if weight == False:
            for roi1 in roi_list:
                for roi2 in roi_list:
                    connectivity = self.connectivity_matrix.loc[roi1,roi2]
                    if connectivity > 0:
                        self.undirected_graph.add_edge(roi1, roi2, weight = 1/connectivity)
        #If weight is not False, generate graph based on weight    
        else:
            for roi1 in roi_list:
                for roi2 in roi_list:
                    weight = self.weight_matrix.loc[roi1,roi2]
                    if weight > 0:
                        self.undirected_graph.add_edge(roi1, roi2, weight = 1/weight)
    
    #Generate a directed graph, weight = false by default to generate graph based on
    #connections, weight != false for graph based on weights of connections
    def generate_directed_graph(self, roi_list, weight = False):
        #Generate a directional graph object
        self.directed_graph = nx.DiGraph()
        #Add nodes of interest
        self.directed_graph.add_nodes_from(roi_list)
        
        #If weight is False (default) generate graph based on connectivity
        if weight == False:
            for roi1 in roi_list:
                for roi2 in roi_list:
                    connectivity = self.connectivity_matrix.loc[roi1,roi2]
                    if connectivity > 0:
                        self.directed_graph.add_edge(roi1,roi2, weight = 1/connectivity)
        #If weight is not False, generate graph based on weight    
        else:
            for roi1 in roi_list:
                for roi2 in roi_list:
                    weight = self.weight_matrix.loc[roi1,roi2]
                    if weight > 0:
                        self.directed_graph.add_edge(roi1,roi2, weight = 1/weight)

#%%

#A simple helper function to calculate eigen centralities and return them in a usable form
def eigen_centralities(graph):
    #Get nodes in graph
    nodes = list(graph)
    #Compute centrality for each node, store in a dict
    eigen_dict = nx.eigenvector_centrality(graph, weight = 'weight')
    #Convert dict to data frame and sort it
    eigen_df = pd.DataFrame.from_dict(eigen_dict, orient='index')
    eigen_df = eigen_df.sort_values(by=0, axis = 0, ascending = False)
    
    #Create ordered numpy array from dict for use in color mapping
    eigen_list = []
    for node in nodes:
        eigen_list.append(eigen_dict[node])
    eigen_array = np.array(eigen_list)
    return eigen_array, eigen_df

#A simple helper function to calculate betweenness centralities and return them in a usable form
def between_centralities(graph):
    #Get nodes in graph
    nodes = list(graph)
    #Compute centrality for each node, store in a dict
    between_dict = nx.betweenness_centrality(graph, weight = 'weight')
    #Convert dict to data frame and sort it
    between_df = pd.DataFrame.from_dict(between_dict, orient='index')
    between_df = between_df.sort_values(by=0, axis = 0, ascending = False)
    between_list = []
    
    #Create ordered numpy array from dict for use in color mapping
    for node in nodes:
        between_list.append(between_dict[node])
    between_array = np.array(between_list)
    return between_array, between_df

#Simple helper function to calculate load centralities in a graph
def reach_centralities(graph):
    #Get nodes in graph
    nodes = list(graph)
    #Compute centrality for each node
    reach_list = []
    for node in nodes:
        reach_list.append(nx.local_reaching_centrality(graph, node, weight = 'weight'))
    reach_array = np.array(reach_list)
    return reach_array   
#%%

