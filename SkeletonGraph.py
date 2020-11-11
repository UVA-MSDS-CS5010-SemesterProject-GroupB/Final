#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)

import bokeh
import hvplot.pandas
import holoviews as hv

import bokeh.palettes
from bokeh.plotting import figure, show, output_notebook
output_notebook()

from dotenv import load_dotenv
load_dotenv()

from neuprint import Client
from neuprint import fetch_synapses, NeuronCriteria as NC, SynapseCriteria as SC
from neuprint import fetch_synapse_connections
from neuprint import fetch_neurons
from neuprint import merge_neuron_properties


import os

authtoken = os.environ.get("api-token")

c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token=authtoken)

# Creates a SkeletonGraph class to generate a skeleton graph for a given set of neurons
class SkeletonGraph:
    #This class has the following attributes: 
    # - client: the neuprint client used to query the DB, which 
    #   - This requires an authenticated token to use, which can be obtained by creating an account at:
    #     https://neuprint-test.janelia.org/
    # - neurons: A set of Neuron Criteria (NC), defaulted to an empty
    # - plot: A plot object, initialized to an empty bokeh plot
    def __init__(self, client, cellType, rois, neurons=None, plot=None, validQuery=False, errorMessage='Unexpected Error.'):        
        self.client = client
        self.cellType = cellType
        self.rois = rois
        if (neurons == None):
            self.neurons = []
        if(plot == None):
            self.plot = figure()
        self.validQuery = validQuery
        self.errorMessage = errorMessage
    
    # Str function outputs a message containing the cell type and number of neurons 
    # inputted as parameters during class initialization
    def __str__(self):
        if (self.validQuery):
            if (len(self.neurons.index) > 0):
                joinString = ", "
                roiString =  joinString.join(self.rois)
                info = "Neuron/Synapse Criteria: \n" + "ROIs: " + roiString + "\n" + "Cell Type: " + self.cellType
                return info
            else:
                return "No neurons meet the specified cell type/ROI combination."
        else:
            return self.errorMessage
        
    def neuron_connection_count(self):
        return len(self.neurons)

    # Generates a scatterplot based on associated tbar synapses for the initial set of neurons
    def generateSkeleton(self):
        try:
            neuron_criteria = NC(status='Traced',type=self.cellType, regex=True)
            tbar_criteria = SC(rois=self.rois, type='pre', primary_only=True)
            tbars = fetch_synapses(neuron_criteria, tbar_criteria)
            
            # Plot the synapse positions in a 2D projection
            # self.plot = figure()
            self.plot.scatter(tbars['x'], tbars['z'])
            self.plot.y_range.flipped = True
    
            neuron_criteria = NC(status='Traced',type=self.cellType, regex=True)
            syn_criteria = SC(rois=self.rois, primary_only=True)
            conns = fetch_synapse_connections(neuron_criteria, None, syn_criteria)
            # conns.head()
            
            # Retrieve the types of the post-synaptic neurons
            post_neurons, _ = fetch_neurons(conns['bodyId_post'].unique())
            self.neurons = post_neurons
            conns = merge_neuron_properties(post_neurons, conns, 'type')
            
            top10_counts = conns['type_post'].value_counts().head(10)
            top10_counts
            
            colormap = dict(zip(top10_counts.index, bokeh.palettes.Category10[10]))
            points = conns.query('type_post in @top10_counts.index').copy()
            points['color'] = points['type_post'].map(colormap)
            
            # self.plot = figure()
            self.plot.scatter(points['x_post'], points['z_post'], color=points['color'])
            self.plot.y_range.flipped = True
    
            # Download some skeletons as DataFrames and attach columns for bodyId and color
            skeletons = []
            # print(conns['bodyId_pre'].unique())
            for i, bodyId in enumerate(conns['bodyId_pre'].unique()):
                s = c.fetch_skeleton(bodyId, format='pandas')
                s['bodyId'] = bodyId
                s['color'] = bokeh.palettes.Accent[6][i]
                skeletons.append(s)
            
            # Combine into one big table for convenient processing
            skeletons = pd.concat(skeletons, ignore_index=True)
            skeletons.head()
            
            # Join parent/child nodes for plotting as line segments below.
            segments = skeletons.merge(skeletons, 'inner',
            left_on=['bodyId', 'rowId'],
            right_on=['bodyId', 'link'],
            suffixes=['_child', '_parent'])
            
            # self.plot = figure()
            self.plot.y_range.flipped = True
            
            # Plot skeleton segments (in 2D)
            self.plot.segment(x0='x_child', x1='x_parent',
                      y0='z_child', y1='z_parent',
                      color='color_child',
                      source=segments)
            
            # Also plot the synapses from the above example
            self.plot.scatter(points['x_post'], points['z_post'], color=points['color'])
            self.validQuery = True
            return self
        except:
            self.errorMessage = "Invalid cell type/ROI combination. Make sure your cell types and ROIs are related and rerun the application."
            return self

    