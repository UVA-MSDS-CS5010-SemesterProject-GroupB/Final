#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 21:56:39 2020

@author: jing
"""
import pandas as pd
import unittest
import Neuron_FetchData
import pandas.testing as pd_testing
from neuprint import Client
from dotenv import load_dotenv
import networkx as nx
import numpy as np
from bokeh.plotting import figure, show, output_notebook

#from neuprint import fetch_roi_completeness
import os
authtoken = os.environ.get("api-token")

#Our user-defined code:
from ConnectivitySuite import *
from SkeletonGraph import *

c = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)

class Neuron_Functions_test(unittest.TestCase):

    # testing the function can remove item from a list contraining a mark string
    def test_make_unwanted_columns_list(self):
        old_list = ["(L)", "1", "(R)"]
        remove_str = "(L)"
        
        new_list, removed_list = Neuron_FetchData.make_unwanted_columns_list(old_list, remove_str)
        self.assertEqual(new_list, ["1", "(R)"])
        self.assertEqual(removed_list, ["(L)"])
    
    #testing to see if the matrix is the same as we expected.
    def test_left_or_right_brain(self):
        test=pd.read_csv("test.csv")
        test_result=pd.read_csv("test_result.csv")
        
        # convert object to boolean
        a=Neuron_FetchData.left_or_right_brain(test)
        a=a.convert_dtypes(infer_objects=True)
        test_result=test_result.convert_dtypes(infer_objects=True)
        
        pd_testing.assert_frame_equal(a, test_result)
    
    # testing connection function with NeuPrint
    def test_connection_set(self):
        c=Neuron_FetchData.connection_setup()
        self.assertIsNotNone(c)
    
    
    def test_Neuron_filter(self):
        test=pd.read_csv("test.csv")
        a, b = Neuron_FetchData.Neuron_filter(test)
        len(a.columns)
        
        # because L Roi removed, column number is changed
        self.assertNotEqual(len(b.columns),len(test.columns))
    
    def test_fetch_PrimaryROI(self):
        s = Neuron_FetchData.fetch_PrimaryROI()
        self.assertIsNotNone(s)
        
    ##########################################################################
    # because the dynamic of database, neuron number might vary with time
    # warning fetching may take very long time
    
    
    def test_is_fetch_neuron(self): # testing if neuron fetch function get right number of neuron
        self.assertNotEqual(Neuron_FetchData.fetch_neurons(), 0)
        
    # testing if dictionary can successfully convert to dataframe    
    # warning converting takes 20 ~30 minutes
    
    
    # def test_is_Dict2df_working(self):
    #     print("testing Dict2df function")
    #     self.assertEqual(Neuron_FetchData.Dict2df(), True)
   
class validity_tests(unittest.TestCase):
    #Test that the AL->LH pathway is enriched for MBONS vs either of the two regions
    def test_LH_MBONS(self):
        #Fetch connection between the regions
        q = """\
        MATCH (n :Neuron)
        WHERE  n.roiInfo CONTAINS 'aL(R)' and n.roiInfo CONTAINS 'LH(R)'
        RETURN n.bodyId AS bodyId, n.name AS name, n.pre AS numpre, n.post AS numpost, n.type AS type
        ORDER BY n.pre + n.post DESC
        """
        all_connected = c.fetch_custom(q)
        #Fetch all members of that connection that are mbons
        q = """\
        MATCH (n :Neuron)
        WHERE  n.roiInfo CONTAINS 'aL(R)' and n.roiInfo CONTAINS 'LH(R)' and n.type CONTAINS 'MBON'
        RETURN n.bodyId AS bodyId, n.name AS name, n.pre AS numpre, n.post AS numpost, n.type AS type
        ORDER BY n.type DESC
        """
        connecting_mbon = c.fetch_custom(q)
        #Fetch all mbons in the LH
        q = """\
        MATCH (n :Neuron)
        WHERE  n.roiInfo CONTAINS 'LH(R)' and n.type CONTAINS 'MBON'
        RETURN  n.roiInfo as ROI
        ORDER BY n.type DESC
        """
        results_all_lh_mbon = c.fetch_custom(q)
        #Fetch all neurons in the LH
        q = """\
        MATCH (n :Neuron)
        WHERE  n.roiInfo CONTAINS 'LH(R)'
        RETURN n.bodyId AS bodyId, n.name AS name, n.pre AS numpre, n.post AS numpost, n.type AS type
        ORDER BY n.type DESC
        """
        results_all_lh = c.fetch_custom(q)
        #Fetch all neurons in the AL
        q = """\
        MATCH (n :Neuron)
        WHERE  n.roiInfo CONTAINS 'aL(R)'
        RETURN n.bodyId AS bodyId, n.name AS name, n.pre AS numpre, n.post AS numpost, n.type AS type
        ORDER BY n.type DESC
        """
        results_all_al = c.fetch_custom(q)
        #Fetch all mbons in the Al
        q = """\
        MATCH (n :Neuron)
        WHERE  n.roiInfo CONTAINS 'aL(R)' and n.type CONTAINS 'MBON'
        RETURN n.bodyId AS bodyId, n.name AS name, n.pre AS numpre, n.post AS numpost, n.type AS type, n.roiInfo as ROI
        ORDER BY n.type DESC
        """
        results_all_al_mbon = c.fetch_custom(q)

        #Calculate proportions of mbons in the AL, LH, and connection between the two
        pct_al_mbon = len(results_all_al_mbon)/len(results_all_al)
        pct_lh_mbon = len(results_all_lh_mbon)/len(results_all_lh)
        pct_connected_mbon = len(connecting_mbon)/len(all_connected)
        #The connection should have the highest proportion of MBONs
        self.assertTrue(pct_connected_mbon>pct_lh_mbon and pct_connected_mbon>pct_al_mbon)
        

        
    
    def test_KC_connectivity(self):
        #The APL neuron should connect to every KC, let's test that
        #Fetch everything that the APL connects to
        q = """\
            MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
            WHERE  a.type CONTAINS 'APL'
            RETURN b.bodyId AS bodyId, b.type AS type
            ORDER BY b.bodyId DESC
            """
        APL = c.fetch_custom(q)
        #Fetch all KCs
        q = """\
            MATCH (n :Neuron)
            WHERE  n.type CONTAINS 'KC'
            RETURN n.bodyId AS bodyId, n.type AS type
            ORDER BY n.type DESC
            """
        kenyons = c.fetch_custom(q)
        
        #For each Kenyon cell, assert that it is in the set of neurons the APL
        #connects to.
        for bodyId in list(kenyons['bodyId']):
            self.assertIn(bodyId, list(APL['bodyId']))
        
    def test_ring_connectivity(self):
        #The APL neuron should connect to every KC, let's test that
        #Fetch everything that the APL connects to
        q = """\
            MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
            WHERE  a.roiInfo CONTAINS 'EBr3' and  a.instance CONTAINS 'ring'
            RETURN b.bodyId AS bodyId, b.type AS type
            ORDER BY b.bodyId DESC
            """
        er_three = c.fetch_custom(q)
        #Fetch all KCs
        q = """\
            MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
            WHERE  a.roiIInfo CONTAINS 'EBr1' and a.instance CONTAINS 'ring'
            RETURN b.bodyId AS bodyId, b.type AS type
            ORDER BY b.bodyId DESC
            """
        er_one = c.fetch_custom(q)
        #For each Kenyon cell, assert that it is in the set of neurons the APL
        #connects to.
        for bodyId in list(er_three['bodyId']):
            self.assertNotIn(bodyId, list(er_one['bodyId']))

#Testing our helper functions on simple cannonical graphs with hand- calculated
#values. 
class helper_connector_function_tests(unittest.TestCase):
    def test_eigen_centralities_values(self):
        #Canonical star graph
        test_graph = nx.star_graph(6)
        #Appropriate values
        test_values = np.array([0.7071065,  0.28867525, 0.28867525, 0.28867525,\
                                0.28867525, 0.28867525, 0.28867525])
        #Calculate with function    
        values, dataframe = eigen_centralities(test_graph)
        #Check
        self.assertEqual(test_values.all(), values.all())
    def test_eigen_centralities_df(self):
        #Canonical star graph
        test_graph = nx.star_graph(6)
        #Appropriate values, note series rounds down and array rounds up
        test = {'nodes':[0, 1, 2, 3, 4, 5, 6],'values':[0.707106,  0.28867525, \
                0.28867525, 0.28867525,0.28867525, 0.28867525,\
                0.28867525]}
        test_df = pd.DataFrame(test)
        #Calculate with function    
        values, dataframe = eigen_centralities(test_graph)
        #Check that the dataframe equals what it should
        self.assertEqual(test_df['values'].all(), dataframe[0].all())
    def test_betweenness_centralities_values(self):
        #Canonical star graph
        test_graph = nx.star_graph(6)
        #Appropriate values
        test_values = np.array([1, 0, 0, 0, 0, 0, 0,])
        #Calculate with function    
        values, dataframe =between_centralities(test_graph)
        #Check
        self.assertEqual(test_values.all(), values.all())
    def test_betweenness_centralities_df(self):
        #Canonical star graph
        test_graph = nx.star_graph(6)
        #Appropriate values, note series rounds down and array rounds up
        test = {'nodes':[0, 1, 2, 3, 4, 5, 6],'values':[1, 0, 0, 0, 0, 0, 0,]}
        test_df = pd.DataFrame(test)
        #Calculate with function    
        values, dataframe = between_centralities(test_graph)
        #Check that the dataframe equals what it should
        self.assertEqual(test_df['values'].all(), dataframe[0].all())
    
    def test_reach_centralities(self):
        #test reach centralities on simple triangle, (hand math for this is hard)
        
        #make our triangle
        test_nodes = [1, 2, 3]
        test_graph = nx.Graph()
        test_graph.add_nodes_from(test_nodes)
        test_graph.add_edge(1,2, weight = 1)
        test_graph.add_edge(1,3, weight = 1)
        test_graph.add_edge(2,3, weight = 1)
        
        #Calculate reach centralities for triangle
        reach = reach_centralities(test_graph)
        
        self.assertEqual(reach.all(), np.array([1,1,1]).all())
        
class connector_class_tests(unittest.TestCase):
    
    def test_init_client(self):
        #Test that client is initialized correctly
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        self.assertEqual(connector.client, client)
        
    def test_init_digraph(self):
        #Test that graph is initialized correctly, no input
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        empty_string = 'No Current graphs initialized for this type'
        self.assertEqual(connector.directed_graph, empty_string)
        
    def test_init_graph(self):
        #Test that graph is initialized correctly, no input
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        empty_string = 'No current graphs initialized for this type'
        self.assertEqual(connector.undirected_graph, empty_string)
        
    def test_init_weights(self):
        #Setup the correct dataframe
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        roi_connectivity = client.fetch_roi_connectivity()
        weights = pd.pivot_table(roi_connectivity, \
                values = 'weight', index='from_roi', columns = 'to_roi',\
                fill_value=0).sort_index(axis=1).sort_index()
        #Setup connector, see if it autocomputed
        connector = Connector(client)
        pd.testing.assert_frame_equal(weights, connector.weight_matrix)
        
    def test_init_connections(self):
        #Setup the correct dataframe
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        roi_connectivity = client.fetch_roi_connectivity()
        connections = pd.pivot_table(roi_connectivity, \
                values = 'count', index='from_roi', columns = 'to_roi',\
                fill_value=0).sort_index(axis=1).sort_index()
        #Setup connector, see if it autocomputed
        connector = Connector(client)
        pd.testing.assert_frame_equal(connections, connector.connectivity_matrix)
        
    def test_digraph_nodes(self):
        #test that the nodes on the client directional graph are set up correctly
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        connector.generate_directed_graph(['AL(R)','LH(R)'])
        nodes = list(connector.directed_graph)
        self.assertEqual(nodes, ['AL(R)','LH(R)'])
        
    def test_graph_nodes(self):
        #test that the nodes on the client directional graph are set up correctly
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        connector.generate_undirected_graph(['AL(R)','LH(R)'])
        nodes = list(connector.undirected_graph)
        self.assertEqual(nodes, ['AL(R)','LH(R)'])
        
    def test_digraph_edges(self):
        #test that the edges of the digraph are computed correctly in the mutually
        #regulated example
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        connector.generate_directed_graph(['AL(R)','LH(R)'])
        edges = list(connector.directed_graph.edges)
        self.assertEqual(edges, [('AL(R)', 'AL(R)'), ('AL(R)', 'LH(R)'),('LH(R)', 'AL(R)'), ('LH(R)', 'LH(R)')])
        
    def test_graph_edges(self):
        #test that the edges of the undirectedgraph are computed correctly in the mutually
        #regulated example. Should have one less since direction is not important
        client = Client('neuprint.janelia.org', dataset='hemibrain:v1.0.1', token=authtoken)
        connector = Connector(client)
        connector.generate_undirected_graph(['AL(R)','LH(R)'])
        edges = list(connector.undirected_graph.edges)
        self.assertEqual(edges, [('AL(R)', 'AL(R)'), ('AL(R)', 'LH(R)'), ('LH(R)', 'LH(R)')])    

class skeletonGraph_class_tests(unittest.TestCase):

    def test_skeletonGraph_init_rois(self):
        #test that the rois list parameter for the class is initialized correctly
        test_rois=['ATL(R)','ATL(L)']
        test_cellType='ATL014'   
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        self.assertEqual(test_graph.rois, ['ATL(R)','ATL(L)'])
    
    def test_skeletonGraph_init_cellType(self):
        #test that the cell type parameter for the class is initialized correctly
        test_rois=['ATL(R)','ATL(L)']
        test_cellType='ATL014'   
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        self.assertEqual(test_graph.cellType, 'ATL014')
        
    def test_skeletonGraph_invalid_query_rois(self):
        #test that inputting an invalid ROI will result in validQuery = False
        test_rois = ['InvalidROI','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='SMP573' 
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        self.assertEqual(test_graph.validQuery, False)
        
    def test_skeletonGraph_invalid_query_rois_errorMessafe(self):
        #test that inputting an invalid ROI will generate the correct error message
        test_rois = ['InvalidROI','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='SMP573' 
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        self.assertEqual(test_graph.errorMessage, "Invalid cell type/ROI combination. Make sure your cell types and ROIs are related and rerun the application.") 
        
    def test_skeletonGraph_invalid_query_celltype(self):
        #test that inputting an invalid cell type will result in validQuery = False
        test_rois = ['SMP(R)','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='InvalidType' 
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        self.assertEqual(test_graph.validQuery, False)
        
    def test_skeletonGraph_invalid_query_celltype_errorMessafe(self):
        #test that inputting an invalid cell type will generate the correct error message
        test_rois = ['SMP(R)','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='InvalidType' 
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        self.assertEqual(test_graph.errorMessage, "Invalid cell type/ROI combination. Make sure your cell types and ROIs are related and rerun the application.")   
        
    def test_skeletonGraph_invalid_query_returns_correct_neuron_count(self):
        #test that an invalid query will return a connection count equal to zero
        test_rois=['InvalidROI','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='InvalidType'
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        connection_count = test_graph.neuron_connection_count()
        self.assertEqual(connection_count, 0)
        
    def test_skeletonGraph_valid_query(self):
        #tests that inputting a valid query will result in the validQuery property being set to True
        test_rois=['SMP(R)','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='SMP573' 
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        self.assertEqual(test_graph.validQuery, True)
        
    def test_skeletonGraph_returns_correct_neuron_count(self):
        #test that a valid query for returning neuron connection information returns the right number of records
        test_rois=['SMP(R)','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='SMP573'
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        connection_count = test_graph.neuron_connection_count()
        self.assertEqual(connection_count, 193)
                            
    def test_skeletonGraph_generates_bokeh_plot(self):
        #test that a valid query will set the plot property on our class to be a bokeh figure
        test_rois=['SMP(R)','SCL(R)','SLP(R)','CRE(R)']
        test_cellType='SMP573'
        test_graph = SkeletonGraph(c, test_cellType, test_rois)
        test_graph = test_graph.generateSkeleton()
        self.assertEqual(type(test_graph.plot), type(figure()))    
        
if __name__ == '__main__':
    unittest.main() 
