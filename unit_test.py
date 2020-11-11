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
#from neuprint import fetch_roi_completeness
import os
authtoken = os.environ.get("api-token")

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
                
    
if __name__ == '__main__':
    unittest.main() 
