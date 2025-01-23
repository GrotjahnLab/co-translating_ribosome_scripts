import glob
import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from pycurv import TriangleGraph, io
from graph_tool import load_graph, GraphView

#patches_graphs_path = "/data1/atty/cotrans_ribo_bin2/OMM-patches_IMM_dist_csv_1/cotrans_patches_membrane/*.gt"
#patches_graphs_path = "/data1/atty/cotrans_ribo_bin2/noncotrans_OMM-patches_IMM_dist_csv/noncotrans_patches_graph/*.gt"
patches_graphs_path = "/data1/atty/cotrans_ribo_bin2/random_OMM-patches_IMM_dist_csv_1/random_patches_membrane/*.gt"
CJ_projected_OMM_graphs_path = "/data1/atty/cotrans_ribo_bin2/CJ_IMM_csv_cleanup_1/expand_CJ_projected_OMM_graph/"

data = []

for patches_graph in glob.glob(patches_graphs_path):
    # Load the cotrans patch graph
    tg1 = TriangleGraph()
    tg1.graph = load_graph(patches_graph)
    
    # Extract the matching part of the filename
    #filename_part = os.path.basename(patches_graph).split('_cotrans_patches')[0]
    #filename_part = os.path.basename(patches_graph).split('_noncotrans_patches')[0]
    filename_part = os.path.basename(patches_graph).split('_random_patches')[0]
    matching_filename = f"{filename_part}_expand_CJ_projected.gt"
    
    # Construct the full path of the matching CJ projected OMM graph
    matching_graph_path = os.path.join(CJ_projected_OMM_graphs_path, matching_filename)
    
    if os.path.exists(matching_graph_path):
        # Load the matching CJ projected OMM graph
        tg2 = TriangleGraph()
        tg2.graph = load_graph(matching_graph_path)
        
        # Fetch the vertex property patch
        patch = tg1.graph.vp.patch
        # Extract the xyz coordinates which have patch = 1
        xyz1 = tg1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
        xyz1_1 = xyz1[patch.a == 1]

        # Fetch the vertex property CJ_projected_OMM
        expanded_CJ_projected_OMM = tg2.graph.vp.expanded_CJ_projected_OMM
        # Extract the xyz coordinates which have CJ_projected_OMM = 1
        xyz2 = tg2.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
        xyz2_1 = xyz2[expanded_CJ_projected_OMM.a == 1]

        # Use KDTree for efficient matching
        tree = cKDTree(xyz2_1)
        distances, i = tree.query(xyz1_1)
        #get the index of the matching triangles
        i = i[distances == 0]
        #create a new vertex property to store the matching triangles
        overlap_with_CJ_projected_OMM = tg1.graph.new_vertex_property("double")
        # if xyz2_1[i] == xyz1, overlap_with_CJ_projected_OMM = 1. Otherwise, overlap_with_CJ_projected_OMM = 0
        for vertex in tg1.graph.vertices():
            vertex_index = int(vertex)
            if np.any(np.all(xyz1[vertex_index] == xyz2_1[i], axis=1)):
                overlap_with_CJ_projected_OMM[vertex] = 1
            else:
                overlap_with_CJ_projected_OMM[vertex] = 0
        # add the new vertex property to the graph
        tg1.graph.vertex_properties["overlap_with_CJ_projected_OMM"] = overlap_with_CJ_projected_OMM
        
        # save the graph
        #directory_path = "/data1/atty/cotrans_ribo_bin2/OMM-patches_IMM_dist_csv_1/cotrans_patches_overlap_with_CJ_membrane/"
        #directory_path = "/data1/atty/cotrans_ribo_bin2/noncotrans_OMM-patches_IMM_dist_csv/noncotrans_patches_overlap_with_CJ_membrane/"
        directory_path = "/data1/atty/cotrans_ribo_bin2/random_OMM-patches_IMM_dist_csv_1/random_patches_overlap_with_CJ_membrane/"
        
        #tg1.graph.save(directory_path + filename_part + "_cotrans_overlap_with_CJ_projected_OMM.gt")
        #tg1.graph.save(directory_path + filename_part + "_noncotrans_overlap_with_CJ_projected_OMM.gt")
        tg1.graph.save(directory_path + filename_part + "_random_overlap_with_CJ_projected_OMM.gt")
        
        #save the graph as a vtp file
        surf = tg1.graph_to_triangle_poly()
        #io.save_vtp(surf, directory_path + filename_part + "_cotrans_overlap_with_CJ_projected_OMM.vtp")
        #io.save_vtp(surf, directory_path + filename_part + "_noncotrans_overlap_with_CJ_projected_OMM.vtp")
        io.save_vtp(surf, directory_path + filename_part + "_random_overlap_with_CJ_projected_OMM.vtp")
        print(f"Overlap with CJ projected OMM for {patches_graph} saved")





##load the cotrans patch graph
#tg1 = TriangleGraph()
#tg1.graph = load_graph(patches_graph)
##Fetch the vertex property patch
#patch = tg1.graph.vp.patch
##extract the xyz coordinates which have patch = 1
#xyz1 = tg1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
#xyz1 = xyz1[patch.a == 1]
#
##load the subgraph
#tg2 = TriangleGraph()
#tg2.graph = load_graph(CJ_projected_OMM_graph)
##Fetch the vertex property CJ_projected_OMM
#CJ_projected_OMM = tg2.graph.vp.CJ_projected_OMM
##extract the xyz coordinates which have CJ_projected_OMM = 1
#xyz2 = tg2.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
#xyz2 = xyz2[CJ_projected_OMM.a == 1]
#
##print the number of coordinates which both in xyz and xyz1
#print(len(np.intersect1d(xyz, xyz1)))
#
