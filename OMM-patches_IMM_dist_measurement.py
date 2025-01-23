import glob
import os
import click
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from pycurv import TriangleGraph, io
from graph_tool import load_graph
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

@click.command()
@click.option('--input_csv', type=click.Path(exists=True), help='Path to the input CSV file')
@click.option('--omm_graph', type=click.Path(exists=True), help='Path to the OMM graph file')
@click.option('--imm_graph', type=click.Path(exists=True), help='Path to the IMM graph file')
@click.option('--tomogram_name', help='Name of the tomogram')
@click.option('--output_csv', default='output.csv', help='Name of the output CSV file')

def OMM_patches_to_IMM(input_csv, omm_graph, imm_graph, tomogram_name, output_csv):
    df = pd.read_csv(input_csv)
    
    #OMM graph
    tg1 = TriangleGraph()
    tg1.graph = load_graph(omm_graph)
    xyz1 = tg1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    #surface1 = omm_graph[:-3] + "_cotrans_patches.vtp"
    #new_omm_graph = omm_graph[:-3] + "_cotrans_patches.gt"
    surface1 = omm_graph[:-3] + "_noncotrans_patches.vtp"
    new_omm_graph = omm_graph[:-3] + "_noncotrans_patches.gt"

    #IMM graph
    tg2 = TriangleGraph()
    tg2.graph = load_graph(imm_graph)
    xyz2 = tg2.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    OMM_center_Index = df['OMM_index']
    center_coordinates = xyz1[OMM_center_Index]

    min_d_1_list = []
    OMM_coordinates_within_range_list = []

    for coordinate in center_coordinates:
        # Calculate the distance between each coordinate and the center
        distances = np.linalg.norm(xyz1 - coordinate, axis=1)
        # Find the coordinates within 150 angstroms of the center
        indices_within_range = np.where(distances <= 150)
        OMM_coordinates_within_range = xyz1[indices_within_range]
        #Measure the distance between the coodinates within 150 angstroms of the center and the IMM
        tree1 = cKDTree(xyz2)
        min_d_1, min_i_1 = tree1.query(OMM_coordinates_within_range, k = 1)
        min_d_1_list.append(min_d_1)
        OMM_coordinates_within_range_list.append(OMM_coordinates_within_range)

    # Combine the min_d_1 values into a single array
    combined_min_d_1 = np.concatenate(min_d_1_list)
    # Combine the OMM coordinates within range into a single array
    combined_OMM_coordinates_within_range = np.concatenate(OMM_coordinates_within_range_list)

    #save combined_OMM_coordinates_within_range as columns of df1
    df1 = pd.DataFrame(combined_OMM_coordinates_within_range, columns=['OMM_coordinate_x', 'OMM_coordinate_y', 'OMM_coordinate_z'])
    #save combined_min_d_1 as a column of df1
    df1['min_d_1'] = combined_min_d_1
    #add tomogram column to df1
    df1['tomogram'] = tomogram_name
    #save df1 as a csv file
    df1.to_csv(output_csv, index=False)

    # Save patch as a vertex property in the OMM graph
    patch = tg1.graph.new_vertex_property("double")
    #convert combined_OMM_coordinates_within_range to tuple
    combined_OMM_coordinates_within_range = [tuple(coord) for coord in combined_OMM_coordinates_within_range]
    for vertex in tg1.graph.vertices():
        OMM_coordinate = tuple(tg1.graph.vp.xyz[vertex])
        #print(OMM_coordinate)
        patch[vertex] = 1 if OMM_coordinate in combined_OMM_coordinates_within_range else 0
        #print(patch[vertex])

    tg1.graph.vertex_properties["patch"] = patch
    surf1 = tg1.graph_to_triangle_poly()
    io.save_vtp(surf1, surface1)
    tg1.graph.save(new_omm_graph)
    print("Patch saved as a vertex property in the OMM graph")

if __name__ == '__main__':
    OMM_patches_to_IMM()
