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

def OMM_patches_to_OMM(input_csv, omm_graph, imm_graph, tomogram_name, output_csv):
    filtered_df = pd.read_csv(input_csv)

    #OMM graph
    tg = TriangleGraph()
    tg.graph = load_graph(omm_graph)
    xyz1 = tg.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    #IMM graph
    tg = TriangleGraph()
    tg.graph = load_graph(imm_graph)
    xyz2 = tg.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    OMM_center_Index = filtered_df['OMM_index']
    center_coordinates = xyz1[OMM_center_Index]

    indicies_within_range_list = []

    for coordinate in center_coordinates:
        # Calculate the distance between each coordinate and the center
        distances = np.linalg.norm(xyz1 - coordinate, axis=1)
        # Find the coordinates within 150 angstroms of the center
        indices_within_range = np.where(distances <= 150)
        # Convert the tuples to a 1xn array
        indices_within_range = np.array(indices_within_range).flatten()
        # Append the array to the list
        indicies_within_range_list.append(indices_within_range)

    # Combine the indices_within_range into a single array
    combined_indices_within_range = np.concatenate(indicies_within_range_list)
    #drop the same values in the combined_indices_within_range
    combined_indices_within_range = np.unique(combined_indices_within_range)
    #drop the combined_indicies_within_range from the xyz1
    OMM_coordinates_outside_range = np.delete(xyz1, combined_indices_within_range, axis=0)
    
    #calculate the distance between the OMM_coordinates_outside_range and the IMM_coordinates
    tree1 = cKDTree(xyz2)
    min_d_1, min_i_1 = tree1.query(OMM_coordinates_outside_range, k = 1)

    #save min_d_1 as a column of df1
    df1 = pd.DataFrame(min_d_1, columns=['min_d_1'])
    #add tomogram column to df1
    df1['tomogram'] = tomogram_name
    #save df1 as a csv file
    df1.to_csv(output_csv, index=False)

if __name__ == '__main__':
    OMM_patches_to_OMM()
