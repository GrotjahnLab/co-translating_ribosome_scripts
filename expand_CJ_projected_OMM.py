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
@click.option('--input_graph', required=True, type=click.Path(exists=True),
              help='Path to the input graph file (CJ_projected_OMM).')
@click.option('--output_directory', required=True, type=click.Path(),
              help='Path to save the modified graph file and surface (expand_CJ_projected_OMM).')

def expand_CJ_projected_OMM(input_graph, output_directory):
    # Load the graph of CJ_projected_OMM
    tg1 = TriangleGraph()
    tg1.graph = load_graph(input_graph)

    # Get the coordinates of the graph
    xyz1 = tg1.graph.vp.xyz.get_2d_array([0, 1, 2]).transpose()

    # Fetch the vertex property CJ_projected_OMM
    CJ_projected_OMM = tg1.graph.vp.CJ_projected_OMM

    # Fetch the coordinates which have CJ_projected_OMM = 1
    xyz2 = xyz1[CJ_projected_OMM.a == 1]

    # Create a new vertex property to store the matching triangles
    expanded_CJ_projected_OMM = tg1.graph.new_vertex_property("double")

    indices_within_range_list = []
    for coordinate in xyz2:
        # Calculate the distance between each coordinate and the CJ_projected_OMM
        distances = np.linalg.norm(xyz1 - coordinate, axis=1)
        # Find the coordinates within the specified distance threshold
        indices_within_range = np.where(distances <= 150)[0]
        indices_within_range_list.append(indices_within_range)

    # Flatten the list of indices arrays
    flattened_indices_within_range_list = [item for sublist in indices_within_range_list for item in sublist]

    # Convert to a numpy array and remove duplicates
    combined_indices_within_range = np.array(flattened_indices_within_range_list)
    unique_combined_indices_within_range = np.unique(combined_indices_within_range)

    # Set the expanded_CJ_projected_OMM property to 1 for the matching triangles
    for vertex in tg1.graph.vertices():
        vertex_index = int(vertex)
        if vertex_index in unique_combined_indices_within_range:
            expanded_CJ_projected_OMM[vertex] = 1
        else:
            expanded_CJ_projected_OMM[vertex] = 0

    # Add the new vertex property to the graph
    tg1.graph.vertex_properties["expanded_CJ_projected_OMM"] = expanded_CJ_projected_OMM

    # get the base of filename of the input graph
    filename_part = os.path.basename(input_graph).split('_CJ_projected_OMM.gt')[0]

    # Save the graph
    tg1.graph.save(output_directory + filename_part + '_expand_CJ_projected.gt')

    # Save the surface
    surf = tg1.graph_to_triangle_poly()
    io.save_vtp(surf, output_directory + filename_part + '_expand_CJ_projected.vtp')
    
    print(f"Saved the expanded CJ_projected_OMM graph and surface to {output_directory}")

if __name__ == '__main__':
    expand_CJ_projected_OMM()
