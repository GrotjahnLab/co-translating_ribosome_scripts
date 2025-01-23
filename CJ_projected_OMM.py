import glob
import os
import click
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from pycurv import TriangleGraph, io
from graph_tool import load_graph

@click.command()
@click.option('--input_omm_graph', required=True, type=str, help='Input OMM graph')
@click.option('--input_cj_imm_csv', required=True, type=str, help='Input CJ_IMM CSV file')
@click.option('--output_omm_surface', required=True, type=str, help='Output OMM surface')
@click.option('--output_omm_graph', required=True, type=str, help='Output OMM graph')

def CJ_projected_OMM(input_omm_graph, input_cj_imm_csv, output_omm_surface, output_omm_graph):
    # OMM graph
    tg1 = TriangleGraph()
    tg1.graph = load_graph(input_omm_graph)

    # Load the CSV file of CJ_IMM
    df = pd.read_csv(input_cj_imm_csv)

    # Use the column named 'OMM_neighbor_index' to get the vertex index from the graph
    OMM_neighbor_index = df['OMM_neighbor_index']
    # Drop the duplicates
    OMM_neighbor_index = OMM_neighbor_index.drop_duplicates()

    # Convert OMM_neighbor_index to a list
    OMM_neighbor_index_list = OMM_neighbor_index.tolist()
    print(f"Number of unique OMM neighbor indices: {len(OMM_neighbor_index_list)}")

    # Create a vertex filter
    vertex_filter = tg1.graph.new_vertex_property("double")
    for vertex in tg1.graph.vertices():
        if vertex in OMM_neighbor_index_list:
            vertex_filter[vertex] = 1
        else:
            vertex_filter[vertex] = 0

    tg1.graph.vertex_properties["CJ_projected_OMM"] = vertex_filter
    surf1 = tg1.graph_to_triangle_poly()
    io.save_vtp(surf1, output_omm_surface)
    tg1.graph.save(output_omm_graph)

if __name__ == '__main__':
    CJ_projected_OMM()


