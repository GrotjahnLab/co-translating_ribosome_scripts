#! /usr/bin/env python
import glob
import os
import click
import starfile
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from pycurv import TriangleGraph, io
from graph_tool import load_graph
from eulerangles import euler2matrix

def process_files(star_file, star_file_1, omm_graph_file, imm_graph_file, tomogram_name, condition):
    # Load star file for the particles centering on ribosome exit tunnel
    df1 = starfile.read(star_file)
    df1['rlnCoordinateX'] = df1['rlnCoordinateX']*df1['rlnPixelSize']
    df1['rlnCoordinateY'] = df1['rlnCoordinateY']*df1['rlnPixelSize']
    df1['rlnCoordinateZ'] = df1['rlnCoordinateZ']*df1['rlnPixelSize']
    #convert the exit tunnel coordinates to numpy array
    a = np.array(df1[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']])
    #read the star file with the original coordinates
    df2 = starfile.read(star_file_1)

    # Load OMM triangle graph
    tg = TriangleGraph()
    tg.graph = load_graph(omm_graph_file)
    xyz1 = tg.graph.vp.xyz.get_2d_array([0,1,2]).transpose()
    # load IMM triangle graph
    tg_1 = TriangleGraph()
    tg_1.graph = load_graph(imm_graph_file)
    xyz2 = tg_1.graph.vp.xyz.get_2d_array([0,1,2]).transpose()

    # Measure the nearest distance between ribosome exit tunnel and OMM
    tree2 = cKDTree(xyz1)
    min_d_2, min_i_2 = tree2.query(a, k = 1)

    # Convert euler angles to rotation matrices
    eulers = -1 * np.array(df2[['rlnAnglePsi', 'rlnAngleTilt', 'rlnAngleRot']])
    rotation_matrices = euler2matrix(eulers,
                                   axes='zyz',
                                   intrinsic=True,
                                   right_handed_rotation=True)
    #convert rotation matrices to normal vectors (the first normal vectors)
    normal_vectors = []
    for rotation_matrix in rotation_matrices:
        normal_vector = rotation_matrix[:, 2]  # Extract third column
        normalized_normal_vector = normal_vector / np.linalg.norm(normal_vector)  # Normalize
        normal_vectors.append(normalized_normal_vector)
    # Convert normal_vectors list to numpy array
    normal_vectors = np.array(normal_vectors)

    #convert rotation matrices to normal vectors (the second normal vectors)
    normal_vectors_1 = []
    for rotation_matrix in rotation_matrices:
        normal_vector_1 = rotation_matrix[:, 1]  # Extract second column
        normalized_normal_vector_1 = normal_vector_1 / np.linalg.norm(normal_vector_1)  # Normalize
        normal_vectors_1.append(normalized_normal_vector_1)
    # Convert normal_vectors list to numpy array
    normal_vectors_1 = np.array(normal_vectors_1)

    #convert rotation matrices to normal vectors (the third normal vectors)
    normal_vectors_2 = []
    for rotation_matrix in rotation_matrices:
        normal_vector_2 = rotation_matrix[:, 0]  # Extract first column
        normalized_normal_vector_2 = normal_vector_2 / np.linalg.norm(normal_vector_2)  # Normalize
        normal_vectors_2.append(normalized_normal_vector_2)
    # Convert normal_vectors list to numpy array
    normal_vectors_2 = np.array(normal_vectors_2)
    
    
    # Get the normal vector of membrane
    angles1 = tg.graph.vp.n_v.get_2d_array([0,1,2]).transpose()

    # Get the 3 relative angles between ribosome and membrane
    neighbor_angles_2to1 = angles1[min_i_2]
    relative_angles_2 = [np.arccos(np.abs(np.dot(normal_vectors[i], neighbor_angles_2to1[i])))*180/np.pi for i in range(len(normal_vectors))]
    relative_angles_2_1 = [np.arccos(np.abs(np.dot(normal_vectors_1[i], neighbor_angles_2to1[i])))*180/np.pi for i in range(len(normal_vectors_1))]
    relative_angles_2_2 = [np.arccos(np.abs(np.dot(normal_vectors_2[i], neighbor_angles_2to1[i])))*180/np.pi for i in range(len(normal_vectors_2))]

    #get the nearest distance between IMM and the OMM which is nearest to the ribosome
    neighbor_OMM_xyz = xyz1[min_i_2]
    # Measure the nearest distance between OMM and IMM
    tree3 = cKDTree(xyz2)
    min_d_3, min_i_3 = tree3.query(neighbor_OMM_xyz, k = 1)

    #get the original ribosome coordinates
    df4 = pd.DataFrame()    
    df4['rlnCoordinateX'] = df2['rlnCoordinateX']*df2['rlnPixelSize']
    df4['rlnCoordinateY'] = df2['rlnCoordinateY']*df2['rlnPixelSize']
    df4['rlnCoordinateZ'] = df2['rlnCoordinateZ']*df2['rlnPixelSize']
    #convert the original ribosome coordinates to numpy array
    b = np.array(df4[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']])
    # Measure the nearest distance between ribosome and OMM
    tree4 = cKDTree(xyz1)
    min_d_4, min_i_4 = tree4.query(b, k = 1)

    #Add the measurements to the dataframe
    df3 = pd.DataFrame()
    df3['rlnMagnification'] = df2['rlnMagnification']
    df3['rlnDetectorPixelSize'] = df2['rlnDetectorPixelSize']
    df3['rlnCoordinateX'] = df2['rlnCoordinateX']
    df3['rlnCoordinateY'] = df2['rlnCoordinateY']
    df3['rlnCoordinateZ'] = df2['rlnCoordinateZ']
    df3['rlnAngleRot'] = df2['rlnAngleRot']
    df3['rlnAngleTilt'] = df2['rlnAngleTilt']
    df3['rlnAnglePsi'] = df2['rlnAnglePsi']
    df3['rlnImageName'] = df2['rlnImageName']
    df3['rlnCtfImage'] = df2['rlnCtfImage']
    df3['rlnRandomSubset'] = df2['rlnRandomSubset']
    df3['rlnPixelSize'] = df2['rlnPixelSize']
    df3['rlnVoltage'] = df2['rlnVoltage']
    df3['rlnSphericalAberration'] = df2['rlnSphericalAberration']
    df3['rlnMicrographName'] = df2['rlnMicrographName']
    relative_angles_2 = pd.Series(relative_angles_2)
    df3['blue'] = relative_angles_2
    relative_angles_2_1 = pd.Series(relative_angles_2_1)
    df3['yellow'] = relative_angles_2_1
    relative_angles_2_2 = pd.Series(relative_angles_2_2)
    df3['red'] = relative_angles_2_2
    min_i_2 = pd.Series(min_i_2)
    df3['OMM_index'] = min_i_2
    min_d_2 = pd.Series(min_d_2)
    df3['ribo_exit_tunnel-OMM_dist'] = min_d_2
    min_i_3 = pd.Series(min_i_3)
    df3['IMM_index'] = min_i_3
    min_d_3 = pd.Series(min_d_3)
    df3['OMM-IMM_dist'] = min_d_3
    df3['tomogram'] = tomogram_name
    df3['condition'] = condition
    df3['ribo_OMM_dist'] = min_d_4

    return df3

@click.command()
@click.option('--star_file', type=click.Path(exists=True), help='Path to the star file with the particles centering on ribosome exit tunnel')
@click.option('--star_file_1', type=click.Path(exists=True), help='Path to the star file with the original coordinates')
@click.option('--omm_graph_file', type=click.Path(exists=True), help='Path to the OMM graph file')
@click.option('--imm_graph_file', type=click.Path(exists=True), help='Path to the IMM graph file')
@click.option('--condition', default='respiratory', help='Condition for the files')
@click.option('--output_file', default='output.csv', help='Name of the output CSV file')
def main(star_file, star_file_1, omm_graph_file, imm_graph_file, condition, output_file):
    star_file = os.path.join(star_file)
    star_file_1 = os.path.join(star_file_1)
    omm_graph_file = os.path.join(omm_graph_file)
    imm_graph_file = os.path.join(imm_graph_file)

    tomogram_name = os.path.splitext(os.path.basename(star_file))[0]
    result = process_files(star_file, star_file_1, omm_graph_file, imm_graph_file, tomogram_name, condition)

    result.index.name = 'ribo_index'
    result.to_csv(output_file, index=True)
    print(f"DataFrame saved to {output_file}")

if __name__ == '__main__':
    main()
