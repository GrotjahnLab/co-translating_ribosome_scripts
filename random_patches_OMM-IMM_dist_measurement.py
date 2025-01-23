import click
import starfile
import random
import numpy as np
from pycurv import TriangleGraph, io
from graph_tool import load_graph
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from scipy.spatial import cKDTree
import pandas as pd

@click.command()
@click.option('--input_starfile', type=click.Path(exists=True), help='Path to the input star file')
@click.option('--omm_graph', type=click.Path(exists=True), help='Path to the OMM graph file')
@click.option('--imm_graph', type=click.Path(exists=True), help='Path to the IMM graph file')
@click.option('--tomogram_name', help='Name of the tomogram')
@click.option('--output_csv', default='output.csv', help='Name of the output CSV file')
def main(input_starfile, omm_graph, imm_graph, tomogram_name, output_csv):
    # Load the OMM graph
    tg1 = TriangleGraph()
    tg1.graph = load_graph(omm_graph)
    xyz1 = tg1.graph.vp.xyz.get_2d_array([0, 1, 2]).transpose()

    surface1 = omm_graph[:-3] + "_random_patches.vtp"
    new_omm_graph = omm_graph[:-3] + "_random_patches.gt"
    
    # Load the IMM graph
    tg2 = TriangleGraph()
    tg2.graph = load_graph(imm_graph)
    xyz2 = tg2.graph.vp.xyz.get_2d_array([0, 1, 2]).transpose()
    
    # Load the starfile
    star = starfile.read(input_starfile)
    num_coordinates = len(star)
    print(f"Number of co-translating ribosomes: {num_coordinates}")

    # Convert the xyz1 to a list of tuples
    set_of_coordinates = [tuple(coord) for coord in xyz1]

    random_coordinates = generate_random_coordinates(set_of_coordinates, num_coordinates)

    # Perform the distance measurement and save results
    combined_OMM_coordinates_within_range, combined_min_d_1, output_csv = random_patches_OMM_IMM_dist(set_of_coordinates, random_coordinates, xyz1, xyz2, tomogram_name, output_csv)

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

def calculate_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def generate_random_coordinates(set_of_coordinates, num_coordinates, min_distance=150):
    if len(set_of_coordinates) < num_coordinates:
        raise ValueError("The provided set of coordinates is too small to select the required number of random coordinates.")

    random_coordinates = []
    remaining_coordinates = set_of_coordinates.copy()

    while len(random_coordinates) < num_coordinates and remaining_coordinates:
        candidate = random.choice(remaining_coordinates)
        if all(calculate_distance(candidate, existing_coord) >= min_distance for existing_coord in random_coordinates):
            random_coordinates.append(candidate)
            remaining_coordinates.remove(candidate)

    if len(random_coordinates) < num_coordinates:
        raise ValueError("It is not possible to select the required number of coordinates with the given minimum distance constraint.")

    return random_coordinates

def random_patches_OMM_IMM_dist(set_of_coordinates, random_coordinates, xyz1, xyz2, tomogram_name, output_csv, min_distance=150):
    min_d_1_list = []
    OMM_coordinates_within_range_list = []

    for i in random_coordinates:
        # Calculate the distance between i and each coordinate in xyz1
        distances = np.linalg.norm(xyz1 - np.array(i), axis=1)
        # Find the coordinates within 150 angstroms of i
        indices_within_range = np.where(distances <= min_distance)
        OMM_coordinates_within_range = xyz1[indices_within_range]
        # Measure the distance between the OMM_coordinates_within_range and the IMM
        tree1 = cKDTree(xyz2)
        min_d_1, min_i_1 = tree1.query(OMM_coordinates_within_range, k=1)
        # Append the min_d_1 to a list
        min_d_1_list.append(min_d_1)
        # Append the OMM_coordinates_within_range to a list
        OMM_coordinates_within_range_list.append(OMM_coordinates_within_range)
    
    # Combine the min_d_1 values into a single array
    combined_min_d_1 = np.concatenate(min_d_1_list)
    # Combine the OMM_coordinates_within_range into a single array
    combined_OMM_coordinates_within_range = np.concatenate(OMM_coordinates_within_range_list)

    # Save combined_min_d_1 as a column of df1
    df1 = pd.DataFrame(combined_min_d_1, columns=['min_d_1'])
    # Save combined_OMM_coordinates_within_range as columns of df1
    df1[['OMM_coordinate_x', 'OMM_coordinate_y', 'OMM_coordinate_z']] = pd.DataFrame(combined_OMM_coordinates_within_range)
    # Add tomogram column to df1
    df1['tomogram'] = tomogram_name
    # Save df1 as a csv file
    df1.to_csv(output_csv, index=False)    
    return combined_OMM_coordinates_within_range, combined_min_d_1, output_csv

#function to plot the OMM_coordinates_within_range as a scatter plot
def random_patches_OMM_IMM_dist_plot(combined_OMM_coordinates_within_range, combined_min_d_1, xyz1, xyz2):
    # Define a list of colors for the gradient
    colors = ['#F64F59', '#C471ED', '#12C2E9']
    # Create a custom color map with a gradient
    cmap = LinearSegmentedColormap.from_list('mycmap', colors)
    normalize = Normalize(vmin=min(combined_min_d_1), vmax=max(combined_min_d_1))
    # Create a list of normalized distance values
    normalized_distances = [normalize(d) for d in combined_min_d_1]

    # Plot the OMM_coordinates_within_range as a scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(combined_OMM_coordinates_within_range[:, 0], combined_OMM_coordinates_within_range[:, 1], combined_OMM_coordinates_within_range[:, 2], c=normalized_distances, cmap=cmap)
    # Plot xyz1 as a scatter plot
    ax.scatter(xyz1[:, 0], xyz1[:, 1], xyz1[:, 2], c='grey', marker='o', s=0.01)
    # Plot xyz2 as a scatter plot
    ax.scatter(xyz2[:, 0], xyz2[:, 1], xyz2[:, 2], c='red', marker='o', s=0.01)
    # Add a color bar to the plot
    cbar = fig.colorbar(ScalarMappable(norm=normalize, cmap=cmap), ax=ax, orientation='vertical')
    cbar.set_label('Distance to IMM')
    # Add labels to the plot
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    # Save the plot as a png file
    #plt.savefig('random_patches_OMM_IMM_dist.png')
    plt.show()

if __name__ == '__main__':
    main()
