import glob
import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from pycurv import TriangleGraph, io
from graph_tool import load_graph, GraphView

#cotrans_path = "/data1/atty/cotrans_ribo_bin2/OMM-patches_IMM_dist_csv_1/cotrans_patches_overlap_with_CJ_membrane/*.gt"
#noncotrans_path = "/data1/atty/cotrans_ribo_bin2/noncotrans_OMM-patches_IMM_dist_csv/noncotrans_patches_overlap_with_CJ_membrane/*.gt"
random_path = "/data1/atty/cotrans_ribo_bin2/random_OMM-patches_IMM_dist_csv_1/random_patches_overlap_with_CJ_membrane/*.gt"

#cotrans_data = []
#noncotrans_data = []
random_data = []

#for cotrans_graph in glob.glob(cotrans_path):
#    # Load the cotrans patch graph
#    tg1 = TriangleGraph()
#    tg1.graph = load_graph(cotrans_graph)
#
#    # Extract the matching part of the filename
#    tomogram_name = os.path.basename(cotrans_graph).split('_labels_OMM.AVV_rh100_cotrans_overlap_with_CJ_projected_OMM.gt')[0]
#    print(tomogram_name)
#
#    #fetch the vertex property patch
#    patch = tg1.graph.vp.patch
#    #extract the area of the triangles which have patch = 1
#    area1 = tg1.graph.vp.area.get_array()
#    area1_patch = area1[patch.a == 1]
#    #extract the area of the triangles which have patch = 1 and overlap with CJ_projected_OMM = 1
#    overlap_with_CJ_projected_OMM = tg1.graph.vp.overlap_with_CJ_projected_OMM
#    #if there is 1 in the overlap_with_CJ_projected_OMM, then get the overlap area
#    if 1 in overlap_with_CJ_projected_OMM:
#        area1_overlap = area1[overlap_with_CJ_projected_OMM.a == 1]
#        #calculate the ratio of the area of the triangles which have patch = 1 and overlap with CJ_projected_OMM = 1
#        ratio1 = np.sum(area1_overlap) / np.sum(area1_patch)
#        print(ratio1)
#    else:
#        #if there is no 1 in the overlap_with_CJ_projected_OMM, then the area1_overlap and ratio is 0
#        area1_overlap = 0
#        ratio1 = 0
#        print(ratio1)
#    #append the area1, area1_overlap, ratio, and tomogram_name to the list
#    cotrans_data.append({
#        "tomogram": tomogram_name,
#        "patches_total_area": np.sum(area1_patch),
#        "patches_overlap_area": np.sum(area1_overlap),
#        "area_ratio": ratio1
#    })
#
#for noncotrans_graph in glob.glob(noncotrans_path):
#    # Load the noncotrans patch graph
#    tg2 = TriangleGraph()
#    tg2.graph = load_graph(noncotrans_graph)
#
#    # Extract the matching part of the filename
#    tomogram_name = os.path.basename(noncotrans_graph).split('_labels_OMM.AVV_rh100_noncotrans_overlap_with_CJ_projected_OMM.gt')[0]
#    print(tomogram_name)
#
#    #fetch the vertex property patch
#    patch = tg2.graph.vp.patch
#    #extract the area of the triangles which have patch = 1
#    area2 = tg2.graph.vp.area.get_array()
#    area2_patch = area2[patch.a == 1]
#    #extract the area of the triangles which have patch = 1 and overlap with CJ_projected_OMM = 1
#    overlap_with_CJ_projected_OMM = tg2.graph.vp.overlap_with_CJ_projected_OMM
#    #if there is 1 in the overlap_with_CJ_projected_OMM, then get the overlap area
#    if 1 in overlap_with_CJ_projected_OMM:
#        area2_overlap = area2[overlap_with_CJ_projected_OMM.a == 1]
#        #calculate the ratio of the area of the triangles which have patch = 1 and overlap with CJ_projected_OMM = 1
#        ratio2 = np.sum(area2_overlap) / np.sum(area2_patch)
#        print(ratio2)
#    else:
#        #if there is no 1 in the overlap_with_CJ_projected_OMM, then the area2_overlap and ratio is 0
#        area2_overlap = 0
#        ratio2 = 0
#        print(ratio2)
#    #append the area2, area2_overlap, ratio, and tomogram_name to the list
#    noncotrans_data.append({
#        "tomogram": tomogram_name,
#        "patches_total_area": np.sum(area2_patch),
#        "patches_overlap_area": np.sum(area2_overlap),
#        "area_ratio": ratio2
#    })
#
#

for random_graph in glob.glob(random_path):
    # Load the random patch graph
    tg3 = TriangleGraph()
    tg3.graph = load_graph(random_graph)

    # Extract the matching part of the filename
    tomogram_name = os.path.basename(random_graph).split('_labels_OMM.AVV_rh100_random_overlap_with_CJ_projected_OMM.gt')[0]
    print(tomogram_name)

    #fetch the vertex property patch
    patch = tg3.graph.vp.patch
    #extract the area of the triangles which have patch = 1
    area3 = tg3.graph.vp.area.get_array()
    area3_patch = area3[patch.a == 1]
    #extract the area of the triangles which have patch = 1 and overlap with CJ_projected_OMM = 1
    overlap_with_CJ_projected_OMM = tg3.graph.vp.overlap_with_CJ_projected_OMM
    #if there is 1 in the overlap_with_CJ_projected_OMM, then get the overlap area
    if 1 in overlap_with_CJ_projected_OMM:
        area3_overlap = area3[overlap_with_CJ_projected_OMM.a == 1]
        #calculate the ratio of the area of the triangles which have patch = 1 and overlap with CJ_projected_OMM = 1
        ratio3 = np.sum(area3_overlap) / np.sum(area3_patch)
        print(ratio3)
    else:
        #if there is no 1 in the overlap_with_CJ_projected_OMM, then the area3_overlap and ratio is 0
        area3_overlap = 0
        ratio3 = 0
        print(ratio3)
    #append the area3, area3_overlap, ratio, and tomogram_name to the list
    random_data.append({
        "tomogram": tomogram_name,
        "patches_total_area": np.sum(area3_patch),
        "patches_overlap_area": np.sum(area3_overlap),
        "area_ratio": ratio3
    })

# Create a dataframe from the collected data
#df_cotrans = pd.DataFrame(cotrans_data)
#df_noncotrans = pd.DataFrame(noncotrans_data)
df_random = pd.DataFrame(random_data)

# Save the dataframe to a CSV file
#cotrans_output_path = "/data1/atty/cotrans_ribo_bin2/OMM-patches_IMM_dist_csv_1/cotrans_patches_overlap_with_CJ_membrane/cotrans_patches_overlap_ratio.csv"
#df_cotrans.to_csv(cotrans_output_path, index=False)
#print(f"Data saved to {cotrans_output_path}")
#noncotrans_output_path = "/data1/atty/cotrans_ribo_bin2/noncotrans_OMM-patches_IMM_dist_csv/noncotrans_patches_overlap_with_CJ_membrane/noncotrans_patches_overlap_ratio.csv"
#df_noncotrans.to_csv(noncotrans_output_path, index=False)
#print(f"Data saved to {noncotrans_output_path}")
random_output_path = "/data1/atty/cotrans_ribo_bin2/random_OMM-patches_IMM_dist_csv_1/random_patches_overlap_with_CJ_membrane/random_patches_overlap_ratio.csv"
df_random.to_csv(random_output_path, index=False)
print(f"Data saved to {random_output_path}")

