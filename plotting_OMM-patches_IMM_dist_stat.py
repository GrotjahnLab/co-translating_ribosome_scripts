import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import glob

# Paths for CSV files
patches_csvfile_list_nonCHX = glob.glob("/data1/atty/nonCHX/OMM-patches_IMM_dist_csv/*.csv")
patches_csvfile_list_CHX = glob.glob("/data1/atty/cotrans_ribo_bin2/OMM-patches_IMM_dist_csv_1/*.csv")

nonpatches_csvfile_list_nonCHX = glob.glob("/data1/atty/nonCHX/outside_OMM-patches_IMM_dist_csv/*.csv")
nonpatches_csvfile_list_CHX = glob.glob("/data1/atty/cotrans_ribo_bin2/outside_OMM-patches_IMM_dist_csv/*.csv")

all_membrane_csvfile_list_nonCHX_1 = glob.glob("/scratch/atty/YTC009_3_morphometrics_3/processing/*_OMM.AVV_rh100.csv")
all_membrane_csvfile_list_nonCHX_2 = glob.glob("/scratch/atty/YTC025_2_morphometrics_1/processing/*_OMM.AVV_rh100.csv")
all_membrane_csvfile_list_CHX_1 = glob.glob("/scratch/atty/YTC041_1_morphometrics/processing/*_OMM.AVV_rh100.csv")
all_membrane_csvfile_list_CHX_2 = glob.glob("/scratch/atty/YTC041_2_morphometrics/processing/*_OMM.AVV_rh100.csv")
all_membrane_csvfile_list_CHX_3 = glob.glob("/scratch/atty/YTC042_2_morphometrics/processing/*_OMM.AVV_rh100.csv")
all_membrane_csvfile_list_CHX_4 = glob.glob("/scratch/atty/YTC043_1_morphometrics/processing/*_OMM.AVV_rh100.csv")
all_membrane_csvfile_list_CHX_5 = glob.glob("/scratch/atty/YTC043_2_morphometrics/processing/*_OMM.AVV_rh100.csv")

# Combine lists
all_membrane_csvfile_list_nonCHX = all_membrane_csvfile_list_nonCHX_1 + all_membrane_csvfile_list_nonCHX_2
all_membrane_csvfile_list_CHX = all_membrane_csvfile_list_CHX_1 + all_membrane_csvfile_list_CHX_2 + all_membrane_csvfile_list_CHX_3 + all_membrane_csvfile_list_CHX_4 + all_membrane_csvfile_list_CHX_5

random_patches_csvfile_list_nonCHX = glob.glob("/data1/atty/nonCHX/random_OMM-patches_IMM_dist_csv/*.csv")
random_patches_csvfile_list_CHX = glob.glob("/data1/atty/cotrans_ribo_bin2/random_OMM-patches_IMM_dist_csv_1/*.csv")

# Function to calculate peak values from CSV files
def calculate_peak_values(csvfile_list, column_name):
    peak_values = []
    for csvfile in csvfile_list:
        df = pd.read_csv(csvfile)
        hist, bin_edges = np.histogram(df[column_name], bins=100)
        i = np.argmax(hist)
        bin_vals = (bin_edges[i], bin_edges[i+1])
        peak_value = np.mean(bin_vals)
        peak_values.append(peak_value)
    return np.array(peak_values)

# Calculate peak values for each dataset
patches_peak_value_list_nonCHX = calculate_peak_values(patches_csvfile_list_nonCHX, 'min_d_1')
patches_peak_value_list_CHX = calculate_peak_values(patches_csvfile_list_CHX, 'min_d_1')

nonpatches_peak_value_list_nonCHX = calculate_peak_values(nonpatches_csvfile_list_nonCHX, 'min_d_1')
nonpatches_peak_value_list_CHX = calculate_peak_values(nonpatches_csvfile_list_CHX, 'min_d_1')

all_membrane_peak_value_list_nonCHX = calculate_peak_values(all_membrane_csvfile_list_nonCHX, 'IMM_dist')
all_membrane_peak_value_list_CHX = calculate_peak_values(all_membrane_csvfile_list_CHX, 'IMM_dist')

random_patches_peak_value_list_nonCHX = calculate_peak_values(random_patches_csvfile_list_nonCHX, 'min_d_1')
random_patches_peak_value_list_CHX = calculate_peak_values(random_patches_csvfile_list_CHX, 'min_d_1')

# Print mean peak values
print("Mean of the peak value for cotrans (nonCHX): ", np.mean(patches_peak_value_list_nonCHX))
print("Mean of the peak value for cotrans (CHX): ", np.mean(patches_peak_value_list_CHX))
print("Mean of the peak value for noncotrans (nonCHX): ", np.mean(nonpatches_peak_value_list_nonCHX))
print("Mean of the peak value for noncotrans (CHX): ", np.mean(nonpatches_peak_value_list_CHX))
print("Mean of the peak value for all membrane (nonCHX): ", np.mean(all_membrane_peak_value_list_nonCHX))
print("Mean of the peak value for all membrane (CHX): ", np.mean(all_membrane_peak_value_list_CHX))
print("Mean of the peak value for random (nonCHX): ", np.mean(random_patches_peak_value_list_nonCHX))
print("Mean of the peak value for random (CHX): ", np.mean(random_patches_peak_value_list_CHX))

# Perform statistical tests and print results
def perform_stat_tests(data1, data2, label1, label2):
    print(f"\nComparing {label1} and {label2}:")
    stat_stars = " "
    u, p_u = stats.mannwhitneyu(data1, data2)
    if p_u < 0.001:
        stat_stars = "****"
    elif p_u < 0.005:
        stat_stars = "***"
    elif p_u < 0.01:
        stat_stars = "**"
    elif p_u < 0.05:
        stat_stars = "*"
    print(f"Mann-Whitney U test p-value: {p_u:.5f}, U statistic: {u}, U test stars: {stat_stars}")
    t, p_t = stats.ttest_ind(data1, data2, equal_var=False)
    print(f"T-test p-value: {p_t:.5f}, T statistic: {t:.5f}")
    ks, p_ks = stats.ks_2samp(data1, data2)
    print(f"KS test p-value: {p_ks:.5f}, KS statistic: {ks}")

perform_stat_tests(patches_peak_value_list_nonCHX, nonpatches_peak_value_list_nonCHX, "cotrans (nonCHX)", "noncotrans (nonCHX)")
perform_stat_tests(patches_peak_value_list_nonCHX, all_membrane_peak_value_list_nonCHX, "cotrans (nonCHX)", "all membrane (nonCHX)")
perform_stat_tests(nonpatches_peak_value_list_nonCHX, all_membrane_peak_value_list_nonCHX, "noncotrans (nonCHX)", "all membrane (nonCHX)")
perform_stat_tests(patches_peak_value_list_nonCHX, random_patches_peak_value_list_nonCHX, "cotrans (nonCHX)", "random (nonCHX)")
perform_stat_tests(nonpatches_peak_value_list_nonCHX, random_patches_peak_value_list_nonCHX, "noncotrans (nonCHX)", "random (nonCHX)")
perform_stat_tests(random_patches_peak_value_list_nonCHX, all_membrane_peak_value_list_nonCHX, "random (nonCHX)", "all membrane (nonCHX)")

perform_stat_tests(patches_peak_value_list_CHX, nonpatches_peak_value_list_CHX, "cotrans (CHX)", "noncotrans (CHX)")
perform_stat_tests(patches_peak_value_list_CHX, all_membrane_peak_value_list_CHX, "cotrans (CHX)", "all membrane (CHX)")
perform_stat_tests(nonpatches_peak_value_list_CHX, all_membrane_peak_value_list_CHX, "noncotrans (CHX)", "all membrane (CHX)")
perform_stat_tests(patches_peak_value_list_CHX, random_patches_peak_value_list_CHX, "cotrans (CHX)", "random (CHX)")
perform_stat_tests(nonpatches_peak_value_list_CHX, random_patches_peak_value_list_CHX, "noncotrans (CHX)", "random (CHX)")
perform_stat_tests(random_patches_peak_value_list_CHX, all_membrane_peak_value_list_CHX, "random (CHX)", "all membrane (CHX)")

# Combine data for plotting
data = [
    patches_peak_value_list_CHX, nonpatches_peak_value_list_CHX,
    all_membrane_peak_value_list_CHX, random_patches_peak_value_list_CHX,
    patches_peak_value_list_nonCHX, nonpatches_peak_value_list_nonCHX,
    all_membrane_peak_value_list_nonCHX, random_patches_peak_value_list_nonCHX
]

labels = [
    'cotrans (CHX)', 'noncotrans (CHX)', 
    'all membrane (CHX)', 'random (CHX)',
    'cotrans (nonCHX)', 'noncotrans (nonCHX)',
    'all membrane (nonCHX)', 'random (nonCHX)'
]

colors = [
    '#75C2F6', '#B4B4B8',
    '#BA704F', '#0A6847',
    '#75C2F6', '#B4B4B8',
    '#BA704F', '#0A6847'
]

# Create violin plot
fig, ax = plt.subplots(figsize=(12, 8))

violin_parts = ax.violinplot(data, showmeans=True, showmedians=True)

# Color the violin parts
for i, body in enumerate(violin_parts['bodies']):
    body.set_facecolor(colors[i])

violin_parts['cmeans'].set_color('black')
violin_parts['cmedians'].set_color('black')
violin_parts['cmedians'].set_linestyle('--')
violin_parts['cmaxes'].set_color('black')
violin_parts['cmins'].set_color('black')
violin_parts['cbars'].set_color('black')

# Scatter plot the data points
for i, d in enumerate(data):
    y = d
    x = np.random.normal(i + 1, 0.04, len(y))
    ax.scatter(x, y, color='black', s=5, alpha=0.5)

ax.set_xticks(np.arange(1, len(labels) + 1))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Peak Values')
ax.set_title('OMM patches-IMM distance peak values')

plt.tight_layout()
plt.show()
