import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy import stats

#load csv files
cotrans_csv = pd.read_csv('/data1/atty/cotrans_ribo_bin2/OMM-patches_IMM_dist_csv_1/cotrans_patches_overlap_with_CJ_membrane/cotrans_patches_overlap_ratio.csv')
#load csv files
non_cotrans_csv = pd.read_csv('/data1/atty/cotrans_ribo_bin2/noncotrans_OMM-patches_IMM_dist_csv/noncotrans_patches_overlap_with_CJ_membrane/noncotrans_patches_overlap_ratio.csv')
#load csv files
random_csv = pd.read_csv('/data1/atty/cotrans_ribo_bin2/random_OMM-patches_IMM_dist_csv_1/random_patches_overlap_with_CJ_membrane/random_patches_overlap_ratio.csv')


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
    print(f"Mann-Whitney U test p-value: {p_u:.8f}, U statistic: {u}, U test stars: {stat_stars}")
    t, p_t = stats.ttest_ind(data1, data2, equal_var=False)
    print(f"T-test p-value: {p_t:.8f}, T statistic: {t:.5f}")
    ks, p_ks = stats.ks_2samp(data1, data2)
    print(f"KS test p-value: {p_ks:.8f}, KS statistic: {ks}")

#perform statistical tests
perform_stat_tests(cotrans_csv['area_ratio'], non_cotrans_csv['area_ratio'], 'cotrans', 'non-cotrans')
perform_stat_tests(cotrans_csv['area_ratio'], random_csv['area_ratio'], 'cotrans', 'random')
perform_stat_tests(non_cotrans_csv['area_ratio'], random_csv['area_ratio'], 'non-cotrans', 'random')

#plot violin plot
fig, ax = plt.subplots()
violin_parts = ax.violinplot([cotrans_csv['area_ratio'], non_cotrans_csv['area_ratio'], random_csv['area_ratio']], showmeans=True, showmedians=True)
violin_parts['bodies'][0].set_facecolor('#75C2F6')
violin_parts['bodies'][1].set_facecolor('#F075AA')
violin_parts['bodies'][2].set_facecolor('#0A6847')
violin_parts['cmeans'].set_color('black')
violin_parts['cmedians'].set_color('black')
violin_parts['cmedians'].set_linestyle('--')
violin_parts['cmaxes'].set_color('black')
violin_parts['cmins'].set_color('black')
violin_parts['cbars'].set_color('black')

#randomly scatter plot data points
cotrans_x = np.random.normal(1, 0.04, size=len(cotrans_csv))
noncotrans_x = np.random.normal(2, 0.04, size=len(non_cotrans_csv))
random_x = np.random.normal(3, 0.04, size=len(random_csv))
ax.scatter(cotrans_x, cotrans_csv['area_ratio'], color='black', alpha=0.5)
ax.scatter(noncotrans_x, non_cotrans_csv['area_ratio'], color='black', alpha=0.5)
ax.scatter(random_x, random_csv['area_ratio'], color='black', alpha=0.5)
plt.xticks([1, 2, 3], ['cotrans', 'non-cotrans', 'random'])
plt.ylabel('The overlap area ratio')
plt.title('The overlap area ratio between ribosome patches and CJ projected OMM')
plt.show()

#print the mean and median of the area ratio
print(f"The mean of the area ratio for cotrans: {np.mean(cotrans_csv['area_ratio'])}")
print(f"The median of the area ratio for cotrans: {np.median(cotrans_csv['area_ratio'])}")
print(f"The mean of the area ratio for non-cotrans: {np.mean(non_cotrans_csv['area_ratio'])}")
print(f"The median of the area ratio for non-cotrans: {np.median(non_cotrans_csv['area_ratio'])}")
print(f"The mean of the area ratio for random: {np.mean(random_csv['area_ratio'])}")
print(f"The median of the area ratio for random: {np.median(random_csv['area_ratio'])}")

