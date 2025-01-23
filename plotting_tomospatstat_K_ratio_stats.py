import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import glob

# paths of the CSV files for cotranslation
cotrans_path_list = glob.glob('/scratch1/users/atty/co-translation/tomospatstat_cotrans/csv_files/*.csv')
# path of the CSV files for non-cotranslation
noncotrans_path_list = glob.glob('/scratch1/users/atty/co-translation/tomospatstat_noncotrans/csv_files/*.csv')

# a function to get the maximun value of K/Kcsv(k) between k=30 and k=40 for each file in the list
def get_max_k(path_list):
    max_k_list = []
    for path in path_list:
        df = pd.read_csv(path)
        df = df[(df['k'] >= 30) & (df['k'] <= 40)]
        max_k = df['K/Kcsr(k)'].max()
        # if the max_k is not NaN, append it to the list
        if not np.isnan(max_k):
            max_k_list.append(max_k)
    return max_k_list


# get the maximum value of K/Kcsv(k) for cotranslation
cotrans_max_k_list = get_max_k(cotrans_path_list)
# get the maximum value of K/Kcsv(k) for non-cotranslation
noncotrans_max_k_list = get_max_k(noncotrans_path_list)

# A function for statistical tests and print results
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
    print(f"T-test p-value: {p_t:.8f}, T statistic: {t:.8f}")
    ks, p_ks = stats.ks_2samp(data1, data2)
    print(f"KS test p-value: {p_ks:.8f}, KS statistic: {ks}")

# perform statistical tests
perform_stat_tests(cotrans_max_k_list, noncotrans_max_k_list, "cotranslation", "non-cotranslation")

violin_data = [cotrans_max_k_list, noncotrans_max_k_list]
violin_labels = ['Co-translation', 'Non-co-translation']
colors = ['#75C2F6', '#F075AA']

# plot violin plot
fig, ax = plt.subplots(figsize=(5,5))
violin_parts = ax.violinplot(violin_data, showmeans=True, showmedians=True)
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
for i, d in enumerate(violin_data):
    y = d
    x = np.random.normal(i + 1, 0.04, len(y))
    ax.scatter(x, y, color='black', alpha=0.5)

ax.set_xticks(np.arange(1, len(violin_labels) + 1))
ax.set_xticklabels(violin_labels)
ax.set_ylim(-1, 45)
ax.set_ylabel('K(r)/Kcsr(r)')
ax.set_title('r = 30-40 nm')

plt.tight_layout()
plt.show()