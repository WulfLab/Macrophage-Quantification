#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:31:27 2022

@author: stevenvandal
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

directory = 'example/directory'  # set this to the directory containing the files
ctrl_files = ['Lin_Ctrl1_Scan1_CD4CD8only.csv', 'Lin_PosCtrlv2_Scan1_CD4CD8only.csv', 'Lin_Ctrl5_Scan1_CD4CD8only.csv']
ola_files = ['Lin_Ola1_Scan1_CD4CD8only.csv', 'Lin_Ola2_Scan1_CD4CD8only.csv', 'Lin_Ola3_Scan1_CD4CD8only.csv']
ola_cd47_files = ['Lin_OlaCD47_1_Scan1_CD4CD8only.csv', 'Lin_OlaCD47_2_Scan1_CD4CD8only.csv', 'Lin_OlaCD47_3_Scan1_CD4CD8only.csv']
mna_cd47_files = ['Lin_MNACD47_1_Scan1_CD4CD8only.csv', 'Lin_MNACD47_2_Scan1_CD4CD8only.csv', 'Lin_MNACD47_4_Scan1_CD4CD8only.csv', 'Lin_MNACD47_5_Scan1_CD4CD8only.csv']
filelist = [ctrl_files, ola_files, ola_cd47_files, mna_cd47_files]
treatment_names = ['Control', 'Ola', 'Ola+CD47', 'MNA+CD47']
nbins=100

# Load the data from the files and rename the columns to be more user-friendly
data = pd.DataFrame()
for i, files in enumerate(filelist):
    df = pd.concat([pd.read_csv(directory + file) for file in files])
    df.insert(0, 'Treatment', treatment_names[i])
    data = data.append(df)
data = data.drop(columns=['Unnamed: 0', 'Path', 'Sample.Name'])
data = data.rename(columns={'Tissue.Category': 'Tissue', 'Cytoplasm.iFN..Opal.520..Mean..Normalized.Counts..Total.Weighting.': 'IFN', 'Cytoplasm.GzB..Opal.570..Mean..Normalized.Counts..Total.Weighting.': 'GZB', 'Cytoplasm.CD8..Opal.620..Mean..Normalized.Counts..Total.Weighting.': 'CD8', 'Cytoplasm.CD4..Opal.690..Mean..Normalized.Counts..Total.Weighting.': 'CD4', 'section': 'Section'})
data = data.loc[data['Tissue'] == 'Tumor']
data = data.drop(columns=['Tissue'])

# Calculate the percentages of cells that are IFN or GZB positive
count_all_cells = data.groupby(['Treatment', 'Phenotype']).size().rename('Cells')
count_ifn_pos_cells = data.loc[data['IFN'] > 1].groupby(['Treatment', 'Phenotype']).size().rename('IFN+ Cells')
count_gzb_pos_cells = data.loc[data['GZB'] > 1].groupby(['Treatment', 'Phenotype']).size().rename('GZB+ Cells')
count_data = pd.concat([count_all_cells, count_ifn_pos_cells, count_gzb_pos_cells], axis=1).reset_index()
count_data['% IFN+ Cells'] = count_data['IFN+ Cells'] / count_data['Cells']
count_data['% GZB+ Cells'] = count_data['GZB+ Cells'] / count_data['Cells']

# Create a dataframe containing lists of IFN or GZB intensities for only the IFN or GZB positive cells
lists_ifn_pos_cells = data.loc[data['IFN'] > 1].groupby(['Treatment', 'Phenotype']).agg({'IFN': list}).reset_index()
lists_gzb_pos_cells = data.loc[data['GZB'] > 1].groupby(['Treatment', 'Phenotype']).agg({'GZB': list}).reset_index()

# Create a dataframe containing mean intensities of only IFN and GZB positive cells
mean_intensities_ifn_pos_cells = data.loc[data['IFN'] > 1].groupby(['Treatment', 'Phenotype']).agg({'IFN': ['mean', 'std']})
mean_intensities_gzb_pos_cells = data.loc[data['GZB'] > 1].groupby(['Treatment', 'Phenotype']).agg({'GZB': ['mean', 'std']})
mean_intensities_data = pd.concat([mean_intensities_ifn_pos_cells, mean_intensities_gzb_pos_cells], axis=1).reset_index()


def get_histogram(lists_of_cells, treatment, phenotype, signal, nbins):
    list_of_cells = lists_of_cells.loc[(lists_ifn_pos_cells['Treatment'] == treatment) & (lists_ifn_pos_cells['Phenotype'] == phenotype)].reset_index()[signal][0]
    counts, bins = np.histogram(list_of_cells, bins=nbins, range=(min(list_of_cells), max(list_of_cells)))
    bins = [(a + b) / 2 for a, b in zip(bins[0:nbins], bins[1:nbins+1])]
    return counts, bins


# CD4 IFN-positive histograms
cd4_ifn_ctrl_counts, cd4_ifn_ctrl_bins = get_histogram(lists_ifn_pos_cells, 'Control', 'CD4+', 'IFN', nbins)
cd4_ifn_ola_counts, cd4_ifn_ola_bins = get_histogram(lists_ifn_pos_cells, 'Ola', 'CD4+', 'IFN', nbins)
cd4_ifn_ola_cd47_counts, cd4_ifn_ola_cd47_bins = get_histogram(lists_ifn_pos_cells, 'Ola+CD47', 'CD4+', 'IFN', nbins)
cd4_ifn_mna_cd47_counts, cd4_ifn_mna_cd47_bins = get_histogram(lists_ifn_pos_cells, 'MNA+CD47', 'CD4+', 'IFN', nbins)

fig = plt.figure()
plt.plot(cd4_ifn_ctrl_bins, cd4_ifn_ctrl_counts, label='Control', color='blue')
plt.plot(cd4_ifn_ola_bins, cd4_ifn_ola_counts, label='Ola', color='red')
plt.plot(cd4_ifn_ola_cd47_bins, cd4_ifn_ola_cd47_counts, label='Ola + CD47', color='orange')
plt.plot(cd4_ifn_mna_cd47_bins, cd4_ifn_mna_cd47_counts, label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Counts')
plt.title('IFN in CD4 Cells')
plt.xlim([1, 35])
plt.ylim([0, 1200])
plt.legend()
# fig.savefig(directory + '/IFN in CD4 Cells (IFN > 1).png', dpi=300)

# CD8 IFN-positive histograms
cd8_ifn_ctrl_counts, cd8_ifn_ctrl_bins = get_histogram(lists_ifn_pos_cells, 'Control', 'CD8+', 'IFN', nbins)
cd8_ifn_ola_counts, cd8_ifn_ola_bins = get_histogram(lists_ifn_pos_cells, 'Ola', 'CD8+', 'IFN', nbins)
cd8_ifn_ola_cd47_counts, cd8_ifn_ola_cd47_bins = get_histogram(lists_ifn_pos_cells, 'Ola+CD47', 'CD8+', 'IFN', nbins)
cd8_ifn_mna_cd47_counts, cd8_ifn_mna_cd47_bins = get_histogram(lists_ifn_pos_cells, 'MNA+CD47', 'CD8+', 'IFN', nbins)

fig = plt.figure()
plt.plot(cd8_ifn_ctrl_bins, cd8_ifn_ctrl_counts, label='Control', color='blue')
plt.plot(cd8_ifn_ola_bins, cd8_ifn_ola_counts, label='Ola', color='red')
plt.plot(cd8_ifn_ola_cd47_bins, cd8_ifn_ola_cd47_counts, label='Ola + CD47', color='orange')
plt.plot(cd8_ifn_mna_cd47_bins, cd8_ifn_mna_cd47_counts, label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Counts')
plt.title('IFN in CD8 Cells')
plt.xlim([1, 30])
plt.ylim([0, 200])
plt.legend()
# fig.savefig(directory + '/IFN in CD8 Cells (IFN > 1).png', dpi=300)

# CD4 GZB-positive histograms
cd4_gzb_ctrl_counts, cd4_gzb_ctrl_bins = get_histogram(lists_gzb_pos_cells, 'Control', 'CD4+', 'GZB', nbins)
cd4_gzb_ola_counts, cd4_gzb_ola_bins = get_histogram(lists_gzb_pos_cells, 'Ola', 'CD4+', 'GZB', nbins)
cd4_gzb_ola_cd47_counts, cd4_gzb_ola_cd47_bins = get_histogram(lists_gzb_pos_cells, 'Ola+CD47', 'CD4+', 'GZB', nbins)
cd4_gzb_mna_cd47_counts, cd4_gzb_mna_cd47_bins = get_histogram(lists_gzb_pos_cells, 'MNA+CD47', 'CD4+', 'GZB', nbins)

fig = plt.figure()
plt.plot(cd4_gzb_ctrl_bins, cd4_gzb_ctrl_counts, label='Control', color='blue')
plt.plot(cd4_gzb_ola_bins, cd4_gzb_ola_counts, label='Ola', color='red')
plt.plot(cd4_gzb_ola_cd47_bins, cd4_gzb_ola_cd47_counts, label='Ola + CD47', color='orange')
plt.plot(cd4_gzb_mna_cd47_bins, cd4_gzb_mna_cd47_counts, label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Counts')
plt.title('GZB in CD4 Cells')
plt.xlim([1, 35])
plt.ylim([0, 1200])
plt.legend()
# fig.savefig(directory + '/GZB in CD4 Cells (GZB > 1).png', dpi=300)

# CD8 GZB-positive histograms
cd8_gzb_ctrl_counts, cd8_gzb_ctrl_bins = get_histogram(lists_gzb_pos_cells, 'Control', 'CD8+', 'GZB', nbins)
cd8_gzb_ola_counts, cd8_gzb_ola_bins = get_histogram(lists_gzb_pos_cells, 'Ola', 'CD8+', 'GZB', nbins)
cd8_gzb_ola_cd47_counts, cd8_gzb_ola_cd47_bins = get_histogram(lists_gzb_pos_cells, 'Ola+CD47', 'CD8+', 'GZB', nbins)
cd8_gzb_mna_cd47_counts, cd8_gzb_mna_cd47_bins = get_histogram(lists_gzb_pos_cells, 'MNA+CD47', 'CD8+', 'GZB', nbins)

fig = plt.figure()
plt.plot(cd8_gzb_ctrl_bins, cd8_gzb_ctrl_counts, label='Control', color='blue')
plt.plot(cd8_gzb_ola_bins, cd8_gzb_ola_counts, label='Ola', color='red')
plt.plot(cd8_gzb_ola_cd47_bins, cd8_gzb_ola_cd47_counts, label='Ola + CD47', color='orange')
plt.plot(cd8_gzb_mna_cd47_bins, cd8_gzb_mna_cd47_counts, label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Counts')
plt.title('GZB in CD8 Cells')
plt.xlim([1, 30])
plt.ylim([0, 200])
plt.legend()
# fig.savefig(directory + '/GZB in CD8 Cells (IFN > 1).png', dpi=300)

# CUMULATIVE DISTRIBUTION PLOTS

# IFN positive CD4 cells
ifn_cd4_ctrl = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'Control') & (lists_ifn_pos_cells['Phenotype'] == 'CD4+')].reset_index()['IFN'][0]
ifn_cd4_ola = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'Ola') & (lists_ifn_pos_cells['Phenotype'] == 'CD4+')].reset_index()['IFN'][0]
ifn_cd4_ola_cd47 = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'Ola+CD47') & (lists_ifn_pos_cells['Phenotype'] == 'CD4+')].reset_index()['IFN'][0]
ifn_cd4_mna_cd47 = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'MNA+CD47') & (lists_ifn_pos_cells['Phenotype'] == 'CD4+')].reset_index()['IFN'][0]

fig = plt.figure()
plt.plot(sorted(ifn_cd4_ctrl), 1.0 * np.arange(len(ifn_cd4_ctrl)) / (len(ifn_cd4_ctrl) - 1), label='Control', color='blue')
plt.plot(sorted(ifn_cd4_ola), 1.0 * np.arange(len(ifn_cd4_ola)) / (len(ifn_cd4_ola) - 1), label='Ola', color='red')
plt.plot(sorted(ifn_cd4_ola_cd47), 1.0 * np.arange(len(ifn_cd4_ola_cd47)) / (len(ifn_cd4_ola_cd47) - 1), label='Ola + CD47', color='orange')
plt.plot(sorted(ifn_cd4_mna_cd47), 1.0 * np.arange(len(ifn_cd4_mna_cd47)) / (len(ifn_cd4_mna_cd47) - 1), label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Cumulative Probability')
plt.title('IFN in CD4 Cells')
plt.xlim([1, 35])
plt.ylim([0, 1])
plt.legend()
fig.savefig(directory + '/IFN in CD4 Cells (IFN > 1) cumulative.png', dpi=300)

# IFN positive CD8 cells
ifn_cd8_ctrl = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'Control') & (lists_ifn_pos_cells['Phenotype'] == 'CD8+')].reset_index()['IFN'][0]
ifn_cd8_ola = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'Ola') & (lists_ifn_pos_cells['Phenotype'] == 'CD8+')].reset_index()['IFN'][0]
ifn_cd8_ola_cd47 = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'Ola+CD47') & (lists_ifn_pos_cells['Phenotype'] == 'CD8+')].reset_index()['IFN'][0]
ifn_cd8_mna_cd47 = lists_ifn_pos_cells.loc[(lists_ifn_pos_cells['Treatment'] == 'MNA+CD47') & (lists_ifn_pos_cells['Phenotype'] == 'CD8+')].reset_index()['IFN'][0]

fig = plt.figure()
plt.plot(sorted(ifn_cd8_ctrl), 1.0 * np.arange(len(ifn_cd8_ctrl)) / (len(ifn_cd8_ctrl) - 1), label='Control', color='blue')
plt.plot(sorted(ifn_cd8_ola), 1.0 * np.arange(len(ifn_cd8_ola)) / (len(ifn_cd8_ola) - 1), label='Ola', color='red')
plt.plot(sorted(ifn_cd8_ola_cd47), 1.0 * np.arange(len(ifn_cd8_ola_cd47)) / (len(ifn_cd8_ola_cd47) - 1), label='Ola + CD47', color='orange')
plt.plot(sorted(ifn_cd8_mna_cd47), 1.0 * np.arange(len(ifn_cd8_mna_cd47)) / (len(ifn_cd8_mna_cd47) - 1), label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Cumulative Probability')
plt.title('IFN in CD8 Cells')
plt.xlim([1, 30])
plt.ylim([0, 1])
plt.legend()
fig.savefig(directory + '/IFN in CD8 Cells (IFN > 1) cumulative.png', dpi=300)

# GZB positive CD4 cells
gzb_cd4_ctrl = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'Control') & (lists_gzb_pos_cells['Phenotype'] == 'CD4+')].reset_index()['GZB'][0]
gzb_cd4_ola = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'Ola') & (lists_gzb_pos_cells['Phenotype'] == 'CD4+')].reset_index()['GZB'][0]
gzb_cd4_ola_cd47 = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'Ola+CD47') & (lists_gzb_pos_cells['Phenotype'] == 'CD4+')].reset_index()['GZB'][0]
gzb_cd4_mna_cd47 = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'MNA+CD47') & (lists_gzb_pos_cells['Phenotype'] == 'CD4+')].reset_index()['GZB'][0]

fig = plt.figure()
plt.plot(sorted(gzb_cd4_ctrl), 1.0 * np.arange(len(gzb_cd4_ctrl)) / (len(gzb_cd4_ctrl) - 1), label='Control', color='blue')
plt.plot(sorted(gzb_cd4_ola), 1.0 * np.arange(len(gzb_cd4_ola)) / (len(gzb_cd4_ola) - 1), label='Ola', color='red')
plt.plot(sorted(gzb_cd4_ola_cd47), 1.0 * np.arange(len(gzb_cd4_ola_cd47)) / (len(gzb_cd4_ola_cd47) - 1), label='Ola + CD47', color='orange')
plt.plot(sorted(gzb_cd4_mna_cd47), 1.0 * np.arange(len(gzb_cd4_mna_cd47)) / (len(gzb_cd4_mna_cd47) - 1), label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Cumulative Probability')
plt.title('GZB in CD4 Cells')
plt.xlim([1, 35])
plt.ylim([0, 1])
plt.legend()
fig.savefig(directory + '/GZB in CD4 Cells (GZB > 1) cumulative.png', dpi=300)

# GZB positive CD8 cells
gzb_cd8_ctrl = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'Control') & (lists_gzb_pos_cells['Phenotype'] == 'CD8+')].reset_index()['GZB'][0]
gzb_cd8_ola = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'Ola') & (lists_gzb_pos_cells['Phenotype'] == 'CD8+')].reset_index()['GZB'][0]
gzb_cd8_ola_cd47 = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'Ola+CD47') & (lists_gzb_pos_cells['Phenotype'] == 'CD8+')].reset_index()['GZB'][0]
gzb_cd8_mna_cd47 = lists_gzb_pos_cells.loc[(lists_gzb_pos_cells['Treatment'] == 'MNA+CD47') & (lists_gzb_pos_cells['Phenotype'] == 'CD8+')].reset_index()['GZB'][0]

fig = plt.figure()
plt.plot(sorted(gzb_cd8_ctrl), 1.0 * np.arange(len(gzb_cd8_ctrl)) / (len(gzb_cd8_ctrl) - 1), label='Control', color='blue')
plt.plot(sorted(gzb_cd8_ola), 1.0 * np.arange(len(gzb_cd8_ola)) / (len(gzb_cd8_ola) - 1), label='Ola', color='red')
plt.plot(sorted(gzb_cd8_ola_cd47), 1.0 * np.arange(len(gzb_cd8_ola_cd47)) / (len(gzb_cd8_ola_cd47) - 1), label='Ola + CD47', color='orange')
plt.plot(sorted(gzb_cd8_mna_cd47), 1.0 * np.arange(len(gzb_cd8_mna_cd47)) / (len(gzb_cd8_mna_cd47) - 1), label='MNA + CD47', color='green')
plt.xlabel('Intensity')
plt.ylabel('Cumulative Probability')
plt.title('GZB in CD8 Cells')
plt.xlim([1, 30])
plt.ylim([0, 1])
plt.legend()
fig.savefig(directory + '/GZB in CD8 Cells (GZB > 1) cumulative.png', dpi=300)

# Perform the Kolmogorov-Smirnov tests
test_results = pd.DataFrame(columns=['Signal', 'Phenotype', 'Treatment', 'KS statistic', 'p-value'])
ifn_cd4_ola_statistic, ifn_cd4_ola_pvalue = stats.ks_2samp(ifn_cd4_ctrl, ifn_cd4_ola)
test_results.loc[len(test_results.index)] = ['IFN', 'CD4', 'Ola vs. Control', ifn_cd4_ola_statistic, ifn_cd4_ola_pvalue]
ifn_cd4_ola_cd47_statistic, ifn_cd4_ola_cd47_pvalue = stats.ks_2samp(ifn_cd4_ctrl, ifn_cd4_ola_cd47)
test_results.loc[len(test_results.index)] = ['IFN', 'CD4', 'Ola + CD47 vs. Control', ifn_cd4_ola_cd47_statistic, ifn_cd4_ola_cd47_pvalue]
ifn_cd4_mna_cd47_statistic, ifn_cd4_mna_cd47_pvalue = stats.ks_2samp(ifn_cd4_ctrl, ifn_cd4_mna_cd47)
test_results.loc[len(test_results.index)] = ['IFN', 'CD4', 'MNA + CD47 vs. Control', ifn_cd4_mna_cd47_statistic, ifn_cd4_mna_cd47_pvalue]

gzb_cd4_ola_statistic, gzb_cd4_ola_pvalue = stats.ks_2samp(gzb_cd4_ctrl, gzb_cd4_ola)
test_results.loc[len(test_results.index)] = ['GZB', 'CD4', 'Ola vs. Control', gzb_cd4_ola_statistic, gzb_cd4_ola_pvalue]
gzb_cd4_ola_cd47_statistic, gzb_cd4_ola_cd47_pvalue = stats.ks_2samp(gzb_cd4_ctrl, gzb_cd4_ola_cd47)
test_results.loc[len(test_results.index)] = ['GZB', 'CD4', 'Ola + CD47 vs. Control', gzb_cd4_ola_cd47_statistic, gzb_cd4_ola_cd47_pvalue]
gzb_cd4_mna_cd47_statistic, gzb_cd4_mna_cd47_pvalue = stats.ks_2samp(gzb_cd4_ctrl, gzb_cd4_mna_cd47)
test_results.loc[len(test_results.index)] = ['GZB', 'CD4', 'MNA + CD47 vs. Control', gzb_cd4_mna_cd47_statistic, gzb_cd4_mna_cd47_pvalue]

ifn_cd8_ola_statistic, ifn_cd8_ola_pvalue = stats.ks_2samp(ifn_cd8_ctrl, ifn_cd8_ola)
test_results.loc[len(test_results.index)] = ['IFN', 'CD8', 'Ola vs. Control', ifn_cd8_ola_statistic, ifn_cd8_ola_pvalue]
ifn_cd8_ola_cd47_statistic, ifn_cd8_ola_cd47_pvalue = stats.ks_2samp(ifn_cd8_ctrl, ifn_cd8_ola_cd47)
test_results.loc[len(test_results.index)] = ['IFN', 'CD8', 'Ola + CD47 vs. Control', ifn_cd8_ola_cd47_statistic, ifn_cd8_ola_cd47_pvalue]
ifn_cd8_mna_cd47_statistic, ifn_cd8_mna_cd47_pvalue = stats.ks_2samp(ifn_cd8_ctrl, ifn_cd8_mna_cd47)
test_results.loc[len(test_results.index)] = ['IFN', 'CD8', 'MNA + CD47 vs. Control', ifn_cd8_mna_cd47_statistic, ifn_cd8_mna_cd47_pvalue]

gzb_cd8_ola_statistic, gzb_cd8_ola_pvalue = stats.ks_2samp(gzb_cd8_ctrl, gzb_cd8_ola)
test_results.loc[len(test_results.index)] = ['GZB', 'CD8', 'Ola vs. Control', gzb_cd8_ola_statistic, gzb_cd8_ola_pvalue]
gzb_cd8_ola_cd47_statistic, gzb_cd8_ola_cd47_pvalue = stats.ks_2samp(gzb_cd8_ctrl, gzb_cd8_ola_cd47)
test_results.loc[len(test_results.index)] = ['GZB', 'CD8', 'Ola + CD87 vs. Control', gzb_cd8_ola_cd47_statistic, gzb_cd8_ola_cd47_pvalue]
gzb_cd8_mna_cd47_statistic, gzb_cd8_mna_cd47_pvalue = stats.ks_2samp(gzb_cd8_ctrl, gzb_cd8_mna_cd47)
test_results.loc[len(test_results.index)] = ['GZB', 'CD8', 'MNA + CD87 vs. Control', gzb_cd8_mna_cd47_statistic, gzb_cd8_mna_cd47_pvalue]

# Save the data files
# mean_intensities_data.to_csv(directory + '/mean intensities (IFN and GZB > 1).csv')
# count_data.to_csv(directory + '/cell counts.csv')
test_results.to_csv(directory + '/test results.csv')
