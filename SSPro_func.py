# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:59:25 2023

@author: leongs
"""
import os
import pandas as pd
pd.options.display.float_format = "{:,.2f}".format
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize, LogNorm
import seaborn as sns
import re
from tabulate import tabulate
from matplotlib.table import Table
import textwrap
import mpld3
import plotly.graph_objects as go

from collections import OrderedDict
import numpy as np
from datetime import datetime, timedelta
import parse_CX_log
# import dataframe_image as dfi
from matplotlib.ticker import MaxNLocator



#%%
#%% read_log_CX
def read_log_CX(log_path, input_path, output_folder, instrument, start_date):
        
    df_dewpointlog=parse_CX_log.parse(input_path, output_folder)
    # read_log_CX(df_log, df_dewpointlog, instrument, start_date)
    df_log = log_parse_CX(log_path,start_date)
    merged_df=[]
    merged_df = pd.merge_asof(df_dewpointlog, df_log, left_on='Timestamp', right_on='Start_Time')
    
    # Filter the merged_df to keep only rows where 'Timestamp' is within 'Start Time' and 'End Time'
    merged_df = merged_df[(merged_df['Timestamp'] >= merged_df['Start_Time']) & (merged_df['Timestamp'] <= merged_df['End_Time'])]
    
    log_stats={}
    sspro_steps = ["1_Cells and Controls", "2_Scan Chip", "3_RT", "4_PCR1", "5_Tagmentation", "6_P5 Index", "7_P7 Index"]
    params = ["ChipTemp", "ChamberTemp","RH"]
    merged_df = merged_df.reset_index()
    for step_name in sspro_steps:
        log_stats[step_name]={}
        for param in params:
            mean = merged_df[merged_df['Steps'] == step_name][param].mean()
            range_max_min = merged_df[merged_df['Steps'] == step_name][param].max() -  merged_df[merged_df['Steps'] == step_name][param].min()
            log_stats[step_name][param] = {'Mean': mean, 'Range': range_max_min}
    
    # Flatten the nested dictionary
    flat_dict = {key: flatten_dict(value) for key, value in log_stats.items()}
    # Create a DataFrame from the flattened dictionary
    df = pd.DataFrame.from_dict(flat_dict, orient='index')
    # func.plot_table_from_df(df, 'CX1018_20231009', base_width=6)
    # dfi.export(df,master_path+"RH log_"+instrument+'.png')
    return df

#%%
def log_parse_CX(log_path,start_date):
    ''' This function parses logs from CX log path, and pulls out the start and end time of each step. Label which timestamp belongs to which step. Output is a dataframe.
    '''
    keywords = ["MasterPlateDispense_FullChip_35",
                "Chip scan complete",
                "MasterPlateDispense_FullChip_35_filtered",
                "MasterPlateDispense_FullChip_35_index_filtered",
                "MasterPlateDispense_FullChip_100_index_filtered"]
    
    sspro_steps = ["1_Cells and Controls", "2_Scan Chip", "3_RT", "4_PCR1", "5_Tagmentation", "6_P5 Index", "7_P7 Index"]
    
    
    date_formats = ('%Y.%m.%d-%H:%M:%S', '%Y.%m.%d.%H:%M:%S', '%Y.%m.%d', '%Y.%m.%d.%H:%M:%S')
    experiment_start_date = datetime.strptime(start_date, date_formats[2])
    experiment_end_date = (experiment_start_date + timedelta(days=2))
    
    # Initialize lists to store the extracted data
    steps = []
    execution_times = []
    end_times =[]
    
    execution_time_pattern = r"Execution time for ([\w_]+): (\d+:\d+:\d+)"
    
    chip_scan_pattern = r"Chip scan complete : (\d+) seconds"
        
    with open(log_path, 'r') as file:
        for line in file:
            if any(keyword in line for keyword in keywords):
                parts = line.split()
                log_datetime = parts[0]
    
                date = None
                for date_format in date_formats:
                    try:
                        date = datetime.strptime(log_datetime, date_format)
                        break
                    except ValueError:
                        continue
    
                if date is None:
                    continue
    
                if experiment_start_date <= date <= experiment_end_date:
                    execution_time_match = re.search(execution_time_pattern, line)
                    if execution_time_match:
                        execution_keyword, execution_time = execution_time_match.groups()
                        execution_time = pd.to_timedelta(execution_time)
                    
                    chip_scan_match = re.search(chip_scan_pattern, line)
                    if chip_scan_match:
                        execution_keyword = "Scan Chip"
                        execution_time_seconds = int(chip_scan_match.group(1))
                        execution_time = timedelta(seconds=execution_time_seconds)
                    steps.append(execution_keyword)         
                    end_time = date
                    end_times.append(end_time)
                    execution_times.append(execution_time)
                        
                        
    df = pd.DataFrame({'End_Time': end_times, 'Execution_Time': execution_times,'Steps':steps})
    # Calculate the "Start_Time" by subtracting "Execution_Time" from "End_Time"
    df['Start_Time'] = df['End_Time'] - df['Execution_Time']
    df['Steps']=sspro_steps
    df.set_index('Start_Time')
    print(df)
    return df



#%%
def plot_table_from_df(df,index,base_width,save_path=None):
    ''' this is used to plot a table from a dataframe 
    '''

    df = df.applymap(lambda x: textwrap.fill(str(x), width=50) if isinstance(x, str) else round(x, 2) if isinstance(x, (int, float)) else x)

    
    if index =='no index':
        df=df.reset_index(drop=True)

    else: df.index.name=index
    
    df.index = [textwrap.fill(str(idx), width=100) for idx in df.index]

    
    fig, ax = plt.subplots(figsize=(base_width,2),dpi=200)

    if index != 'no index':
        table = ax.table(cellText=df.values, colLabels=df.columns, rowLabels=df.index, cellLoc='center', loc='center')
        w, h = table[0, 1].get_width(), table[0, 1].get_height()
        table.add_cell(0, -1, w, h, text=df.index.name)
    else:     
       table = ax.table(cellText=df.values, colLabels=df.columns, cellLoc='center', loc='center')


    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(3,3)
    # table.auto_set_column_width(col=list(range(len(table.columns))))
    ax.axis('off')


    
    if save_path:
        plt.savefig(save_path,bbox_inches='tight',dpi=200)
        
        # convert fig into html using mpld3
        html_fig = mpld3.fig_to_html(fig)
        
        #write html content to a file
        with open(save_path + ".html", "w") as html_file:
            html_file.write(html_fig)
    plt.show()
    plt.close()



def plot_table_from_df_ordered(df, index, base_width, sort_column="no order", top_rows=20):
    ''' This is used to plot a table from a dataframe.
    '''
    

    # Check if a specific column for sorting is specified
    if sort_column != "no order" and sort_column in df.columns:
        df = df.sort_values(by=sort_column, ascending=False)
    
    df = df.applymap(lambda x: round(x, 2) if isinstance(x, (int, float)) else x)
    df.index.name = index
    
    if index == 'no index':
        df = df.reset_index()
    df=df[0:top_rows]
    fig, ax = plt.subplots(figsize=(base_width, 4), dpi=200)
    table = ax.table(cellText=df.values, colLabels=df.columns, rowLabels=df.index, cellLoc='center', loc='center')
    
    # Adds the index header name
    w, h = table[0, 1].get_width(), table[0, 1].get_height()
    table.add_cell(0, -1, w, h, text=df.index.name)
    
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(3, 3)
    ax.axis('off')

    plt.show()
            
        

#%% SSPRO SPEC1 OLD
def sspro_spec1_defunct(df1):
    negCtrl_BarcodedReadsMean=df1[df1['Sample']=='Neg_Ctrl']['Barcoded_Reads'].mean()
    negCtrl_BarcodedReadsSD=df1[df1['Sample']=='Neg_Ctrl']['Barcoded_Reads'].std()
    negCtrl_avg3SD = negCtrl_BarcodedReadsMean + 3*negCtrl_BarcodedReadsSD
    negCtrl_max = df1[df1['Sample']=='Neg_Ctrl']['Barcoded_Reads'].max()
    negCtrl_min = df1[df1['Sample']=='Neg_Ctrl']['Barcoded_Reads'].min()


    posCtrl_BarcodedReadsMean=df1[df1['Sample']=='Pos_Ctrl']['Barcoded_Reads'].mean()
    posCtrl_BarcodedReadsSD=df1[df1['Sample']=='Pos_Ctrl']['Barcoded_Reads'].std()
    posCtrl_avg3SD = posCtrl_BarcodedReadsMean + 3*posCtrl_BarcodedReadsSD
    posCtrl_max = df1[df1['Sample']=='Pos_Ctrl']['Barcoded_Reads'].max()
    posCtrl_min = df1[df1['Sample']=='Pos_Ctrl']['Barcoded_Reads'].min()
    
    

    sample_df = df1[(df1['Sample'].str.lower() == 'sample') | (df1['Sample'] == 'PBMC')] 
    sample_BarcodedReadsMean = sample_df['Barcoded_Reads'].mean()
    sample_BarcodedReadsSD = sample_df['Barcoded_Reads'].std()
    sample_avg3SD = posCtrl_BarcodedReadsMean + 3*posCtrl_BarcodedReadsSD
    sample_max = sample_df['Barcoded_Reads'].max()
    sample_min = sample_df['Barcoded_Reads'].min()

    result_df = pd.DataFrame({
        'Stats':['Mean','StdDev','AVE+3SD','Max','Min'],
        'Neg_Ctrl': [negCtrl_BarcodedReadsMean, negCtrl_BarcodedReadsSD, negCtrl_avg3SD,negCtrl_max,negCtrl_min],
        'Pos_Ctrl': [posCtrl_BarcodedReadsMean,posCtrl_BarcodedReadsSD,posCtrl_avg3SD,posCtrl_max,posCtrl_min],
        'Cells': [sample_BarcodedReadsMean,sample_BarcodedReadsSD,sample_avg3SD,sample_max,sample_min],
        })
     
    result_df=result_df.set_index('Stats').T
    
    crit1 = negCtrl_BarcodedReadsMean < 0.05* posCtrl_BarcodedReadsMean
    if crit1==True: 
        crit1='Pass'
    else: crit1='Fail'
    crit2 = negCtrl_avg3SD < min(posCtrl_BarcodedReadsMean, posCtrl_BarcodedReadsSD, posCtrl_max, posCtrl_min)
    if crit2==True: 
        crit2='Pass'
    else: crit2='Fail'
    print("crit1: %s, crit2: %s" % (crit1, crit2))
    
    
    # plt.figure(figsize=(4, 3))
    # plt.boxplot([df1[df1['Sample']=='Neg_Ctrl']['Barcoded_Reads'], df1[df1['Sample']=='Pos_Ctrl']['Barcoded_Reads']], labels=['Neg_Ctrl', 'Pos_Ctrl'])
    # plt.title('CX1018_Raw')
    # plt.xlabel('Sample')
    # plt.ylabel('Barcoded_Reads')
    # plt.show()

    
    return result_df, crit1, crit2
#%% SSPRO SPEC1
def sspro_spec1(df1):
    negCtrl_BarcodedReadsMean=df1[df1['Well_Type']=='Neg_Ctrl']['Barcoded_Reads'].mean()
    negCtrl_BarcodedReadsSD=df1[df1['Well_Type']=='Neg_Ctrl']['Barcoded_Reads'].std()
    negCtrl_avg3SD = negCtrl_BarcodedReadsMean + 3*negCtrl_BarcodedReadsSD
    negCtrl_max = df1[df1['Well_Type']=='Neg_Ctrl']['Barcoded_Reads'].max()
    negCtrl_min = df1[df1['Well_Type']=='Neg_Ctrl']['Barcoded_Reads'].min()


    posCtrl_BarcodedReadsMean=df1[df1['Well_Type']=='Pos_Ctrl']['Barcoded_Reads'].mean()
    posCtrl_BarcodedReadsSD=df1[df1['Well_Type']=='Pos_Ctrl']['Barcoded_Reads'].std()
    posCtrl_avg3SD = posCtrl_BarcodedReadsMean + 3*posCtrl_BarcodedReadsSD
    posCtrl_max = df1[df1['Well_Type']=='Pos_Ctrl']['Barcoded_Reads'].max()
    posCtrl_min = df1[df1['Well_Type']=='Pos_Ctrl']['Barcoded_Reads'].min()

    result_df = pd.DataFrame({
        'Stats':['Mean','StdDev','AVE+3SD','Max','Min'],
        'Neg_Ctrl': [negCtrl_BarcodedReadsMean, negCtrl_BarcodedReadsSD, negCtrl_avg3SD,negCtrl_max,negCtrl_min],
        'Pos_Ctrl': [posCtrl_BarcodedReadsMean,posCtrl_BarcodedReadsSD,posCtrl_avg3SD,posCtrl_max,posCtrl_min]
        })
     
    result_df=result_df.set_index('Stats').T
    
    crit1 = negCtrl_BarcodedReadsMean < 0.05* posCtrl_BarcodedReadsMean
    crit2 = negCtrl_avg3SD < min(posCtrl_BarcodedReadsMean, posCtrl_BarcodedReadsSD, posCtrl_max, posCtrl_min)
    
    print("crit1: %s, crit2: %s" % (crit1, crit2))
    
    
    plt.figure(figsize=(4, 3))
    plt.boxplot([df1[df1['Well_Type']=='Neg_Ctrl']['Barcoded_Reads'], df1[df1['Well_Type']=='Pos_Ctrl']['Barcoded_Reads']], labels=['Neg_Ctrl', 'Pos_Ctrl'])
    plt.title('CX1018_Raw')
    plt.xlabel('Sample')
    plt.ylabel('Barcoded_Reads')
    plt.show()

    
    return result_df



#%% SSPRO SPEC2

def sspro_spec2(df1, column):
    quantile_list = [1, 0.975, 0.95, 0.925, 0.7, 0.5, 0.15, 0.125, 0.1, 0.075, 0.0]
    quantiles = [(f"{q * 100}%", df1[column].quantile(q)) for q in quantile_list]

    # # just select samples (no control)
    # df1= df1[(df1["Well_Type"] == "PBMC") | (df1["Well_Type"] == "sample") ]
    
    # Create a new DataFrame with quantiles as rows
    quantile_df = pd.DataFrame(quantiles, columns=['Quantile', column])

    range_val_1 = quantile_df.loc[quantile_df['Quantile'] == '92.5%', column].values[0] / quantile_df.loc[quantile_df['Quantile'] == '7.5%', column].values[0]

    range_val_2 = quantile_df.loc[quantile_df['Quantile'] == '95.0%', column].values[0] / quantile_df.loc[quantile_df['Quantile'] == '10.0%', column].values[0]
    
    range_val_3 = quantile_df.loc[quantile_df['Quantile'] == '97.5%', column].values[0] / quantile_df.loc[quantile_df['Quantile'] == '12.5%', column].values[0]
    
    range_val_4 = quantile_df.loc[quantile_df['Quantile'] == '100%', column].values[0] / quantile_df.loc[quantile_df['Quantile'] == '15.0%', column].values[0]


    val_70 = quantile_df.loc[quantile_df['Quantile'] == '70.0%', column].values[0]

    failspec = []
    for i, range_val in enumerate([range_val_4, range_val_3, range_val_1, range_val_2]):
        if range_val > 10:
            failspec.append(f'range_val_{i+1}')
            
        # If any variable is greater than 10, print the list of variables
    if failspec:
        print("> 10-fold difference:", failspec)
        
    print("100/15: %s, 97.5/12.5: %s, 92.5/7.5 : %s, 95.0/10.0: %s, 70.0: %s" % (range_val_4, range_val_3, range_val_1, range_val_2, val_70))
    

    return quantile_df, range_val_1, range_val_2,range_val_3,range_val_4, val_70

#%%


def sspro_spec3(df1_list, column):
    quantile_list = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    
    # Create a list of 'Sample' values from df1_list
    sample_values = [df['Sample'].values[0] for df in df1_list]
    # sample_values = [df['Sample'].unique()[0] for df in df1_list]

    
    # Initialize a dictionary to store quantile values for each 'Sample'
    quantile_data = {'Quantile': []}
    quantile_data['Quantile'].extend([f"{q * 100}%" for q in quantile_list])

    # Populate the dictionary with quantile values for each 'Sample'
    for i, sample in enumerate(sample_values):
        quantile_key = f"%Intergenic_Reads\nSample {sample}"
        quantile_values = [df1_list[i][column].quantile(q) for q in quantile_list]
        quantile_data[quantile_key] = quantile_values

    # Create the final DataFrame from the dictionary
    result_df = pd.DataFrame(quantile_data)

    tenth_percentile_list = []
    for df in df1_list:
        sample = df['Sample'].values[0]
        tenth_percentile = (df[column][df[column] <= 10].count() / len(df[column]))
        print(f"{sample}: {tenth_percentile * 100}%")
        tenth_percentile_list.append({'Sample': sample, 'Tenth_Percentile': tenth_percentile * 100})
    tenth_percentile_df=pd.DataFrame(tenth_percentile_list)
    return result_df,tenth_percentile_df


#%%


def sspro_spec4(df1):
    barcoded_reads = df1['Barcoded_Reads'].sum()
    barcoded_reads = df1['Trimmed_Reads'].sum()
    mapped_reads = df1['Mapped_Reads'].sum()
    unmapped_reads = df1['Unmapped_Reads'].sum()
    percent_unmapped_reads=unmapped_reads*100/barcoded_reads
    
    result_df = pd.DataFrame({
        'Stats':['Sample \nBarcoded Reads','Sample \nTrimmed Reads','Sample \nMapped Reads','Sample \nUnmapped Reads','% Unmapped Reads'],
        'Values':[barcoded_reads,barcoded_reads,mapped_reads,unmapped_reads,percent_unmapped_reads]
        })
    result_df=result_df.set_index('Stats').T
    pd.set_option('display.float_format', '{:.2f}'.format)
    return result_df

def sscomparison_spec1(df,shasta,CX,shasta_100K,CX_100K):
    CX_no_genes_median = df[df['Sample'] == CX]['No_of_Genes'].median()
    shasta_no_genes_median = df[df['Sample'] == shasta]['No_of_Genes'].median()
    CX100K_no_genes_median = df[df['Sample'] == CX_100K]['No_of_Genes'].median()
    shasta100K_no_genes_median = df[df['Sample'] == shasta_100K]['No_of_Genes'].median()
    
    result_df = pd.DataFrame({
        'Stats':[f'No_of_Genes \n {shasta}_median',f'No_of_Genes \n {CX}_median',f'{shasta} to {CX} ratio'],
        'Raw':[shasta_no_genes_median,CX_no_genes_median,shasta_no_genes_median/CX_no_genes_median],
        '100K':[shasta100K_no_genes_median,CX100K_no_genes_median,shasta100K_no_genes_median/CX100K_no_genes_median]
        })
    result_df=result_df.set_index('Stats').T

    return result_df

def sscomparison_spec1_nocx(df,shasta,shasta_100K):
    shasta_no_genes_median = df[df['Sample'] == shasta]['No_of_Genes'].median()
    shasta100K_no_genes_median = df[df['Sample'] == shasta_100K]['No_of_Genes'].median()
    
    result_df = pd.DataFrame({
        'Stats':[f'No_of_Genes \n {shasta}_median'],
        'Raw':[shasta_no_genes_median],
        '100K':[shasta100K_no_genes_median]
        })
    result_df=result_df.set_index('Stats').T
    print(result_df)



    return result_df


def sscomparison_spec1_doublet_old(df,Beta1human,Beta1mouse,CXhuman,CXmouse):
    shasta_human_no_genes_median = df[df['Sample'] == Beta1human]['No_of_Genes'].median()
    shasta_mouse_no_genes_median = df[df['Sample'] == Beta1mouse]['No_of_Genes'].median()
    CX_human_no_genes_median = df[df['Sample'] == CXhuman]['No_of_Genes'].median()
    CX_mouse_no_genes_median = df[df['Sample'] == CXmouse]['No_of_Genes'].median()
    
    result_df = pd.DataFrame({
        'Stats':['No_of_Genes \n Beta2_median','No_of_Genes \n CX1018_median','Beta1 to CX1018 ratio'],
        'K562 at 100K':[shasta_human_no_genes_median,CX_human_no_genes_median,shasta_human_no_genes_median/CX_human_no_genes_median],
        '3T3 at 100K':[shasta_mouse_no_genes_median,CX_mouse_no_genes_median, shasta_mouse_no_genes_median/CX_mouse_no_genes_median]
        })
    result_df=result_df.set_index('Stats').T

    return result_df


def sscomparison_spec1_doublet(df,shasta,cx):
    shasta_no_genes_median = df[df['Instrument'] == shasta]['No_of_Genes'].median()
    CX_no_genes_median = df[df['Instrument'] == cx]['No_of_Genes'].median()
    
    result_df = pd.DataFrame({
        'Stats':['No_of_Genes \n Beta2_median','No_of_Genes \n CX1018_median','Beta2 to CX1018 ratio'],
        'All cells at 100K':[shasta_no_genes_median,CX_no_genes_median,shasta_no_genes_median/CX_no_genes_median]
        })
    result_df=result_df.set_index('Stats').T

    return result_df


#%%
def plot_dataframe(df, column_widths=None, fontsize=10, decimal_places=2):
    fig, ax = plt.subplots(figsize=(10, 6),dpi=300)

    if column_widths:
        ax.axis('tight')
        ax.axis('off')
        table_data = []
        column_headers = [''] + list(df.columns)
        table_data.append(column_headers)  # Add the column headers
        table_data.extend(
            [([index] + [f'{value:.{decimal_places}f}' if isinstance(value, (int, float)) else value for value in row]) for
             index, row in df.iterrows()])

        table = ax.table(cellText=table_data, cellLoc='center', loc='center', colWidths=column_widths)
        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
    else:
        ax.axis('off')
        table_data = [[index] + row for index, row in df.iterrows()]  # Include the index as the first column
        column_headers = [''] + list(df.columns)
        table = ax.table(cellText=table_data, colLabels=column_headers, cellLoc='center', loc='center')

    plt.show()


#%%


def search_files_in_directory(directory_path):
    try:
        files_in_directory = os.listdir(directory_path)
        found_files = [os.path.join(directory_path, file) for file in files_in_directory if "WellList" in file] #well list path as str       
        if not found_files:
            return None  # No files with "WellList" found
        df = pd.read_csv(found_files[0],sep='\t')
        return df
    except FileNotFoundError:
        return None  # Directory not found


def cleanup_barcode(df):
    if 'Barcode' in df.columns:
       df['Barcode'] = df['Barcode'].str.replace('+', '',regex=False)
       
#%%

def gen_heatmap_individual_test(df_list, param, group_column, colour, base_width):
    # Initialize variables to store global minimum and maximum values
    global_min = float('inf')
    global_max = float('-inf')

    # Iterate through all dataframes in df_list to find global min and max
    for df in df_list:
        min_val = df[param].min()
        max_val = df[param].max()
        if min_val < global_min:
            global_min = min_val
        if max_val > global_max:
            global_max = max_val

    # for df in df_list:
    #     df.fillna(0, inplace=True)
    new_df = pd.DataFrame([(i, j) for i in range(72) for j in range(72)], columns=['Row', 'Col'])

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        for group_value in unique_groups:
            fig, ax = plt.subplots()
            merged_df = pd.merge(new_df, df[df[group_column] == group_value], on=['Row', 'Col'], how='left')

            heatmap_data = merged_df.pivot(index='Row', columns='Col', values=param)

            # Take the logarithm of the data for log scaling
            if colour == 'log':  # if log yes, plot log
                # Set the same log scale limits for all plots
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 72, 72, 0], origin='upper', norm=LogNorm(vmin=global_min, vmax=global_max))
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
            else:
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0])
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Value of {param}')

            ax.set_xlabel('Col')
            ax.set_ylabel('Row')
            plt.suptitle(f'{param}')

            # Add y-axis marks
            y_ticks = list(range(0, 75, 10))  # You can customize these values as needed
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_ticks)

            # Add x-axis marks
            x_ticks = list(range(0, 75, 10))
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticks)

            plt.show()
            plt.close()  # Close the figure to release resources

def gen_heatmap_individual(df_list, param, group_column, colour, base_width):
    # for df in df_list:
    #     df.fillna(0, inplace=True)
    new_df = pd.DataFrame([(i, j) for i in range(72) for j in range(72)], columns=['Row', 'Col'])

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        for group_value in unique_groups:
            fig, ax = plt.subplots()
            merged_df = pd.merge(new_df, df[df[group_column] == group_value], on=['Row', 'Col'], how='left')

            heatmap_data = merged_df.pivot(index='Row', columns='Col', values=param)


            # Take the logarithm of the data for log scaling
            if colour == 'log':  #if log yes, plot log
                # heatmap_data = np.log1p(heatmap_data)
                im = ax.imshow(heatmap_data, cmap='bwr',extent=[0, 72, 72, 0], origin='upper', norm=LogNorm())
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
            else: 
                # im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 73, 73, 0], origin='upper')
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0,72,72, 0])
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Value of {param}')
                # ax.legend().set_visible(False)

            
            ax.set_xlabel('Col')
            ax.set_ylabel('Row')
            plt.suptitle(f'{param}')
            
            # Add y-axis marks
            y_ticks = list(range(0,75,5))  # You can customize these values as needed
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_ticks)

            # Add x-axis marks
            x_ticks = list(range(0, 75, 5))
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticks)

            plt.show()
            plt.close()  # Close the figure to release resources


#%%
def gen_heatmap_individual_filtered(df_list, param, group_column, colour, base_width):
    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        for group_value in unique_groups:
            fig, ax = plt.subplots()
            
            #filter for intergenic reads > 0.2
            filtered_df = df[(df[group_column] == group_value) & (df[param] > 20)]
            heatmap_data = filtered_df.pivot(index='Row', columns='Col', values=param)
            
            if colour == 'log':
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 71, 71, 0], origin='upper', norm=LogNorm())
                ax.set_title(f'{group_column}: {group_value}')
                cbar = fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
            else: 
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 71, 71, 0], origin='upper')
                ax.set_title(f'{group_column}: {group_value}')
                cbar = fig.colorbar(im, ax=ax, label=f'Value of {param} > 20%')
                ax.legend().set_visible(False)
            
            cbar.ax.yaxis.labelpad = 10

            ax.set_xlabel('Col')
            ax.set_ylabel('Row')
            plt.suptitle(f'{param}')

            plt.show()
            plt.close()
#%%
def gen_heatmap_filtered(df_list, param, group_column, colour, base_width):
    # Calculate the number of rows and columns for subplots
    num_rows = 2  # You can adjust the number of rows as needed
    num_cols = 2  # You can adjust the number of columns as needed

    # Create a single figure with subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(base_width * num_cols, base_width * num_rows))
    fig.subplots_adjust(hspace=0.1)  # Adjust the vertical spacing between subplots

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        # Determine the subplot position based on the index
        row_position = index // num_cols
        col_position = index % num_cols

        ax = axes[row_position, col_position]

        for group_value in unique_groups:
            # Apply >20% filtering
            filtered_df = df[(df[group_column] == group_value) & (df[param] > 20)]
            heatmap_data = filtered_df.pivot(index='Row', columns='Col', values=param)
            
            if colour == 'log':
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 72, 72, 0], origin='upper', norm=LogNorm())
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
            else: 
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0], origin='upper')
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Value of {param} > 20%')
            
            ax.set_xlabel('Col', fontsize=15)
            ax.set_ylabel('Row', fontsize=15)
            plt.suptitle(f'{param} > 20%', fontsize=20, y=0.9)

    # Close any unused subplots
    for i in range(len(df_list), num_rows * num_cols):
        row_position = i // num_cols
        col_position = i % num_cols
        fig.delaxes(axes[row_position, col_position])

    plt.show()



#%%

def NO_gen_heatmap(df_list, param, group_column, colour, base_width):
    # Step 1: Calculate global min and max
    global_min = float('inf')
    global_max = float('-inf')
    
    for df in df_list:
        unique_groups = df[group_column].unique()
        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)
            data_min = heatmap_data.min().min()
            data_max = heatmap_data.max().max()
            if data_min < global_min:
                global_min = data_min
            if data_max > global_max:
                global_max = data_max
    
    # Step 2: Create figure and subplots
    num_rows = 2  # You can adjust the number of rows as needed
    num_cols = 2  # You can adjust the number of columns as needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(base_width * num_cols, base_width * num_rows))
    fig.subplots_adjust(hspace=0.1)  # Adjust the vertical spacing between subplots

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        # Determine the subplot position based on the index
        row_position = index // num_cols
        col_position = index % num_cols

        ax = axes[row_position, col_position]

        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)
            
            if colour == 'log':
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 72, 72, 0], origin='upper', norm=LogNorm(vmin=global_min, vmax=global_max))
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
            else: 
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0], vmin=global_min, vmax=global_max)
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Value of {param}')
            
            ax.set_xlabel('Col', fontsize=15)
            ax.set_ylabel('Row', fontsize=15)
            plt.suptitle(f'{param}', fontsize=20, y=0.9)

    # Close any unused subplots
    for i in range(len(df_list), num_rows * num_cols):
        row_position = i // num_cols
        col_position = i % num_cols
        fig.delaxes(axes[row_position, col_position])

    plt.show()




#%%
def OLD_gen_heatmap(df_list, param, group_column, colour, base_width):
    # Calculate the number of rows and columns for subplots
    num_rows = 2  # You can adjust the number of rows as needed
    num_cols = 2  # You can adjust the number of columns as needed

    # Create a single figure with subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(base_width * num_cols, base_width * num_rows))
    fig.subplots_adjust(hspace=0.1)  # Adjust the vertical spacing between subplots

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        # Determine the subplot position based on the index
        row_position = index // num_cols
        col_position = index % num_cols

        ax = axes[row_position, col_position]

        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)
            
            if colour == 'log':
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 72, 72, 0], origin='upper', norm=LogNorm())
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
                # cbar = fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
                # cbar.ax.set_aspect('auto')  # Set colorbar aspect ratio to auto
            else: 
                # im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0], origin='upper')
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0])
                # im = ax.imshow(heatmap_data, cmap=colour)
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Value of {param}')
                # cbar = fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
                # cbar.ax.set_aspect('auto')  # Set colorbar aspect ratio to auto
            
            ax.set_xlabel('Col',fontsize=15)
            ax.set_ylabel('Row',fontsize=15)
            # ax.titlepad = 1 #padding between title and plots
            plt.suptitle(f'{param}',fontsize=20,y=0.9)

    # Close any unused subplots
    for i in range(len(df_list), num_rows * num_cols):
        row_position = i // num_cols
        col_position = i % num_cols
        fig.delaxes(axes[row_position, col_position])

    plt.show()
    #%%

def gen_heatmap(df_list, param, group_column, colour, base_width=8):
    # Step 1: Calculate global min and max
    global_min = float('inf')
    global_max = float(0)

    for df in df_list:
        unique_groups = df[group_column].unique()
        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)

            # Ensure the data is numeric and handle non-numeric or NaN values
            if heatmap_data.empty:
                continue
            
            heatmap_data = heatmap_data.to_numpy()  # Convert to numpy array for min/max calculation
            
            data_min = np.nanmin(heatmap_data)
            data_max = np.nanmax(heatmap_data)

            if data_min < global_min:
                global_min = data_min
            if data_max > global_max:
                global_max = data_max

    # Step 2: Calculate the number of rows and columns for subplots
    num_rows = (len(df_list) + 1) // 2  # Use ceiling division for rows
    num_cols = min(len(df_list), 2)  # Use two columns or adjust based on your preference

    # Create a single figure with subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(base_width * num_cols, base_width * num_rows))
    fig.subplots_adjust(hspace=0.3)  # Adjust the vertical spacing between subplots

    # Flatten axes in case of multiple rows/columns
    axes = axes.flatten()

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        ax = axes[index]  # Get the current subplot axis

        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)

            # Ensure the data is numeric and handle non-numeric or NaN values
            if heatmap_data.empty:
                continue
            
            heatmap_data = heatmap_data.to_numpy()  # Convert to numpy array for plotting
            
            # Handling logarithmic scaling
            if colour == 'log':
                norm = LogNorm(vmin=max(global_min, 1e-10), vmax=global_max)
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 72, 72, 0], origin='upper', norm=norm)
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Log Value of {param}')
            else:
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0], vmin=global_min, vmax=global_max)
                ax.set_title(f'{group_column}: {group_value}')
                fig.colorbar(im, ax=ax, label=f'Value of {param}')
            
            ax.set_xlabel('Col', fontsize=15)
            ax.set_ylabel('Row', fontsize=15)

    # Close any unused subplots
    for i in range(len(df_list), num_rows * num_cols):
        fig.delaxes(axes[i])

    # Set the title for the whole figure
    plt.suptitle(f'{param}', fontsize=20, y=0.9)

    # # Save the plot to a buffer and encode it as base64
    # img_buf = BytesIO()
    # fig.savefig(img_buf, format='png', bbox_inches='tight', pad_inches=0)
    # img_buf.seek(0)
    # img_base64 = base64.b64encode(img_buf.getvalue()).decode('utf-8')
    plt.show()

    plt.close(fig)  # Close the figure to free memory




#%%
def plot_barchart(df_list, param, row_or_col, group_column,base_width):
    # Iterate through each DataFrame in the list
    for index, df in enumerate(df_list):
        # Extract unique values of 'row_or_col' for the current DataFrame
        unique_values = df[row_or_col].unique()
        
        # Initialize empty lists to store the x-axis and y-axis data
        x_values = []
        y_values = []

        for value in unique_values:
            # Filter the DataFrame for the current 'row_or_col' value
            filtered_df = df[df[row_or_col] == value]
            
            # Calculate the mean of 'Barcoded_Reads' for the filtered DataFrame
            barcoded_read_mean = filtered_df['Barcoded_Reads'].mean()
            # Append the 'row_or_col' value and the mean of 'Barcoded_Reads' to the lists
            x_values.append(value)
            y_values.append(barcoded_read_mean)

        # Create a bar chart for the current DataFrame
        x_ticks = np.arange(len(x_values))  
        # plt.figure(figsize=(8, 6))  # Adjust the figure size if needed
        # plt.figure(figsize=(base_width*0.75, base_width)) 

        plt.bar(x_ticks, y_values, width=0.5, align='center')
        
        # Set the x-axis tick labels to match the 'row_or_col' values
# Adjust the x-axis tick spacing for legibility
        plt.gca().xaxis.set_major_locator(plt.MultipleLocator(base=10))
        plt.gca().yaxis.set_major_locator(plt.MultipleLocator(base=50000))

        # Add labels and title
        plt.xlabel(row_or_col)
        plt.tick_params(axis='x', direction='out', length=7, width=2)
        plt.ylabel(f'Mean {param}')
        plt.title(f'{df[group_column][0]}: Mean {param} vs. {row_or_col}')
        plt.show()


#%%

def plot(df, target_cols, group_column, specific_order, p_type='violin', base_width=10, base_height=6, rotate=0, hue='Sample'):

    ncol = len(target_cols)
    if ncol == 1:
        fig, axs = plt.subplots()
    else:
        fig, axs = plt.subplots(ncols=ncol)    

    for i, column in enumerate(target_cols):
        ax = None
        if p_type == 'violin':
            if ncol == 1:
                ax = sns.violinplot(data=df, x=group_column, y=column, width=0.5, hue=hue,inner='box')
            else:
                ax = sns.violinplot(data=df, x=group_column, y=column, ax=axs[i],inner='box')
        elif p_type == 'box':
            if ncol == 1:
                ax = sns.boxplot(data=df, x=group_column, y=column)
            else:
                ax = sns.boxplot(data=df, x=group_column, y=column, ax=axs[i])

        if ax is None:
            continue
        
        # Remove the legend
        # ax.legend([], [], frameon=False)

        # medians = df.groupby("Sample")[column].median()
        medians = df.groupby(group_column)[column].median()
        medians = medians.reindex(specific_order)

        # Add median values
        ax.set_title(f'{column}')
        
        for xtick, median in enumerate(medians):
            ax.text(xtick + 0.45, median + 0.01, f"{median:.2f}", horizontalalignment='right', color='black', fontsize=12)


        if rotate == 0:
            ax.set_xticklabels(ax.get_xticklabels())
        else:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotate, ha='right')

    fig.set_figwidth(ncol * base_width)
    fig.set_figheight(base_height)

#%%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def failplot(df, target_cols, group_column, specific_order=None,
         p_type='violin', base_width=6, base_height=6,
         rotate=0, hue='Sample'):

    ncol = len(target_cols)
    fig, axs = plt.subplots(ncols=ncol if ncol > 1 else 1, figsize=(ncol * base_width, base_height))
    
    if ncol == 1:
        axs = [axs]

    for i, column in enumerate(target_cols):
        ax = axs[i]

        # ----- Plotting -----
        if p_type == 'violin':
            sns.violinplot(
                data=df,
                x=group_column,
                y=column,
                hue=hue if hue in df.columns else None,
                order=specific_order,
                width=0.9,
                dodge=True,
                inner=None,  # remove internal median/box
                ax=ax
            )
        elif p_type == 'box':
            sns.boxplot(
                data=df,
                x=group_column,
                y=column,
                hue=hue if hue in df.columns else None,
                order=specific_order,
                ax=ax
            )
        else:
            raise ValueError("p_type must be 'violin' or 'box'")

        # ----- Rotate x labels -----
        if rotate != 0:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotate, ha='right')

        # ----- Overlay median with pointplot (always correct) -----
        if hue in df.columns:
            sns.pointplot(
                data=df,
                x=group_column,
                y=column,
                hue=hue,
                order=specific_order,
                dodge=0.8,
                join=False,
                estimator=np.median,
                ci=None,
                markers="_",
                scale=1.5,
                color='black',
                ax=ax
            )
        else:
            sns.pointplot(
                data=df,
                x=group_column,
                y=column,
                order=specific_order,
                join=False,
                estimator=np.median,
                ci=None,
                markers="_",
                scale=1.5,
                color='black',
                ax=ax
            )

        # ----- Title -----
        ax.set_title(column)

    plt.tight_layout()
    plt.show()

#%%


def plot_by_fraction(input_df, target_cols, base_col, group_column, specific_order, p_type='violin', base_width=10, rotate=0):
    # output_path="C:/Users/leongs/OneDrive - Takara Bio USA, Inc/Desktop/test/output_figs/Alpha1_NewSoftware_SSPro_PBMC/"

    ncol = len(target_cols)
    if ncol == 1:
        fig, axs = plt.subplots()
    else:
        fig, axs = plt.subplots(ncols=ncol)

    for i, target_col in enumerate(target_cols):
        # Create dataframe by fraction
        df = pd.DataFrame(data=[], columns=['Fraction', 'Exp_Key'])
        frac_list = (input_df[target_col] / input_df[base_col]).values.tolist()
        
        #this line sets the % value
        # frac_list = [x * 100 for x in frac_list]

        data = {
            'Fraction': frac_list,
            group_column: input_df[group_column].values.tolist(), 
        }
        # df = df.append(pd.DataFrame(data), ignore_index=True)
        df = pd.concat([df,pd.DataFrame(data)], ignore_index=True)
        df[["Instrument", 'Read Depth']]= df["Sample"].str.split("_", expand = True)
        
        ax = None
        if p_type == 'violin':
            if ncol == 1:
                ax = sns.violinplot(data=df, x=group_column, y='Fraction', width=0.4, hue ='Instrument')
            else:
                ax = sns.violinplot(data=df, x=group_column, y='Fraction',  hue ='Instrument', ax=axs[i])
        elif p_type == 'box':
            if ncol == 1:
                ax = sns.boxplot(data=df, x=group_column, y='Fraction')
            else:
                ax = sns.boxplot(data=df, x=group_column, y='Fraction', ax=axs[i])

        if ax is None:
            continue
             
        ax.legend([], [], frameon=False)
        
        df["Sample"] = pd.Categorical(df["Sample"], categories=specific_order, ordered=True)
        medians = df.groupby('Sample')['Fraction'].median()
        #comment to remove values of median
        for xtick, median in enumerate(medians):
            ax.text(xtick + 0.45, median + 0.01, f"{median:.2f}", horizontalalignment='right', color='black', fontsize=12)
        # ax.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Median')])

        # ax.plot([-0.5, 4.5],[cutoff, cutoff], color='red', linewidth=1)
        ax.set_title(target_col)
        if rotate == 0:
            ax.set_xticklabels(ax.get_xticklabels())
        else:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotate, ha='right')
    
    # fig.suptitle(f'Fractions based on {base_col}')
    fig.set_figwidth(ncol * base_width)
    # fig.savefig(output_path+target_col + ' - Fraction')


#%%
def read_log(log_path,dewpointlog_path,dewpointlog_path2,instrument,experiment_start_date):

    '''This function is used to read log file, extracts date, time and step information
    Then it matches with dataframe returned from log_parse.'''
    #if "Aborted" is between "Start" and "Start" -- remove first "Start"
    # OR no "Aborted" between any "Start"
    # OR 3 rows -- Start,Chip Shift, Completed
    
    date_time_pattern = r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}.\d{3} -\d{2}:\d{2}'
    valid_words=["Cells and Controls","Scan Chip","RT","PCR","Tagmentation","P5 Index","P7 Index"] #only valid for SSPro
    valid_start_pattern = r'Started (' + '|'.join(re.escape(word) for word in valid_words) + r') step'
    valid_end_pattern = r'Completed (' + '|'.join(re.escape(word) for word in valid_words) + r') step'
    start_times=[]
    end_times=[]    
    process_names=[]

    with open(log_path, 'r') as log_file:
        in_process = False  # Flag to indicate whether we are inside a process
        process_name = None  # Stores the current process name
    
        for line in log_file:
            # Search date and time pattern in each line
            match_datetime = re.search(date_time_pattern, line)
            if match_datetime:
                date_time = match_datetime.group(0)  # Store date_time when it encounters
    
            if re.search(r'\[INF] Aborted |\[ERR\] ErrorItem|#########################', line):
                if in_process:
                    # Abort detected inside a process, remove the last started entry
                    start_times.pop()
                    process_names.pop()
                in_process = False
                process_name = None
            elif re.search(valid_end_pattern, line):
                if in_process:
                    end_times.append(date_time)
                    process_names[-1] = process_name
                in_process = False
                process_name = None
            elif re.search(valid_start_pattern, line):
                match_process_name = re.search(valid_start_pattern, line)
                if match_process_name:
                    current_process_name = match_process_name.group(1)
                    in_process = True
                    process_name = current_process_name
                    start_times.append(date_time)
                    process_names.append(current_process_name)
            elif re.search(r'\[INF\] Chip Shift', line):
                # Handle the case where "Started" is followed by "Chip Shift"
                if in_process:
                    continue
    
            # Check if we're processing a "Started" line
            elif not in_process and re.search(valid_start_pattern, line):
                start_times.append(date_time)
                process_names.append(process_name)
    
    # Ensure that the last started process is added if not completed
    if in_process and process_name:
        end_times.append(date_time)
        process_names[-1] = process_name
    start_date=pd.to_datetime(experiment_start_date,utc=True)

    log_data = pd.DataFrame({'Start Time': start_times, 'End Time': end_times, 'Step Name': process_names})
    log_data['Start Time'] = pd.to_datetime(log_data['Start Time'])
    log_data['End Time'] = pd.to_datetime(log_data['End Time'])

    log_data = log_data[log_data['Start Time'] >= start_date]

    # print(str(len(start_times)) + " " + str(len(end_times)) + " "+ str(len(process_names))) 
    
    dewpointlog1 = log_parse(dewpointlog_path,instrument).sort_values(by='Timestamp')
    
    dewpointlog2 = log_parse(dewpointlog_path2,instrument).sort_values(by='Timestamp')
    
    dewpointlog = pd.concat([dewpointlog1,dewpointlog2])
    dewpointlog['Timestamp'] = dewpointlog['Timestamp'].dt.tz_localize(None)
    log_data['Start Time'] = log_data['Start Time'].dt.tz_localize(None)
    log_data['End Time'] = log_data['End Time'].dt.tz_localize(None)


    # Merge dewpointlog with log_data based on 'Timestamp'
    merged_df = pd.merge_asof(dewpointlog, log_data, left_on='Timestamp', right_on='Start Time')

    # Filter the merged_df to keep only rows where 'Timestamp' is within 'Start Time' and 'End Time'
    merged_df = merged_df[(merged_df['Timestamp'] >= merged_df['Start Time']) & (merged_df['Timestamp'] <= merged_df['End Time'])]

    merged_df = merged_df.reset_index()
    # # Add a new column 'Step Name' to dewpointlog with values from log_data
    # dewpointlog['Step Name'] = merged_df['Step Name']

    filtered_logs = merged_df[merged_df['Step Name'].str.contains('|'.join(valid_words),case=False,na=False)]

    log_stats={}

    params = ["ChamberTempValue", "ChipTempValue","RHValue","HumidifierPower","DehumidifierPower","LeftPlateTempValue"]
    key_mapping = {"Cells and Controls":"1_Cells and Controls","Scan Chip":"2_Scan Chip","RT":"3_RT","PCR":"4_PCR1","Tagmentation":"5_Tagmentation","P5 Index":"6_P5 Index","P7 Index":"7_P7 Index"}


    for step_name in valid_words:
        log_stats[step_name]={}
        for param in params:
            mean = filtered_logs[filtered_logs['Step Name'] == step_name][param].mean()
            range_max_min = filtered_logs[filtered_logs['Step Name'] == step_name][param].max() -  filtered_logs[filtered_logs['Step Name'] == step_name][param].min()
            log_stats[step_name][param] = {'Mean': mean, 'Range': range_max_min}
            # log_stats[step_name][param] = {'Mean': mean}


    log_stats_new = {key_mapping.get(key, key): value for key, value in log_stats.items()}
    # Flatten the nested dictionary
    flat_dict = {key: flatten_dict(value) for key, value in log_stats_new.items()}
    
    # Create a DataFrame from the flattened dictionary
    df = pd.DataFrame.from_dict(flat_dict, orient='index')
    
    return df

#%%
def flatten_dict(d, parent_key='', sep='_'):
    items = {}
    for key, value in d.items():
        new_key = f'{parent_key}{sep}{key}' if parent_key else key
        if isinstance(value, dict):
            items.update(flatten_dict(value, new_key, sep=sep))
        else:
            items[new_key] = value
    return items

#%%
'''
This function is used to parse dewpointlog for RH and temp,

'''
def log_parse(dewpointlog_path,instrument):
    colnames = ["Timestamp","Firmware time",
                "ChamberTempSetPoint",
                "ChamberTempValue",
          "ChamberTempPower",
          "HumidifierSetPoint",
          "HumidifierValue",
          "HumidifierPower",
          "ChipTempSetPoint",
          "ChipTempValue",
          "ChipTempPower",
          "LeftPlateTempSetPoint",
          "LeftPlateTempValue",
          "LeftPlateTempPower",
          "RightPlateTempSetPoint",
          "RightPlateTempValue",
          "RightPlateTempPower",
          "DehumidifierSetPoint",
          "DehumidifierValue",
          "DehumidifierPower","NA"]
    dewpointlog = pd.read_csv(dewpointlog_path,delimiter='\t',header=None)
    dewpointlog.columns=colnames
    dewpointlog['Instrument']=instrument
    dewpointlog.rename(columns={'HumidifierValue':'RHValue'},inplace=True) #HumidifierValue is same as DehumidifierValue is same as RH
    # dewpointlog.reset_index(inplace=True)
    dewpointlog['Timestamp'] = pd.to_datetime(dewpointlog['Timestamp'])
    return dewpointlog    
    
def create_scatterplot(join_df,x,y,title):
    plt.figure(figsize=(8, 6))
    plt.scatter(join_df[x],join_df[y],s=20)
    plt.xlabel(x)
    plt.ylabel(y)    
    plt.title(title)
    plt.show()
    
