# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 15:14:10 2023

@author: leongs
"""

import pandas as pd
pd.options.display.float_format = "{:,.4f}".format
import matplotlib.pyplot as plt
import seaborn as sns
import SSPro_func as func
sns.set_theme(style="whitegrid")
# import dataframe_image as dfi
import numpy as np

html_strings = []

#after unzipping the zip file from Mappa
# folder=r"C:/Users/leongs/OneDrive - Takara Bio USA, Inc/XuanRaymondExp/Shasta/SSPro_Beta2_CX1070_Modified_Lplot_20250130/For Sam/"
folder=r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"


path0=folder+"241209_SSPro_Beta2" #shasta raw
# path1=folder+r"539_beta2_human_100k" #cx raw
path2=folder+"250325_SSPro_25nLDoubleDisp_50K" #shasta 100K
# path3=folder+r"541_beta2_mouse_100k" #cx 100K

title = "250325 25nL dispense" #change

shasta = 'Beta-2'
cx = 'cx-1018'

# specific_order = ['Swap_Raw','CX1018_Raw','Swap_100K','CX1018_100K'] #beta2
# specific_order = [shasta+'_Raw',cx+'_Raw',shasta+'_100K',cx+'_100K']



specific_order = [shasta+'_Raw',shasta+'_50K']
 

sample_name=['Sample','sample'] #Shasta, CX




'''should not need anymore user input after here'''

df0 = pd.read_csv(path0+'/analysis_stats.csv')
# df1 = pd.read_csv(path1+'/analysis_stats.csv')
df2 = pd.read_csv(path2+'/analysis_stats.csv')
# df3 = pd.read_csv(path3+'/analysis_stats.csv')

#%% SSPRO spec 1


# df0['Sample'].replace({'Neg_Ctrl': 'Temp_Label', 'Pos_Ctrl': 'Neg_Ctrl'}, inplace=True)
# df0['Sample'].replace({'Temp_Label': 'Pos_Ctrl'}, inplace=True)


sspro_spec1_df = func.sspro_spec1_defunct(df0)[0] #shasta
# sspro_spec1_df = func.sspro_spec1_defunct(df1)[0] #CX
func.plot_table_from_df(sspro_spec1_df, 'Shasta',base_width=3)


#%%
# paths = [path0,path1,path2,path3]
# dfs=[df0,df1,df2,df3]

paths = [path0,path2]
dfs=[df0,df2]

for df in dfs:
    df = df.rename(columns={'Sample': 'Well_Type'},inplace=True)
well_lists = []
for path in paths:
    well_list = func.search_files_in_directory(path)
    func.cleanup_barcode(well_list)
    well_list['Source']=path
    well_lists.append(well_list)
    
for df in well_lists:  
    df = df.rename(columns={'Sample': 'Well_Type'})

well_lists[0]['Sample']=df0['Sample']=specific_order[0]
# well_lists[1]['Sample']=df1['Sample']=specific_order[1]
# well_lists[2]['Sample']=df2['Sample']=specific_order[2]
# well_lists[3]['Sample']=df3['Sample']=specific_order[3]


well_lists[1]['Sample']=df2['Sample']=specific_order[1]


 #%% SSPRO spec 2: Passing Criteria: The number of reads between samples (PBMCs): 85% of cells should be in the range of an order of magnitude.
# Moved this below sectin "XUAN_CODE" after merging dataframes into df on 2024/11/13. Does not allow user to remove specific row.

beta2_spec2_original=df0[(df0["Well_Type"] == "Sample") | (df0["Well_Type"] == "sample") | (df0["Well_Type"] == "PBMC")] #Beta2
# cx_spec2=df1[(df1["Well_Type"] == "PBMC")| (df1["Well_Type"] == "sample") ] #CX1018

sspro_spec2_df = func.sspro_spec2(beta2_spec2_original,'Barcoded_Reads')
# sspro_spec2_df = func.sspro_spec2(cx_spec2,'Barcoded_Reads')
func.plot_table_from_df(sspro_spec2_df[0],'no index',base_width=2)


#%% 

# Merge well_lists and dfs based on 'Barcode' and 'Sample' columns
list_of_df=[]
for i in range(len(dfs)):
   
    # df_test = pd.merge(dfs[i], well_lists[i][['Barcode', 'Row', 'Col', 'Sample','Signal1']], on=['Barcode', 'Sample'])
    df_test = pd.merge(dfs[i], well_lists[i][['Barcode', 'Row', 'Col', 'Sample','Signal1']], on=['Barcode', 'Sample']) #one off

    list_of_df.append(df_test)
df = pd.concat(list_of_df, ignore_index=True)

# now drop unnecessary samples
for df in list_of_df: 
    df.drop(df[df['Well_Type'] == 'Non_sample'].index, inplace=True)              # type: ignore
    df.drop(df[df['Well_Type'] == 'Pos_Ctrl'].index, inplace=True)                # type: ignore
    df.drop(df[df['Well_Type'] == 'Neg_Ctrl'].index, inplace=True)

for df in list_of_df:     
    df['%Intergenic_Reads']=df['Intergenic_Reads']*100/df['Barcoded_Reads'] #CHECK IF THIS IS CORREC
    df['%Ribosomal']=df['Ribosomal_Reads']*100/df['Barcoded_Reads']
    df['%Mitochondrial']=df['Mitochondrial_Reads']*100/df['Barcoded_Reads']
    df['Total_Exon_Reads'] = df['Exon_Reads'] + df['Ambiguous_Exon_Reads']
    df['Total_Intron_Reads'] = df['Intron_Reads'] + df['Ambiguous_Intron_Reads']
    df['i7'] = df['Barcode'].str[:8]
    df['i5'] = df['Barcode'].str[-8:]
    
    
    

#%% sspro spec4: : % unmapped for PBMCs should be <10%. Controls are excluded. 

sspro_spec4_df = func.sspro_spec4(list_of_df[0][list_of_df[0]["Well_Type"] == sample_name[0]])  #beta1 raw
func.plot_table_from_df(sspro_spec4_df, 'no index', base_width=5)

sspro_spec4_df = func.sspro_spec4(list_of_df[1][list_of_df[1]["Well_Type"] == sample_name[1]]) # CX raw
func.plot_table_from_df(sspro_spec4_df, 'no index', base_width=5)



#%% XUAN CODE
# collapse list_of_df
df = pd.concat(list_of_df, ignore_index=True) #use this with Xuan code

df.drop(df[df['Sample'] == 'Non_sample'].index, inplace=True)                     # type: ignore
df.drop(df[df['Sample'] == 'Pos_Ctrl'].index, inplace=True)                       # type: ignore
df.drop(df[df['Sample'] == 'Neg_Ctrl'].index, inplace=True)                       # type: ignore


df.drop(df[df['Well_Type'] == 'Non_sample'].index, inplace=True)                     # type: ignore
df.drop(df[df['Well_Type'] == 'Pos_Ctrl'].index, inplace=True)                       # type: ignore
df.drop(df[df['Well_Type'] == 'Neg_Ctrl'].index, inplace=True)   
df[['Instrument', 'Read Depth']] = df.Sample.str.split("_", expand = True)

#df = df.sort_values('Barcoded_Reads')
df['Barcoded_Reads'] = df['Barcoded_Reads'].astype(float)
df['Mapped_Reads']= df['Mapped_Reads'].astype(float)

#%%
# Side by side spec2: %Barcoded Reads
# # fraction of barcoded reads = barcoded reads/demux reads (excluding ctrls, including samples and nonsamples)
# barcoded_sum=df.groupby(['Sample'])['Barcoded_Reads'].sum().sort_index()
# barcoded_sum/demux_reads_sum['Sum']


# df = df.merge(total_reads, on='Sample', suffixes=('', '_total_reads'))
# df['Fraction_Barcoded_Reads'] = df['Barcoded_Reads']/df['Total_Reads']
# df.groupby(['Sample'])['Total_Reads'].mean()


# side by side comparison spec 4 
# df['Total_Exon_Reads'] = df['Exon_Reads'] + df['Ambiguous_Exon_Reads']
# df['%Exon_Reads'] = (df['Exon_Reads'] + df['Ambiguous_Exon_Reads'])/df['Barcoded_Reads']
# df.groupby(['Sample'])['%Exon_Reads'].mean().sort_index(ascending=False)
# df.groupby(['Sample'])['%Exon_Reads'].median().sort_index(ascending=False)


# # Insert additional column
# df['Total_Intron_Reads'] = df['Intron_Reads'] + df['Ambiguous_Intron_Reads']
# df['Group'] = df['Sample']
# #df = df.sort_values('Sample')

# df.groupby(['Sample'])['No_of_Genes'].mean()
# df.groupby(['Sample'])['No_of_Genes'].median()
# # df['ExonFraction'] = df['Exon_Reads']/df['Trimmed_Reads']
# # df.groupby(['Sample'])['ExonFraction'].median()

#%%
#use this one for side by side comparison spec 6
df['Ribosomal+Mitochondrial_Reads'] = df['Ribosomal_Reads']+df['Mitochondrial_Reads']
df['Ribosomal+Mitochondrial_Fraction'] = df['Ribosomal+Mitochondrial_Reads']/df['Barcoded_Reads']
df.groupby(['Sample'])['Ribosomal+Mitochondrial_Fraction'].mean().sort_index(ascending=False)
df.groupby(['Sample'])['Ribosomal+Mitochondrial_Fraction'].median().sort_index(ascending=False)


df['%Mitochondrial']=df['Mitochondrial_Reads']*100/df['Barcoded_Reads']
# df['%Ribosomal']=df['Ribosomal_Reads']*100/df['Barcoded_Reads']
df['%Intergenic_Reads']=df['Intergenic_Reads']*100/df['Barcoded_Reads']


#%%
# side by side comparison spec 4 
# df['ExonFraction'] = (df['Exon_Reads'] + df['Ambiguous_Exon_Reads'])/df['Barcoded_Reads']
# df.groupby(['Sample'])['ExonFraction'].median().sort_index(ascending=False)
# df.groupby(['Sample'])['ExonFraction'].mean().sort_index(ascending=False)



#%% # side by side comparison spec 4 and 5
#group then sum USE THIS ONE

exonReads = df.groupby('Sample').apply(
    lambda group: (group['Exon_Reads'].sum() + group['Ambiguous_Exon_Reads'].sum()) / group['Barcoded_Reads'].sum()
)*100 #Percent
print(exonReads.sort_index(ascending=False).round(1))

# Calculate Ribosomal + Mitochondrial Reads and Ribosomal + Mitochondrial Fraction per group
riboMitoSum = df.groupby('Sample').apply(
    lambda group: (group['Ribosomal_Reads'].sum() + group['Mitochondrial_Reads'].sum()) / group['Barcoded_Reads'].sum()
)*100
print(riboMitoSum.sort_index(ascending=False).round(1))



#%% SSPRO spec 2: Passing Criteria: The number of reads between samples (PBMCs): 85% of cells should be in the range of an order of magnitude.
# Test function after removing specific failed rows and columns

beta2_spec2 = df[((df["Well_Type"] == "Sample") | 
                  (df["Well_Type"] == "sample") | 
                  (df["Well_Type"] == "PBMC")) & 
                 (df['Group'] == shasta+'_Raw')]
cx_spec2=df[((df["Well_Type"] == "PBMC")| (df["Well_Type"] == "sample")) & (df['Group']==cx+'_Raw')] #CX1018

sspro_spec2_df = func.sspro_spec2(beta2_spec2,'Barcoded_Reads')
# sspro_spec2_df = func.sspro_spec2(cx_spec2,'Barcoded_Reads')
func.plot_table_from_df(sspro_spec2_df[0],'no index',base_width=2)


#%%
#USING XUAN df for sspro spec 3 and side by side comparison spec 1, dropped all non-samples
spec3_list =[df[df['Sample']==specific_order[0]],df[df['Sample']==specific_order[1]]] #beta1 raw and cx raw
sspro_spec3_df = func.sspro_spec3(spec3_list,'%Intergenic_Reads')
func.plot_table_from_df(sspro_spec3_df[0],'no index',base_width=2.1)


#sspro spec 3 for 100K reads
spec3_list =[df[df['Sample']==specific_order[2]],df[df['Sample']==specific_order[3]]] #beta1 100k and cx 100k
sspro_spec3_df = func.sspro_spec3(spec3_list,'%Intergenic_Reads')
func.plot_table_from_df(sspro_spec3_df[0],'no index',base_width=2.1)

# spec3_list =[df[df['Sample']==specific_order[0]]] #beta1 raw only
# sspro_spec3_df = func.sspro_spec3(spec3_list,'%Intergenic_Reads')
# func.plot_table_from_df(sspro_spec3_df[0],'no index',base_width=2.5)


# NOTE: MUST BE IN THIS ORDER: df, shasta,shasta_100K
sscomp_spec1_df = func.sscomparison_spec1_nocx(df,specific_order[0],specific_order[1])


#%%    
# posCtrl = df[df['Well_Type']=='Pos_Ctrl']
# list_of_df = [group for _, group in df.groupby('Sample')]
# test=posCtrl[posCtrl['Sample'].isin(['CX1018_Raw','Beta1_Raw'])]
# test=df[df['Sample'].isin(specific_order[:2])]



test=df[df['Sample'].isin([specific_order[0]])]

# func.plot(test, ['Barcoded_Reads'], 'Sample',specific_order[:-2], base_width=5)
func.plot(test, ['Barcoded_Reads'],'Sample',[specific_order[0]], base_width=3.5, base_height=4.5)

test=df[df['Sample'].isin(specific_order[2::])]



func.plot(df, ['No_of_Genes'], 'Sample',specific_order, base_width=6,base_height=4)

test=df[df['Sample'].isin([specific_order[2],specific_order[3]])]

func.plot(df, ['No_of_Genes'], 'Sample',[specific_order[2],specific_order[3]], base_width=4,base_height=4)


func.plot_by_fraction(df, ['Mapped_Reads'], 'Barcoded_Reads', 'Sample', specific_order, base_width=4)
func.plot_by_fraction(df, ['Uniquely_Mapped_Reads'], 'Mapped_Reads', 'Sample',specific_order, base_width=10)
func.plot_by_fraction(df, ['Total_Exon_Reads', 'Total_Intron_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=5)
func.plot_by_fraction(df, ['Ribosomal+Mitochondrial_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=5)
func.plot_by_fraction(df, ['Total_Exon_Reads'], 'Barcoded_Reads', 'Sample', specific_order,base_width=10)
func.plot_by_fraction(df, ['Mitochondrial_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=6)
func.plot_by_fraction(df, ['Ribosomal_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=6)
func.plot_by_fraction(df, ['Intergenic_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=6)

func.plot_by_fraction(df, ['Ribosomal+Mitochondrial_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=6)


#%%

func.gen_heatmap_individual(list_of_df,'Barcoded_Reads','Sample','log', base_width=20)

func.gen_heatmap(list_of_df,'Barcoded_Reads','Sample','log', base_width=10)
func.gen_heatmap_individual(list_of_df,'Barcoded_Reads','Sample','autumn', base_width=10)

func.gen_heatmap(list_of_df,'No_of_Genes','Sample','log', base_width=7)
func.gen_heatmap_individual(list_of_df,'No_of_Genes','Sample','log', base_width=7)
func.gen_heatmap(list_of_df,'%Intergenic_Reads','Sample','bwr',base_width=7)
func.gen_heatmap(list_of_df,'%Ribosomal','Sample','bwr',base_width=7)
func.gen_heatmap(list_of_df,'%Mitochondrial','Sample','bwr',base_width=7)
func.plot_barchart(list_of_df,'Barcoded_Reads','Row','Sample',base_width=5)

func.gen_heatmap(list_of_df,'Total_Intron_Reads','Sample','log', base_width=7)


#%% Shasta RH logs
import SSPro_func as func
experiment_start_date='2024-04-03'
master_path = "S:/XuanLi/Shasta_SSPro/147279C_20240403102500_SSPro_Beta2_CellRTSwap/"
log_path = master_path + "log-202404.txt"
dewpointlog_path= master_path + "DewPointLog-20240403.tsv"
dewpointlog_path2= master_path + "DewPointLog-20240404.tsv"

instrument="Beta2"
# test=func.log_parse(dewpointlog_path,instrument)
beta2=func.read_log(log_path,dewpointlog_path,dewpointlog_path2,instrument,experiment_start_date)
func.plot_table_from_df(beta2, instrument, base_width=10)
# dfi.export(beta2,instrument+'.png')


#%% CX RH logs 
import SSPro_func as func
master_path = "S:/XuanLi/Shasta_SSPro/2024.02.12.10.58-136271_CX1018_SSPro_Test6/Logs/"
log_path = master_path + "Debug_0212.log" #this is for instrument logs 
start_date = "2024.02.12"  # Replace with your desired start date in this format 2023.10.09
input_path = master_path+"20240212_CX1018_teraterm.log" #this is for teraterm logs
output_folder=master_path+"/test"
instrument = "0212_CX1018"
df=func.read_log_CX(log_path, input_path, output_folder, instrument, start_date)
func.plot_table_from_df(df, instrument, base_width=6)
dfi.export(df,master_path+"RH log_"+instrument+'.png')


#%% i7 and i5 reads distribution for shasta
# i7_list = list_of_df[0].groupby(['i7'])['Barcoded_Reads'].mean().reset_index()
# barcodes = list_of_df[0][['i7','i5','Row','Col','Barcoded_Reads']]
# # AGGCCAAG = barcodes[barcodes['i5']=="AGGCCAAG"]
# # func.plot_table_from_df(AGGCCAAG, 'no index', base_width=3)

# # i7_list = list_of_df[0][list_of_df[0]["Well_Type"]=="PBMC"].groupby('i7')['Barcoded_Reads'].mean().reset_index()

# plt.figure(figsize=(15, 5))  # Width: 8 inches, Height: 6 inches
# plt.bar(i7_list['i7'],i7_list['Barcoded_Reads'])
# plt.xlabel('i7')
# plt.xticks(rotation=90)
# plt.ylabel('Mean Value of Barcoded_Reads')
# plt.title('i7 distribution per barcode',fontsize=16)
# plt.show()

# i7_list = i7_list.sort_values(by='Barcoded_Reads',ascending=False)
# func.plot_table_from_df(pd.concat([i7_list.head(),i7_list.tail()]), 'no index', base_width=2)


# # i5_list = list_of_df[0][list_of_df[0]["Well_Type"]=="PBMC"].groupby('i5')['Barcoded_Reads'].mean().reset_index()
# i5_list = list_of_df[0].groupby('i5')['Barcoded_Reads'].mean().reset_index()

# plt.figure(figsize=(15, 5))  # Width: 8 inches, Height: 6 inches
# plt.bar(i5_list['i5'],i5_list['Barcoded_Reads'])
# plt.xlabel('i5')
# plt.xticks(rotation=90)
# plt.ylabel('Mean Value of Barcoded_Reads')
# plt.title('i5 distribution per barcode',fontsize=16)
# plt.show()

# i5_list = i5_list.sort_values(by='Barcoded_Reads',ascending=False)
# func.plot_table_from_df(pd.concat([i5_list.head(),i5_list.tail()]), 'no index', base_width=2)

 #%% signal1 vs reads

# beta1 = df[df['Sample'] == 'Beta1_Raw']

# # Determine the maximum value in 'Barcoded_Reads' column
# max_barcoded_reads = beta1['Barcoded_Reads'].max()

# # Bin 'Barcoded_Reads' column using NumPy with bin size of 2000
# bins = np.arange(0, int(max_barcoded_reads) + 2001, 2000)
# beta2test3['bin'] = pd.cut(beta2test3['Barcoded_Reads'], bins=bins, right=False)

# # Calculate the mean value of 'Signal1' for each bin
# mean_signal1_by_bin = beta2test3.groupby('bin')['Signal1'].mean()

# # Plot as a bar chart
# plt.figure(figsize=(10, 6))
# bars = plt.bar(range(1, len(mean_signal1_by_bin) + 1), mean_signal1_by_bin, color='steelblue', edgecolor='none', width=0.8)  # Adjust width if needed
# # for bar, bin_number in zip(bars, range(1, len(mean_signal1_by_bin) + 1)):
# #     plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{bin_number}', ha='center', va='bottom')
# plt.xlabel('Bins')
# plt.ylabel('Signal1 Mean')
# plt.title('Bar Chart of Signal1 by Barcoded Reads Bins')
# plt.xticks(range(1, len(mean_signal1_by_bin) + 1), rotation=45)
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# # Plot as boxplot
# plt.figure(figsize=(10, 6))
# boxprops = dict(color='steelblue', linewidth=1.5)
# medianprops = dict(color='orange', linewidth=2)
# plt.boxplot([beta2test3[beta2test3['bin'] == b]['Signal1'] for b in mean_signal1_by_bin.index], labels=range(1,len(mean_signal1_by_bin.index)+1), boxprops=boxprops, medianprops=medianprops)
# plt.xlabel('Bins')
# plt.ylabel('Signal1')
# plt.title('Box Plot of Signal1 by Barcoded Reads Bins')
# plt.xticks(rotation=45)
# plt.grid(True)
# plt.tight_layout()
# plt.show()
# plt.show()