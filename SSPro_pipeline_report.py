# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 15:14:10 2023

@author: leongs
"""
import seaborn as sns

import pandas as pd
pd.options.display.float_format = "{:,.3f}".format
import matplotlib.pyplot as plt
import SSPro_func as func
import SSPro_report_func as reportfunc
sns.set_theme(style="whitegrid")
# import dataframe_image as dfi
import numpy as np
import plotly.graph_objs as go
from matplotlib.colors import ListedColormap, Normalize, LogNorm
import os
import re
#%%

html_strings = []

#after unzipping the zip file from Mappa
folder=r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"
folder1=r"C:/Users/leongs/OneDrive - Takara Bio USA, Inc/0. Projects/Shasta instrument/Hammerhead syringe/SSPro double hammerhead/"

path0=folder1+r"SSPro 20nL stacked validation hammerhead_noDownsampling" #v2 hammer raw
path1=folder+"241209_SSPro_Beta2" #v1 raw
path2=folder1+r"SSPro 20nL stacked validation hammerhead_downsampling50K" #v2 hammer 50K
path3=folder+"241209_SSPro_Beta2_50K" #v1 50K

title = "SSPro hammerhead validation 2025-12-09" #change

masterpath = folder1 + "SSPro_Reports/"
shasta = 'v2'
cx = 'v1'

raw_path = title
notes = "Downsampled to 50k"
cx_total_cells = 1
cx_pass = 1
shasta_total_cells = 1
shasta_pass = 1

# cx_pass_02 = 1
# shasta_pass_02 = 1

gene_sensitivity=0.8 #used to be 0.9 for cx vs shasta

date = title.split("_")[0]
filename = title + ".html"
savepath = masterpath + title + "/"+ filename
image_folder = masterpath+ title+"/images/"
os.makedirs(image_folder, exist_ok=True)




#%%
specific_order = [shasta+'_Raw',cx+'_Raw',shasta+'_50K',cx+'_50K'] #v2 raw, v1 raw, v2 50k, v1 50k


'''should not need anymore user input after here'''
df0 = pd.read_csv(path0+'/analysis_stats.csv')
df1 = pd.read_csv(path1+'/analysis_stats.csv')
df2 = pd.read_csv(path2+'/analysis_stats.csv')
df3 = pd.read_csv(path3+'/analysis_stats.csv')

demux_df0 = pd.read_csv(path0+'/demux_counts_all.csv',names=['Barcode','Well_Type','Reads'])
demux_df1 = pd.read_csv(path1+'/demux_counts_all.csv',names=['Barcode','Well_Type','Reads'])
demux_df2 = pd.read_csv(path2+'/demux_counts_all.csv',names=['Barcode','Well_Type','Reads'])
demux_df3 = pd.read_csv(path3+'/demux_counts_all.csv',names=['Barcode','Well_Type','Reads'])

#demux: this drops controls but keeps non samples
[dataframe.drop(dataframe[dataframe['Well_Type'].isin(['Neg_Ctrl', 'Pos_Ctrl'])].index, inplace=True) for dataframe in [demux_df0, demux_df1, demux_df2, demux_df3]]

demux_df0['Sample'] = specific_order[0]
demux_df1['Sample'] = specific_order[1]
demux_df2['Sample'] = specific_order[2]
demux_df3['Sample'] = specific_order[3]

demux_list_of_df=[demux_df0,demux_df1,demux_df2,demux_df3]
demux_reads_sum = [df['Reads'].sum() for df in [demux_df0, demux_df1, demux_df2, demux_df3]]
demux_reads_sum = pd.DataFrame({'Sum': demux_reads_sum, 'Sample': specific_order})

demux_reads_sum=demux_reads_sum.set_index('Sample')

# df0['Sample'].replace({'Neg_Ctrl': 'Temp_Label', 'Pos_Ctrl': 'Neg_Ctrl'}, inplace=True)
# df0['Sample'].replace({'Temp_Label': 'Pos_Ctrl'}, inplace=True)

# SSPro Spec 1
sspro_spec1_df_shasta = func.sspro_spec1_defunct(df0)
sspro_spec1_df_CX = func.sspro_spec1_defunct(df1)



paths = [path0,path1,path2,path3]
dfs=[df0,df1,df2,df3]

# paths = [path0,path1]
# dfs=[df0,df1]

for df in dfs:
    df = df.rename(columns={'Sample': 'Well_Type'},inplace=True)
well_lists = []
for path in paths:
    well_list = func.search_files_in_directory(path)
    func.cleanup_barcode(well_list)
    # well_list['Source']=path
    well_lists.append(well_list)
    
for df in well_lists:  
    df = df.rename(columns={'Sample': 'Well_Type'})

well_lists[0]['Sample']=df0['Sample']=specific_order[0]
well_lists[1]['Sample']=df1['Sample']=specific_order[1]
well_lists[2]['Sample']=df2['Sample']=specific_order[2]
well_lists[3]['Sample']=df3['Sample']=specific_order[3]

# SSPRO spec 2: Passing Criteria: The number of reads between samples (PBMCs): 85% of cells should be in the range of an order of magnitude. &#8804;

sspro_spec2_shasta=df0[(df0["Well_Type"] == "Sample") | (df0["Well_Type"] == "PBMC_HMN621887") | (df0["Well_Type"] == "PBMC") | (df0["Well_Type"] == "sample")] #v2
sspro_spec2_CX=df1[(df1["Well_Type"] == "sample") | (df1["Well_Type"] == "Sample")|  (df1["Well_Type"] == "PBMC") | (df1["Well_Type"] == "HMN621887")] #v1


sspro_spec2_df_shasta = func.sspro_spec2(sspro_spec2_shasta,'Barcoded_Reads')
sspro_spec2_df_shasta[0].set_index('Quantile',inplace=True)
sspro_spec2_df_CX = func.sspro_spec2(sspro_spec2_CX,'Barcoded_Reads')
sspro_spec2_df_CX[0].set_index('Quantile',inplace=True)

if sspro_spec2_df_shasta[3] <= 10: 
    spec2_shasta = "Pass"
else: spec2_shasta = "Fail"


if sspro_spec2_df_CX[3] <= 10: 
    spec2_cx = "Pass"
else: spec2_cx = "Fail"

sspro_spec2_df = pd.concat([sspro_spec2_df_shasta[0],sspro_spec2_df_CX[0]],axis=1)
sspro_spec2_df.columns = [shasta+' Barcoded_Reads', cx+' Barcoded_Reads']



# Merge well_lists and dfs based on 'Barcode' and 'Sample' columns
list_of_df=[]
for i in range(len(dfs)):
   
    df_test = pd.merge(dfs[i], well_lists[i][['Barcode', 'Row', 'Col', 'Sample']], on=['Barcode', 'Sample'])
    list_of_df.append(df_test)
df = pd.concat(list_of_df, ignore_index=True)

# Side by side spec2: %Barcoded Reads
# # fraction of barcoded reads = barcoded reads/demux reads (excluding ctrls, including samples and nonsamples)
barcoded_sum=df.groupby(['Sample'])['Barcoded_Reads'].sum().sort_index()
sscomp_spec2 = barcoded_sum / demux_reads_sum['Sum']
sscomp_spec2 = sscomp_spec2[sscomp_spec2.index.str.endswith('_Raw')]
sscomp_spec2 = pd.DataFrame(sscomp_spec2).T
sscomp_spec2[f"{shasta} to {cx} ratio"] = sscomp_spec2[f"{shasta}_Raw"]/sscomp_spec2[f"{cx}_Raw"]
sscomp_spec2['Spec'] = sscomp_spec2[f"{shasta} to {cx} ratio"].apply(lambda x: "Fail" if x < 0.9 else "Pass")

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

# collapse list_of_df
df = pd.concat(list_of_df, ignore_index=True) #use this with Xuan code

#sspro spec4: : % unmapped for PBMCs should be <10%. Controls are excluded. 

# sspro_spec4_df_shasta = func.sspro_spec4(list_of_df[0][list_of_df[0]["Well_Type"] == specific_order[0]])  #beta1 raw
# sspro_spec4_df_cx = func.sspro_spec4(list_of_df[1][list_of_df[1]["Well_Type"] == specific_order[1]]) # CX raw

sspro_spec4_df_shasta = func.sspro_spec4(list_of_df[0][list_of_df[0]["Sample"] == f"{shasta}_Raw"])  #beta1 raw
sspro_spec4_df_cx = func.sspro_spec4(list_of_df[1][list_of_df[1]["Sample"] == f"{cx}_Raw"]) # CX raw


sspro_spec4 = pd.concat([sspro_spec4_df_shasta,sspro_spec4_df_cx],keys=[shasta,cx])
sspro_spec4['Spec'] = sspro_spec4['% Unmapped Reads'].apply(lambda x: "Fail" if x > 10 else "Pass")

sscomp_spec3 = pd.DataFrame(columns=[cx, shasta, f'{shasta} to {cx} Ratio'], 
                            index=['Total no of cells', 'Passed Cogent QC', '% Passing Cogent QC'])
sscomp_spec3.loc['Total no of cells',cx] = cx_total_cells #change
sscomp_spec3.loc['Passed Cogent QC',cx] = cx_pass #change
sscomp_spec3.loc['% Passing Cogent QC',cx] = sscomp_spec3.loc['Passed Cogent QC',cx]/sscomp_spec3.loc['Total no of cells',cx]


sscomp_spec3.loc['Total no of cells',shasta] = shasta_total_cells #change
sscomp_spec3.loc['Passed Cogent QC',shasta] = shasta_pass #change
sscomp_spec3.loc['% Passing Cogent QC',shasta] = sscomp_spec3.loc['Passed Cogent QC',shasta]/sscomp_spec3.loc['Total no of cells',shasta] 


sscomp_spec3.loc['% Passing Cogent QC',f'{shasta} to {cx} Ratio'] = sscomp_spec3.loc['% Passing Cogent QC',shasta]/sscomp_spec3.loc['% Passing Cogent QC',cx]

sscomp_spec3['Spec'] = sscomp_spec3[f'{shasta} to {cx} Ratio'].apply(lambda x: '' if pd.isna(x) else ("Pass" if x > 0.9 else "Fail"))
 
#0.2 intergenic reads
# sscomp_spec3_02 = pd.DataFrame(columns=['CX1018', 'Shasta', 'Shasta to CX-1018 Ratio (20% intergenic reads)'], 
#                             index=['Total no of cells', 'Passed Cogent QC', '% Passing Cogent QC'])
# sscomp_spec3_02.loc['Total no of cells','CX1018'] = cx_total_cells #change
# sscomp_spec3_02.loc['Passed Cogent QC','CX1018'] = cx_pass_02 #change
# sscomp_spec3_02.loc['% Passing Cogent QC','CX1018'] = sscomp_spec3_02.loc['Passed Cogent QC','CX1018']/sscomp_spec3_02.loc['Total no of cells','CX1018']


# sscomp_spec3_02.loc['Total no of cells','Shasta'] = shasta_total_cells #change
# sscomp_spec3_02.loc['Passed Cogent QC','Shasta'] = shasta_pass_02 #change
# sscomp_spec3_02.loc['% Passing Cogent QC','Shasta'] = sscomp_spec3_02.loc['Passed Cogent QC','Shasta']/sscomp_spec3_02.loc['Total no of cells','Shasta'] 


# sscomp_spec3_02.loc['% Passing Cogent QC','Shasta to CX-1018 Ratio (20% intergenic reads)'] = sscomp_spec3_02.loc['% Passing Cogent QC','Shasta']/sscomp_spec3_02.loc['% Passing Cogent QC','CX1018']

# sscomp_spec3_02['Spec'] = sscomp_spec3_02['Shasta to CX-1018 Ratio (20% intergenic reads)'].apply(lambda x: '' if pd.isna(x) else ("Pass" if x > 0.9 else "Fail"))
 

# XUAN CODE

df.drop(df[df['Sample'] == 'Non_sample'].index, inplace=True)                     # type: ignore
df.drop(df[df['Sample'] == 'Pos_Ctrl'].index, inplace=True)                       # type: ignore
df.drop(df[df['Sample'] == 'Neg_Ctrl'].index, inplace=True)                       # type: ignore
df[['Instrument', 'Read Depth']] = df.Sample.str.split("_", expand = True)

#df = df.sort_values('Barcoded_Reads')
df['Barcoded_Reads'] = df['Barcoded_Reads'].astype(float)
df['Mapped_Reads']= df['Mapped_Reads'].astype(float)


# Side by side spec2: %Barcoded Reads
# fraction of barcoded reads = barcoded reads/demux reads (excluding ctrls, including samples and nonsamples)
barcoded_sum=df.groupby(['Sample'])['Barcoded_Reads'].sum().sort_index()

# df['Total_Exon_Reads'] = df['Exon_Reads'] + df['Ambiguous_Exon_Reads']
# df['%Exon_Reads'] = (df['Exon_Reads'] + df['Ambiguous_Exon_Reads'])/df['Barcoded_Reads']

sscomp_spec4_df = df.groupby('Sample').apply(
    lambda group: (group['Exon_Reads'].sum() + group['Ambiguous_Exon_Reads'].sum()) / group['Barcoded_Reads'].sum()
).to_frame(name='%Exon_Reads')

sscomp_spec4_df["Reads"] = [sample.split('_')[-1] for sample in sscomp_spec4_df.index]
sscomp_spec4_df["SampleName"] = [sample.split('_')[0] for sample in sscomp_spec4_df.index]
sscomp_spec4 = sscomp_spec4_df.pivot_table(index='Reads', columns='SampleName', values='%Exon_Reads')
sscomp_spec4[f"{shasta} to {cx} Ratio"] = sscomp_spec4[shasta]/sscomp_spec4[cx] #hardcoded to be col0/col1
sscomp_spec4["Spec"] = sscomp_spec4[f"{shasta} to {cx} Ratio"].apply(lambda x: "Fail" if x < gene_sensitivity else "Pass")

# Insert additional column
df['Total_Intron_Reads'] = df['Intron_Reads'] + df['Ambiguous_Intron_Reads']
df['Group'] = df['Sample']
#df = df.sort_values('Sample')

df.groupby(['Sample'])['No_of_Genes'].mean()
df.groupby(['Sample'])['No_of_Genes'].median()
# df['ExonFraction'] = df['Exon_Reads']/df['Trimmed_Reads']
# df.groupby(['Sample'])['ExonFraction'].median()
df['Ribosomal+Mitochondrial_Reads'] = df['Ribosomal_Reads']+df['Mitochondrial_Reads']
df['Ribosomal+Mitochondrial_Fraction'] = df['Ribosomal+Mitochondrial_Reads']/df['Trimmed_Reads']
#%%
#use this one for side by side comparison spec 5
sscomp_spec5_df = df.groupby('Sample').apply(
    lambda group: (group['Ribosomal_Reads'].sum() + group['Mitochondrial_Reads'].sum()) / group['Barcoded_Reads'].sum()
).to_frame(name='Ribosomal+Mitochondrial_Fraction')
sscomp_spec5_df["Reads"] = [sample.split('_')[-1] for sample in sscomp_spec5_df.index]
sscomp_spec5_df["SampleName"] = [sample.split('_')[0] for sample in sscomp_spec5_df.index]
sscomp_spec5 = sscomp_spec5_df.pivot_table(index='Reads', columns='SampleName',values='Ribosomal+Mitochondrial_Fraction')
sscomp_spec5[f"{shasta} to {cx} Ratio"] = sscomp_spec5[shasta]/sscomp_spec5[cx]
sscomp_spec5["Spec"] = sscomp_spec5[f"{shasta} to {cx} Ratio"].apply(lambda x: '' if pd.isna(x) else ("Fail" if x > 1.1 else "Pass"))


df['ExonFraction'] = (df['Exon_Reads'] + df['Ambiguous_Exon_Reads'])/df['Barcoded_Reads']
df.groupby(['Sample'])['ExonFraction'].median()
df.groupby(['Sample'])['ExonFraction'].mean()

df.groupby(['Sample'])['Ribosomal+Mitochondrial_Fraction'].median()
df['%Mitochondrial']=df['Mitochondrial_Reads']*100/df['Barcoded_Reads']
# df['%Ribosomal']=df['Ribosomal_Reads']*100/df['Barcoded_Reads']
df['%Intergenic_Reads']=df['Intergenic_Reads']*100/df['Barcoded_Reads']

beta2_spec2 = df[((df["Well_Type"] == "Sample") | 
                  (df["Well_Type"] == "sample") | 
                  (df["Well_Type"] == "PBMC")) & 
                 (df['Group'] == shasta+'_Raw')]
cx_spec2=df[((df["Well_Type"] == "PBMC")| (df["Well_Type"] == "Sample") | (df["Well_Type"] == "sample")) & (df['Group']==cx+'_Raw')] #CX1018

# sspro_spec2_df = func.sspro_spec2(beta2_spec2,'Barcoded_Reads')
# # sspro_spec2_df = func.sspro_spec2(cx_spec2,'Barcoded_Reads')
# func.plot_table_from_df(sspro_spec2_df[0],'no index',base_width=2)


#USING XUAN df for sspro spec 3 and side by side comparison spec 1, dropped all non-samples
spec3_list_raw =[df[df['Sample']==specific_order[0]],df[df['Sample']==specific_order[1]]] #beta1 raw and cx raw
sspro_spec3_df_raw = func.sspro_spec3(spec3_list_raw,'%Intergenic_Reads')

#sspro spec 3 for 100K reads
spec3_list_100K =[df[df['Sample']==specific_order[2]],df[df['Sample']==specific_order[3]]] #beta1 100k and cx 100k
sspro_spec3_df_100K = func.sspro_spec3(spec3_list_100K,'%Intergenic_Reads')

sspro_spec3_df_raw[1].set_index('Sample',inplace=True)
sspro_spec3_df_100K[1].set_index('Sample',inplace=True)


sspro_spec3_intergenic = pd.concat([sspro_spec3_df_raw[1], sspro_spec3_df_100K[1]],axis=0)
spec3_ratio_raw= sspro_spec3_df_raw[1].loc[specific_order[0]]/sspro_spec3_df_raw[1].loc[specific_order[1]]
spec3_ratio_100K= sspro_spec3_df_100K[1].loc[specific_order[2]]/sspro_spec3_df_100K[1].loc[specific_order[3]]
sspro_spec3_intergenic['Spec'] = sspro_spec3_intergenic['Tenth_Percentile'].apply(lambda x: 'Fail' if x < 70 else 'Pass')
                                                                                
# NOTE: MUST BE IN THIS ORDER: df,CX,shasta,CX_100K,shasta_100K
sscomp_spec1_df = func.sscomparison_spec1(df,specific_order[1],specific_order[0],specific_order[3],specific_order[2])

sscomp_spec1_df['Spec'] = None

# sscomp_spec1_df['Spec'] = sscomp_spec1_df['Shasta to CX1018 ratio'].apply(lambda x: "Fail" if x < 0.9 else "Pass")
sscomp_spec1_df.iloc[1, sscomp_spec1_df.columns.get_loc('Spec')] = "Fail" if sscomp_spec1_df.iloc[1, sscomp_spec1_df.columns.get_loc(f'{shasta}_Raw to {cx}_Raw ratio')] < gene_sensitivity else "Pass" #adjusted to 80% for hammerhead validation


#IMAGE GENERATION
target_cols=['No_of_Genes']
target_cols_fraction = ['Ribosomal_Reads','Mitochondrial_Reads','Intergenic_Reads','Total_Exon_Reads','Total_Intron_Reads']

images_violin = reportfunc.plot(df, target_cols, 'Sample',specific_order,save_path = image_folder)
images_BarcodedReads = reportfunc.plotBarcodedReads(df,specific_order,save_path = image_folder)

images_fraction = reportfunc.plot_by_fraction_to_html(df, target_cols_fraction, 'Barcoded_Reads', 'Sample', specific_order, save_path = image_folder)  

heat_barcoded = reportfunc.gen_heatmap(list_of_df, 'Barcoded_Reads', 'Sample', 'log',image_folder +"Barcoded.png")
heat_inter = reportfunc.gen_heatmap(list_of_df, '%Intergenic_Reads', 'Sample', 'bwr',image_folder +"Intergenic.png")
heat_mito = reportfunc.gen_heatmap(list_of_df, '%Mitochondrial', 'Sample', 'bwr',image_folder +"Mitochondrial.png")
heat_ribo = reportfunc.gen_heatmap(list_of_df, '%Ribosomal', 'Sample', 'bwr',image_folder +"Ribosomal.png")
heat_gene = reportfunc.gen_heatmap(list_of_df, 'No_of_Genes', 'Sample', 'log',image_folder +"Genes.png")




#%%


with open(savepath,'w') as _file:
    _file.write("<h1>"+title+"</h1>")
    _file.write(notes)
    #spec 1
    _file.write("<h2>SSPro Spec 1: Background noise</h2>")
    _file.write("<h3>"+shasta+"</h3>")
    _file.write("The average number of reads from Neg Ctrl wells has to be <5% of that of Pos Ctrl wells (2pg RNA) &#8594; " + str(sspro_spec1_df_shasta[1]))
    _file.write("<br> "+ "The AVE.+3SD of Neg Ctrl wells should not exceed any of Pos Ctrl wells &#8594; " + str(sspro_spec1_df_shasta[2]))
    _file.write(sspro_spec1_df_shasta[0].to_html() + "\n\n")
    
    _file.write("<h3>" + cx + "</h3>")    
    _file.write("The average number of reads from Neg Ctrl wells has to be <5% of that of Pos Ctrl wells (2pg RNA) &#8594; " + str(sspro_spec1_df_CX[1]))
    _file.write("<br> "+ "The AVE.+3SD of Neg Ctrl wells should not exceed any of Pos Ctrl wells &#8594; " + str(sspro_spec1_df_CX[2]))
    _file.write(sspro_spec1_df_CX[0].to_html() + "\n\n<br>")
    #spec 2
    _file.write("<h2>SSPro Spec 2: Variation</h2>")
    _file.write("The number of reads between samples (PBMCs): 85% of cells should be in the range of an order of magnitude ( &#8804;10-fold difference)")
    _file.write(sspro_spec2_df.to_html())
    _file.write(shasta + " (97.5% ~ 12.5%): " +'%.2f'%sspro_spec2_df_shasta[3]+"-fold difference  &#8594; " +spec2_shasta + "<br>")
    _file.write(cx+ " (97.5% ~ 12.5%): " +'%.2f'%sspro_spec2_df_CX[3]+"-fold difference  &#8594; " +spec2_cx + "<br>")
    #spec 3
    _file.write("<h2>SSPro Spec 3: gDNA contamination</h2>")
    _file.write("% genomic DNA for PBMCs: 70% of cells should have <10% intergenic reads.")
    _file.write("<h3>Raw</h3>")
    _file.write(sspro_spec3_df_raw[0].to_html())
    _file.write("<h3>Downsampling to 100K</h3>")
    _file.write(sspro_spec3_df_100K[0].to_html())
    _file.write("%cells have <10% intergenic reads:")
    _file.write(sspro_spec3_intergenic.to_html())
    #spec 4
    _file.write("<h2>SSPro Spec 4: Unmapped Reads</h2>")
    _file.write("% unmapped for PBMCs should be <10%. Controls are excluded.")
    _file.write(sspro_spec4.to_html())
    _file.write("<h2>SSPro Spec 5: Seurat Cell Type Annotation</h2>")
    _file.write("Clear separation of CD4+ and CD8+ T cells clusters (UMAP) for PBMCs")
    _file.write("<h2>SSPro Spec 6: Seurat Cell Type Annotation</h2>")
    _file.write("Identification of clusters (UMAP) of CD4+ T cells, CD8+ T cells, NK cells, Monocytes, B cells for PBMCs samples.")
    #insert seurat data
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write(f"<h3>Seurat:{shasta}_Raw</h3>")
    _file.write(f'<img src="{image_folder}Seurat_Raw_Shasta.png" style="width:500px;"/>')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write(f"<h3>Seurat: {cx}</h3>")
    _file.write(f'<img src="{image_folder}Seurat_Raw_CX.png" style="width:500px;" />')
    _file.write("</div>")
    _file.write("</div>")

    _file.write("<h2>SSPro Side-by-side Spec 1: </h2>")
    _file.write(f"Sensitivity: > {gene_sensitivity*100}% compared to SSPro {cx}. Controls are excluded.")
    _file.write(sscomp_spec1_df.to_html())
    #no of genes
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write(f"<h3>{target_cols[0]}</h3>")
    _file.write(f'<img src="{images_violin[0]}"  />')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write("<h3>No_of_Genes</h3>")
    _file.write(f'<img src="{heat_gene}"/>')
    _file.write("</div>")
    _file.write("</div>")

    _file.write("<h2>SSPro Side-by-side Spec 2: </h2>")
    _file.write("% Barcoded reads: > 90% compared to SSPro ICELL8cx.")
    _file.write(sscomp_spec2.to_html(index=False))
    _file.write("<h2>SSPro Side-by-side Spec 3: </h2>")
    _file.write("% Cells passing Cogent DS filter: > 90% compared to SSPro ICELL8cx.")
    _file.write(sscomp_spec3.to_html())
    # _file.write(sscomp_spec3_02.to_html())
    _file.write("<h2>SSPro Side-by-side Spec 4: </h2>")
    _file.write("% Exon reads: > 90% compared to SSPro ICELL8cx.")
    _file.write(sscomp_spec4.to_html())
    _file.write("<h2>SSPro Side-by-side Spec 5: </h2>")
    _file.write("%rRNA + Mitochondrial RNA in PBMCs: < 10% higher compared to SSPro ICELL8cx. (Excluding Pos Ctrls)")
    _file.write(sscomp_spec5.to_html())
    _file.write("<h1>Details</h1>")


# barcoded reads
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write("<h3>Barcoded Reads</h3>")
    _file.write(f'<img src="{images_BarcodedReads}"/>')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write("<h3>Barcoded Reads</h3>")
    _file.write(f'<img src="{heat_barcoded}"/>')
    _file.write("</div>")
    _file.write("</div>")

# intergenic
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write(f"<h3>{target_cols_fraction[2]}</h3>")
    _file.write(f'<img src="{images_fraction[2]}" />')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write(f"<h3>{target_cols_fraction[2]}</h3>")
    _file.write(f'<img src="{heat_inter}"/>')
    _file.write("</div>")
    _file.write("</div>")

# mitochondrial
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write(f"<h3>{target_cols_fraction[1]}</h3>")
    _file.write(f'<img src="{images_fraction[1]}" />')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write(f"<h3>{target_cols_fraction[1]}</h3>")
    _file.write(f'<img src="{heat_mito}"/>')
    _file.write("</div>")
    _file.write("</div>")



# ribosomal
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write(f"<h3>{target_cols_fraction[0]}</h3>")
    _file.write(f'<img src="{images_fraction[0]}" />')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write(f"<h3>{target_cols_fraction[0]}</h3>")
    _file.write(f'<img src="{heat_ribo}"/>')
    _file.write("</div>")
    _file.write("</div>")


# intron and exon
   
    _file.write('<div style="display: flex;">')
    # Iterate over the images and their corresponding titles
    for title_intron, image_path in zip(target_cols_fraction[3:], images_fraction[3:]):
        _file.write("<div>")
        _file.write(f"<h3>{title_intron}</h3>")  # Write the title
        _file.write(f'<img src="{image_path}" />')  # Write the image
        _file.write("</div>")
    
    _file.write("</div>")

#insert cogentDS data
    _file.write('<div style="display: flex;">')

    _file.write("<div>")
    _file.write(f"<h3>CogentDS:{shasta}_Raw</h3>")
    _file.write(f'<img src="{image_folder}Cogent_Raw_Shasta.png" style="width:500px;"/>')
    _file.write("</div>")

    
    _file.write("<div>")
    _file.write(f"<h3>CogentDS: {cx}_Raw</h3>")
    _file.write(f'<img src="{image_folder}Cogent_Raw_CX.png" style="width:500px;"/>')
    _file.write("</div>")
    _file.write("</div>")
   





# #%%    
# # posCtrl = df[df['Well_Type']=='Pos_Ctrl']
# # list_of_df = [group for _, group in df.groupby('Sample')]
# # test=posCtrl[posCtrl['Sample'].isin(['CX1018_Raw','Beta1_Raw'])]
# # test=df[df['Sample'].isin(specific_order[:2])]

# test=df[df['Sample'].isin([specific_order[0],specific_order[1]])]

# # func.plot(test, ['Barcoded_Reads'], 'Sample',specific_order[:2], base_width=5)
# func.plot(test, ['Barcoded_Reads'],'Sample',[specific_order[0],specific_order[1]], base_width=8) #only beta1

# test=df[df['Sample'].isin(specific_order[2::])]

# func.plot(df, ['No_of_Genes'], 'Sample',specific_order, base_width=5)
# func.plot_by_fraction(df, ['Mapped_Reads'], 'Barcoded_Reads', 'Sample', specific_order, base_width=4)
# func.plot_by_fraction(df, ['Uniquely_Mapped_Reads'], 'Mapped_Reads', 'Sample',specific_order, base_width=10)
# func.plot_by_fraction(df, ['Total_Exon_Reads', 'Total_Intron_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=7)
# func.plot_by_fraction(df, ['Ribosomal+Mitochondrial_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=4)
# func.plot_by_fraction(df, ['Total_Exon_Reads'], 'Barcoded_Reads', 'Sample', specific_order,base_width=10)
# func.plot_by_fraction(df, ['Mitochondrial_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=5)
# func.plot_by_fraction(df, ['Ribosomal_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=5)
# func.plot_by_fraction(df, ['Intergenic_Reads'], 'Barcoded_Reads', 'Sample',specific_order, base_width=5)

#%%

# func.gen_heatmap_individual(list_of_df,'Barcoded_Reads','Sample','log', base_width=10)

# func.gen_heatmap(list_of_df,'Barcoded_Reads','Sample','log', base_width=10)
# func.gen_heatmap_individual(list_of_df,'Barcoded_Reads','Sample','autumn', base_width=10)

# func.gen_heatmap_individual(list_of_df,'No_of_Genes','Sample','log', base_width=7)

# func.gen_heatmap(list_of_df,'%Intergenic_Reads','Sample','bwr',base_width=7)
# func.gen_heatmap(list_of_df,'%Ribosomal','Sample','bwr',base_width=7)
# func.gen_heatmap(list_of_df,'%Mitochondrial','Sample','bwr',base_width=7)
# func.plot_barchart(list_of_df,'Barcoded_Reads','Row','Sample',base_width=5)

# func.gen_heatmap(list_of_df,'Total_Intron_Reads','Sample','log', base_width=7)



#%%
# experiment_start_date='2023-12-18'
# master_path = "S:/XuanLi/Shasta_SSPro/141037C_20231218111254_Beta1_SSProTest4/Logs/"
# log_path = master_path + "log-202312.txt"
# dewpointlog_path= master_path + "DewPointLog-20231218.tsv"
# dewpointlog_path2= master_path + "DewPointLog-20231219.tsv"

# instrument="1218_Beta1"
# # test=func.log_parse(dewpointlog_path,instrument)
# beta1=func.read_log(log_path,dewpointlog_path,dewpointlog_path2,instrument,experiment_start_date)
# func.plot_table_from_df(beta1, instrument, base_width=15)
# dfi.export(beta1,master_path+instrument+'.png')


# #%% CX RH logs 
# import SSPro_func as func
# master_path = "S:/XuanLi/Shasta_SSPro/2024.02.14.09.54-136199_CX1018_SSPro_Test2/"
# log_path = master_path + "Debug_0214.log" #this is for instrument logs 
# start_date = "2024.02.14"  # Replace with your desired start date in this format 2023.10.09
# input_path = master_path+"20240214_CX1018_teraterm.log" #this is for teraterm logs
# output_folder=master_path+"/test"
# instrument = "0214_CX1018"
# df=func.read_log_CX(log_path, input_path, output_folder, instrument, start_date)
# func.plot_table_from_df(df, instrument, base_width=6)
# dfi.export(df,master_path+"RH log_"+instrument+'.png')




