# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:28:15 2024

@author: leongs
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import seaborn as sns
import base64
from io import BytesIO
import os
import numpy as np


def plotBarcodedReads(df,specific_order,save_path):
    specific_order = [specific_order[0],specific_order[1]]
    test=df[df['Sample'].isin(specific_order)]
    path = plot(test, ['Barcoded_Reads'],'Sample',[specific_order[0],specific_order[1]], base_width=10,base_height=5,save_path=save_path)
    return path[0]




def plot(df, target_cols, group_column, specific_order, p_type='violin', base_width=10, base_height=10, rotate=0, save_path=None):
    # Create a directory to save the images if save_path is provided
    if save_path:
        os.makedirs(save_path, exist_ok=True)

    image_paths = []  # Store the paths of saved images

    for column in target_cols:
        # Create a new figure for each plot
        fig, ax = plt.subplots(figsize=(base_width, base_height))

        if p_type == 'violin':
            sns.violinplot(data=df, x=group_column, y=column, width=0.4, hue='Instrument', ax=ax)
        # elif p_type == 'box':
        #     sns.boxplot(data=df, x=group_column, y=column, ax=ax)

        df["Sample"] = pd.Categorical(df["Sample"], categories=specific_order, ordered=True)
        medians = df.groupby("Sample")[column].median().reindex(specific_order)

        # Plot median values (optional)
        for xtick, median in enumerate(medians):
            # ax.text(xtick, median, f"{median:.2f}", horizontalalignment='center', color='black', fontsize=10)
            ax.text(xtick+0.3, median, f"{median:.2f}", horizontalalignment='right', color='black', fontsize=10)


        ax.set_title(f'{column}')
        
        if rotate == 0:
            ax.set_xticklabels(ax.get_xticklabels())
        else:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotate, ha='right')

        fig.tight_layout()
        
        # Save the figure as a PNG file if save_path is provided
        if save_path:
            filename = os.path.join(save_path, f'{column}.png')
            fig.savefig(filename, format='png', bbox_inches='tight', pad_inches=0)
            image_paths.append(filename)  # Append the path to the list
        
        plt.close(fig)  # Close the figure to free memory

    return image_paths  # Return the list of image paths


#%%
def plot_by_fraction_to_html(input_df, target_cols, base_col, group_column, specific_order, p_type='violin', base_width=6, rotate=0, title=" Title", save_path=None):
    image_paths = []  # Store the paths of saved images

    # Create a directory to save the images if save_path is provided
    if save_path:
        os.makedirs(save_path, exist_ok=True)

    for target_col in target_cols:
        # Create a new figure for each target column
        fig, ax = plt.subplots()
        
        # Create dataframe by fraction
        df = pd.DataFrame(data=[], columns=['Fraction', 'Exp_Key'])
        frac_list = (input_df[target_col] / input_df[base_col]).values.tolist()
        data = {
            'Fraction': frac_list,
            group_column: input_df[group_column].values.tolist(), 
        }
        df = pd.concat([df, pd.DataFrame(data)], ignore_index=True)
        df[["Instrument", 'Read Depth']] = df["Sample"].str.split("_", expand=True)

        ax = sns.violinplot(data=df, x=group_column, y='Fraction', width=0.4, hue='Instrument', ax=ax)

        df["Sample"] = pd.Categorical(df["Sample"], categories=specific_order, ordered=True)
        medians = df.groupby('Sample')['Fraction'].median()
        for xtick, median in enumerate(medians):
            ax.text(xtick+0.5, median, f"{median:.2f}", horizontalalignment='right', color='black', fontsize=10)

        ax.set_title(target_col + ' - Fraction')
        if rotate != 0:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=rotate, ha='right')

        fig.set_figwidth(base_width)
        if save_path:
            filename = os.path.join(save_path, f'{target_col}.png')
            fig.savefig(filename, format='png')
            image_paths.append(filename)  # Append the path to the list

        plt.close(fig)  # Close the figure to free memory
    
    return image_paths  # Return the list of image paths


#%%


def gen_heatmap(df_list, param, group_column, colour, filename, base_width=8):
    heatmap_images = []

    # Step 1: Calculate global min and max
    global_min = float('inf')
    global_max = float(0)
    # global_max = float('-inf')

    for df in df_list:
        unique_groups = df[group_column].unique()
        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)

            # Ensure the data is numeric and handle non-numeric or NaN values
            if heatmap_data.empty:
                continue
            
            heatmap_data = heatmap_data.to_numpy()  # Convert to numpy array for min/max calculation
            # heatmap_data[np.isnan(heatmap_data)] = -np.inf
            
            data_min = np.nanmin(heatmap_data)
            data_max = np.nanmax(heatmap_data)

            if data_min < global_min:
                global_min = data_min
            if data_max > global_max:
                global_max = data_max

    # Step 2: Calculate the number of rows and columns for subplots
    num_rows = (len(df_list) + 1) // 2
    num_cols = min(len(df_list), 2)

    # Create a single figure with subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(base_width * num_cols, base_width * num_rows))
    fig.subplots_adjust(hspace=0.3)  # Adjust the vertical spacing between subplots

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        # Determine the subplot position based on the index
        if num_rows == 1:
            ax = axes[index % num_cols]
        else:
            ax = axes[index // num_cols, index % num_cols]

        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)

            # Ensure the data is numeric and handle non-numeric or NaN values
            if heatmap_data.empty:
                continue
            
            heatmap_data = heatmap_data.to_numpy()  # Convert to numpy array for plotting
            # heatmap_data[np.isnan(heatmap_data)] = -np.inf
            
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
        if num_rows == 1:
            fig.delaxes(axes[i % num_cols])
        else:
            fig.delaxes(axes[i // num_cols, i % num_cols])

    plt.suptitle(param, fontsize=20, y=0.9)

    # Save the plot to a buffer and encode it as base64
    img_buf = BytesIO()
    fig.savefig(img_buf, format='png', bbox_inches='tight', pad_inches=0)
    img_buf.seek(0)
    img_base64 = base64.b64encode(img_buf.getvalue()).decode('utf-8')
    heatmap_images.append(img_base64)
    
    # Save the plot to a file
    plt.savefig(filename, bbox_inches='tight', pad_inches=0)
    plt.close(fig)

    return filename


#%%
def OLD_gen_heatmap(df_list, param, group_column, colour, filename, base_width=8):
    heatmap_images = []

    # Calculate the number of rows and columns for subplots
    num_rows = (len(df_list) + 1) // 2
    num_cols = min(len(df_list), 2)

    # Create a single figure with subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(base_width * num_cols, base_width * num_rows))

    fig.subplots_adjust(hspace=0.3)  # Adjust the vertical spacing between subplots

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        # Determine the subplot position based on the index
        if num_rows == 1:
            ax = axes[index % num_cols]
        else:
            ax = axes[index // num_cols, index % num_cols]

        for group_value in unique_groups:
            heatmap_data = df[df[group_column] == group_value].pivot(index='Row', columns='Col', values=param)

            if colour == 'log':
                im = ax.imshow(heatmap_data, cmap='bwr', extent=[0, 72, 72, 0], origin='upper', norm=LogNorm())
            else:
                im = ax.imshow(heatmap_data, cmap=colour, extent=[0, 72, 72, 0], origin='upper')

            ax.set_title(f'{group_column}: {group_value}')
            fig.colorbar(im, ax=ax, label=f'{"Log Value of " if colour == "log" else ""}{param}')

            ax.set_xlabel('Col', fontsize=15)
            ax.set_ylabel('Row', fontsize=15)

    # Close any unused subplots
    for i in range(len(df_list), num_rows * num_cols):
        if num_rows == 1:
            fig.delaxes(axes[i % num_cols])
        else:
            fig.delaxes(axes[i // num_cols, i % num_cols])

    plt.suptitle(param, fontsize=20, y=0.9)
    
    img_buf = BytesIO()
    fig.savefig(img_buf, format='png')
    img_buf.seek(0) 

    img_base64 = base64.b64encode(img_buf.getvalue()).decode('utf-8')
    heatmap_images.append(img_base64)
    plt.savefig(filename, bbox_inches='tight', pad_inches=0)  # Save the image without white space
    plt.close(fig)

    return filename



def gen_heatmap_individual_test_to_base64(df_list, param, group_column, colour):
    # Initialize variables to store global minimum and maximum values
    global_min = float('inf')
    global_max = float('-inf')
    heatmap_images = []

    for df in df_list:
        # Find the minimum and maximum values within the current DataFrame
        min_val = df[param].min()
        max_val = df[param].max()
        
        # Update global minimum and maximum values if necessary
        global_min = min(global_min, min_val)
        global_max = max(global_max, max_val)

    new_df = pd.DataFrame([(i, j) for i in range(72) for j in range(72)], columns=['Row', 'Col'])

    for index, df in enumerate(df_list):
        unique_groups = df[group_column].unique()

        for group_value in unique_groups:
            fig, ax = plt.subplots()
            merged_df = pd.merge(new_df, df[df[group_column] == group_value], on=['Row', 'Col'], how='left')

            heatmap_data = merged_df.pivot(index='Row', columns='Col', values=param)

            if colour == 'log':
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

            y_ticks = list(range(0, 75, 10))
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_ticks)

            x_ticks = list(range(0, 75, 10))
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticks)

            img_buf = BytesIO()
            fig.savefig(img_buf, format='png')
            img_buf.seek(0) 

            img_base64 = base64.b64encode(img_buf.getvalue()).decode('utf-8')
            heatmap_images.append(img_base64)

            plt.close(fig)

    return heatmap_images
