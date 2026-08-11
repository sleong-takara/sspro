# SSPro Script README

This folder contains SSPro analysis scripts for comparing Shasta SSPro / hammerhead syringe runs against ICELL8cx or baseline Shasta runs, generating SSPro validation metrics, plotting QC summaries, generating heatmaps, parsing instrument/RH logs, and optionally creating an HTML-style report.

> **Important:** These scripts are currently written as notebook-style / Spyder-style linear analysis scripts with hard-coded input paths near the top of each file. They are not command-line tools yet. For each new run, edit the configuration block at the top of the script, then run the script from the folder that contains the helper modules.

---

## Folder contents

| Script | Role |
|---|---|
| `SSPro_pipeline.py` | Main exploratory SSPro analysis pipeline for comparing Shasta/CX-style runs, generating SSPro specs, plots, heatmaps, and RH/log summaries. |
| `SSPro_pipeline_report.py` | Report-oriented pipeline that calculates SSPro specs and comparison specs, saves plots to an image folder, and writes an HTML report file. |
| `SSPro_pipeline_nocx.py` | Variant of the main pipeline for Shasta-only or no-CX comparisons. Useful when there is no CX1018 comparator dataset. |
| `SSPro_func.py` | Shared utility functions for SSPro specs, data cleanup, plots, heatmaps, barcharts, and instrument/RH log parsing. |
| `SSPro_report_func.py` | Shared report plotting functions used by `SSPro_pipeline_report.py` to save images and generate report-ready plot outputs. |

Compiled files such as `*.pyc` may exist in `__pycache__`. They are generated automatically by Python and do not need to be run directly.

---

## Recommended environment

### Python packages

The scripts import the following packages:

```bash
pip install pandas numpy matplotlib seaborn plotly mpld3 tabulate
```

Some older or optional code references `dataframe_image`. Install it only if you plan to enable those commented-out export lines:

```bash
pip install dataframe_image
```

### Local helper modules

Keep these helper scripts in the same folder as the pipeline scripts, or make sure the folder is on `PYTHONPATH`:

```text
SSPro_func.py
SSPro_report_func.py
parse_CX_log.py   # imported by SSPro_func.py
```

If `parse_CX_log.py` is not available and you are not parsing CX logs, you can still review the analysis sections, but importing `SSPro_func.py` may fail unless that dependency is present.

---

## Expected input files

Most analysis paths point to unzipped Mappa/Cogent-style output folders. The scripts expect files such as:

```text
analysis_stats.csv
WellList*.txt or files containing "WellList" in the name
demux_counts_all.csv       # required by SSPro_pipeline_report.py
instrument log files       # only for log parsing sections
DewPointLog-*.tsv          # only for RH/dewpoint parsing sections
```

The exact folder layout is controlled by variables such as `folder`, `folder1`, `path0`, `path1`, `path2`, and `path3` near the top of each pipeline script.

---

## General workflow

1. Activate the Python environment.
2. Open the script you want to use.
3. Edit the configuration variables at the top of the script.
4. Confirm the expected CSV/log files exist in each input folder.
5. Run the script from the script folder so local imports resolve.
6. Review generated plots, printed tables, heatmaps, and report outputs.

Example:

```bash
cd "C:/Users/leongs/OneDrive - Takara Bio USA, Inc/2. Scripts/SSPro_Script"
python SSPro_pipeline.py
```

On Linux / WSL, use the mounted equivalent path, for example:

```bash
cd "/mnt/c/Users/leongs/OneDrive - Takara Bio USA, Inc/2. Scripts/SSPro_Script"
python SSPro_pipeline.py
```

---

# Script details

## 1. `SSPro_pipeline.py`

### Purpose

`SSPro_pipeline.py` is the main interactive SSPro analysis workflow. It compares multiple SSPro output folders, calculates SSPro spec tables, merges well-list metadata with analysis statistics, produces read/gene/fraction plots, generates chip heatmaps, and includes sections for parsing Shasta and CX RH/log files.

### Main outputs

Depending on which sections are run, the script can produce:

- SSPro spec 1 background-noise summary table.
- SSPro spec 2 read-variation quantile table.
- SSPro spec 3 intergenic-read summary.
- SSPro spec 4 unmapped-read summary.
- Side-by-side comparison metrics such as median detected genes.
- Violin plots / barcoded-read plots.
- Fraction plots for mapped, uniquely mapped, exon, intron, ribosomal, mitochondrial, and intergenic reads.
- 72 x 72 chip heatmaps for barcoded reads, detected genes, intergenic reads, ribosomal reads, mitochondrial reads, and intron reads.
- RH / temperature log summary tables when log parsing sections are enabled.

### Key configuration variables to edit

At the top of the file, edit:

```python
folder = r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"
folder1 = r"C:/Users/leongs/OneDrive - Takara Bio USA, Inc/0. Projects/Shasta instrument/Hammerhead syringe/20260323_SSPro_S2016/"

path0 = folder1 + r"20260323_Hammerhead SSPro rerun_S2016"          # Shasta / hammerhead raw
path1 = folder + "241209_SSPro_Beta2"                               # comparator raw
path2 = folder1 + r"20260323_Hammerhead SSPro rerun_S2016_100K"     # Shasta / hammerhead downsampled
path3 = folder + "241209_SSPro_Beta2_100K"                           # comparator downsampled

shasta = "v2"
cx = "v1"
specific_order = [shasta + "_Raw", cx + "_Raw", shasta + "_100K", cx + "_100K"]
sample_name = ["PBMC", "Sample"]
```

### Usage example: compare Shasta hammerhead vs beta baseline

```python
folder = r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"
folder1 = r"C:/Users/leongs/OneDrive - Takara Bio USA, Inc/0. Projects/Shasta instrument/Hammerhead syringe/20260323_SSPro_S2016/"

path0 = folder1 + r"20260323_Hammerhead SSPro rerun_S2016"
path1 = folder + "241209_SSPro_Beta2"
path2 = folder1 + r"20260323_Hammerhead SSPro rerun_S2016_100K"
path3 = folder + "241209_SSPro_Beta2_100K"

shasta = "v2"
cx = "v1"
specific_order = ["v2_Raw", "v1_Raw", "v2_100K", "v1_100K"]
sample_name = ["PBMC", "Sample"]
```

Run:

```bash
python SSPro_pipeline.py
```

### Notes and cautions

- This script is section-based. If using Spyder / VS Code cells, run sections in order unless you know which variables have already been created.
- The code assumes the `analysis_stats.csv` columns include fields such as `Sample`, `Barcode`, `Barcoded_Reads`, `Mapped_Reads`, `Unmapped_Reads`, `No_of_Genes`, `Ribosomal_Reads`, `Mitochondrial_Reads`, `Intergenic_Reads`, `Exon_Reads`, `Ambiguous_Exon_Reads`, `Intron_Reads`, and `Ambiguous_Intron_Reads`.
- The heatmap functions assume chip coordinates are available as `Row` and `Col` in the well list / merged dataframe.
- Several cleanup / row-removal sections are commented as "DO NOT RUN" and should remain untouched unless you intentionally need to exclude failed wells.
- The log parsing sections require separate instrument log files and are independent from the main analysis_stats workflow.

---

## 2. `SSPro_pipeline_report.py`

### Purpose

`SSPro_pipeline_report.py` is the report-generation version of the SSPro analysis workflow. It calculates validation metrics, builds plot images using `SSPro_report_func.py`, and writes an HTML report to a run-specific folder.

### Main outputs

- HTML report file, for example:

```text
<masterpath>/<title>/<title>.html
```

- Plot image folder, for example:

```text
<masterpath>/<title>/images/
```

- Report-ready images for:
  - `No_of_Genes`
  - `Barcoded_Reads`
  - read fractions such as ribosomal, mitochondrial, intergenic, exon, and intron reads
  - heatmaps for barcoded reads, intergenic reads, mitochondrial reads, ribosomal reads, and detected genes

### Key configuration variables to edit

```python
folder = r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"
folder1 = r"C:/Users/leongs/OneDrive - Takara Bio USA, Inc/0. Projects/Shasta instrument/Hammerhead syringe/SSPro double hammerhead/"

path0 = folder1 + r"SSPro 20nL stacked validation hammerhead_noDownsampling"
path1 = folder + "241209_SSPro_Beta2"
path2 = folder1 + r"SSPro 20nL stacked validation hammerhead_downsampling50K"
path3 = folder + "241209_SSPro_Beta2_50K"

title = "SSPro hammerhead validation 2025-12-09"
masterpath = folder1 + "SSPro_Reports/"
shasta = "v2"
cx = "v1"
notes = "Downsampled to 50k"
gene_sensitivity = 0.8
```

Also update pass/fail metadata values manually if needed:

```python
cx_total_cells = 1
cx_pass = 1
shasta_total_cells = 1
shasta_pass = 1
```

### Usage example: generate an HTML validation report

```python
title = "SSPro hammerhead validation 2025-12-09"
masterpath = folder1 + "SSPro_Reports/"
shasta = "v2"
cx = "v1"
gene_sensitivity = 0.8
notes = "Downsampled to 50k"
```

Run:

```bash
python SSPro_pipeline_report.py
```

Expected output:

```text
SSPro_Reports/
└── SSPro hammerhead validation 2025-12-09/
    ├── SSPro hammerhead validation 2025-12-09.html
    └── images/
        ├── Barcoded.png
        ├── Intergenic.png
        ├── Mitochondrial.png
        ├── Ribosomal.png
        └── Genes.png
```

### Required input files

In addition to `analysis_stats.csv` and `WellList` files, this script reads:

```text
demux_counts_all.csv
```

for each run folder.

### Notes and cautions

- This script writes output folders using `os.makedirs(image_folder, exist_ok=True)`.
- `specific_order` is constructed as `[shasta+'_Raw', cx+'_Raw', shasta+'_50K', cx+'_50K']`.
- The script contains manually populated values for cell counts and pass counts. Update them before using the report for a final validation summary.
- Some ratio column names in the current code are dynamically formatted and may need review if `shasta` or `cx` labels change.

---

## 3. `SSPro_pipeline_nocx.py`

### Purpose

`SSPro_pipeline_nocx.py` is a Shasta-only / no-CX version of the SSPro analysis workflow. It is useful for cases where you want to evaluate one Shasta run against its downsampled version or do not have a CX1018 comparator dataset.

### Main outputs

Depending on sections run, the script can produce:

- SSPro spec 1 background-noise table.
- SSPro spec 2 read-variation table for Shasta data.
- SSPro spec 3 intergenic summary.
- SSPro spec 4 unmapped-read summary.
- Shasta raw vs downsampled detected-gene comparison using `sscomparison_spec1_nocx`.
- Violin plots and fraction plots.
- Barcoded-read and detected-gene heatmaps.
- RH/log summaries if those sections are enabled.

### Key configuration variables to edit

```python
folder = r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"
path0 = folder + "241209_SSPro_Beta2"                # Shasta raw
path2 = folder + "250325_SSPro_25nLDoubleDisp_50K"   # Shasta downsampled

title = "250325 25nL dispense"
shasta = "Beta-2"
cx = "cx-1018"       # retained as a label, but CX paths are commented out
specific_order = [shasta + "_Raw", shasta + "_50K"]
sample_name = ["Sample", "sample"]
```

### Usage example: Shasta raw vs 50K downsampled comparison

```python
folder = r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/"
path0 = folder + "241209_SSPro_Beta2"
path2 = folder + "250325_SSPro_25nLDoubleDisp_50K"

title = "250325 25nL dispense"
shasta = "Beta-2"
specific_order = ["Beta-2_Raw", "Beta-2_50K"]
sample_name = ["Sample", "sample"]
```

Run:

```bash
python SSPro_pipeline_nocx.py
```

### Notes and cautions

- Keep CX-related code commented unless you restore valid `path1` and `path3` inputs.
- Some plotting lines still reference `specific_order[2]` and `specific_order[3]`, which are not available in a two-condition no-CX setup. Review those sections before running the full script end to end.
- Best use is section-by-section execution, especially after `df`, `list_of_df`, and `specific_order` are created.

---

## 4. `SSPro_func.py`

### Purpose

`SSPro_func.py` contains shared helper functions used by the SSPro pipelines. It centralizes SSPro spec calculations, comparison summaries, plotting, heatmap generation, well-list cleanup, and RH/log parsing.

### Important functions

#### Spec and comparison functions

| Function | What it does |
|---|---|
| `sspro_spec1_defunct(df1)` | Calculates background-noise metrics using `Sample` labels, including negative-control mean, SD, mean + 3SD, and positive-control comparisons. |
| `sspro_spec1(df1)` | Similar background-noise calculation using `Well_Type` labels. |
| `sspro_spec2(df1, column)` | Calculates quantiles for a read-count-like column to evaluate read spread / fold range. |
| `sspro_spec3(df1_list, column)` | Summarizes intergenic-read-related quantiles for one or more dataframes. |
| `sspro_spec4(df1)` | Calculates percent unmapped reads from summed read columns. |
| `sscomparison_spec1(df, shasta, CX, shasta_100K, CX_100K)` | Compares median detected genes between Shasta/CX raw and downsampled conditions. |
| `sscomparison_spec1_nocx(df, shasta, shasta_100K)` | Shasta-only detected-gene comparison. |
| `sscomparison_spec1_doublet(...)` | Compares gene sensitivity in doublet-style human/mouse or instrument comparisons. |

#### File and dataframe helpers

| Function | What it does |
|---|---|
| `search_files_in_directory(directory_path)` | Finds the first file containing `WellList` in a directory and reads it as a tab-delimited dataframe. |
| `cleanup_barcode(df)` | Removes `+` characters from the `Barcode` column. |
| `flatten_dict(d, parent_key='', sep='_')` | Flattens nested dictionaries. |

#### Plotting and heatmaps

| Function | What it does |
|---|---|
| `plot_table_from_df(df, index, base_width, save_path=None)` | Renders a dataframe as a matplotlib table. |
| `plot_table_from_df_ordered(...)` | Renders a sorted / ordered dataframe table. |
| `plot(...)` | Generates violin-style plots for selected columns. |
| `plot_by_fraction(...)` | Plots fractions of target read classes relative to a base column. |
| `plot_barchart(...)` | Generates row/column barcharts by group. |
| `gen_heatmap(...)` | Generates multi-condition chip heatmaps. |
| `gen_heatmap_individual(...)` | Generates individual chip heatmaps. |
| `gen_heatmap_filtered(...)` | Generates filtered heatmaps. |

#### Instrument / RH log helpers

| Function | What it does |
|---|---|
| `log_parse(dewpointlog_path, instrument)` | Parses Shasta dewpoint/RH logs and standardizes RH-related columns. |
| `read_log(log_path, dewpointlog_path, dewpointlog_path2, instrument, experiment_start_date)` | Reads Shasta logs with dewpoint/RH information. |
| `log_parse_CX(log_path, start_date)` | Parses CX logs using selected workflow keywords. |
| `read_log_CX(log_path, input_path, output_folder, instrument, start_date)` | Reads and summarizes CX log / TeraTerm style inputs. |
| `create_scatterplot(join_df, x, y, title)` | Generates a simple scatterplot. |

### Usage example: calculate SSPro spec 2 manually

```python
import pandas as pd
import SSPro_func as func

path = r"Z:/Xuan-Z/Shasta/20241209_SSPro_Beta2_NewVol_15cPCR1/241209_SSPro_Beta2"
df = pd.read_csv(path + "/analysis_stats.csv")
df.rename(columns={"Sample": "Well_Type"}, inplace=True)

sample_df = df[df["Well_Type"].isin(["Sample", "sample", "PBMC"])]
spec2_table, *_ = func.sspro_spec2(sample_df, "Barcoded_Reads")

func.plot_table_from_df(spec2_table, "SSPro Spec 2", base_width=2)
```

### Usage example: generate a barcoded-read heatmap

```python
import SSPro_func as func

well_list = func.search_files_in_directory(path)
func.cleanup_barcode(well_list)
well_list["Sample"] = "Beta-2_Raw"

func.gen_heatmap([well_list], "Barcoded_Reads", "Sample", "log", base_width=8)
```

### Usage example: parse a Shasta dewpoint log

```python
import SSPro_func as func

rh_df = func.log_parse(
    r"S:/XuanLi/Shasta_SSPro/147279C_20240403102500_SSPro_Beta2_CellRTSwap/DewPointLog-20240403.tsv",
    instrument="Beta2"
)

print(rh_df.head())
```

---

## 5. `SSPro_report_func.py`

### Purpose

`SSPro_report_func.py` contains report-specific plotting utilities. These functions are used by `SSPro_pipeline_report.py` to save plot images into a report image folder.

### Important functions

| Function | What it does |
|---|---|
| `plotBarcodedReads(df, specific_order, save_path)` | Generates and saves a barcoded-read plot for the first two entries in `specific_order`. |
| `plot(df, target_cols, group_column, specific_order, ...)` | Saves violin-style plots for target columns. |
| `plot_by_fraction_to_html(input_df, target_cols, base_col, group_column, specific_order, ...)` | Saves fraction plots and returns image paths for report use. |
| `gen_heatmap(df_list, param, group_column, colour, filename, base_width=8)` | Saves report-ready heatmap images. |
| `OLD_gen_heatmap(...)` | Older heatmap helper retained for compatibility. |
| `gen_heatmap_individual_test_to_base64(...)` | Generates heatmap output intended for base64/report embedding. |

### Usage example: save report images

```python
import pandas as pd
import SSPro_report_func as reportfunc

image_folder = r"C:/Users/leongs/OneDrive - Takara Bio USA, Inc/0. Projects/Shasta instrument/Hammerhead syringe/SSPro double hammerhead/SSPro_Reports/demo/images/"
specific_order = ["v2_Raw", "v1_Raw", "v2_50K", "v1_50K"]

# df should already contain Sample labels and calculated read fractions
reportfunc.plot(df, ["No_of_Genes"], "Sample", specific_order, save_path=image_folder)
reportfunc.plotBarcodedReads(df, specific_order, save_path=image_folder)
reportfunc.plot_by_fraction_to_html(
    df,
    ["Ribosomal_Reads", "Mitochondrial_Reads", "Intergenic_Reads", "Total_Exon_Reads", "Total_Intron_Reads"],
    "Barcoded_Reads",
    "Sample",
    specific_order,
    save_path=image_folder
)
```

### Usage example: save a heatmap image

```python
reportfunc.gen_heatmap(
    list_of_df,
    "Barcoded_Reads",
    "Sample",
    "log",
    image_folder + "Barcoded.png"
)
```

---

## Common run examples

### Example A: Full Shasta vs CX validation analysis

Use `SSPro_pipeline.py` when you have four datasets:

```text
Shasta raw
CX or baseline raw
Shasta downsampled
CX or baseline downsampled
```

Command:

```bash
cd "C:/Users/leongs/OneDrive - Takara Bio USA, Inc/2. Scripts/SSPro_Script"
python SSPro_pipeline.py
```

### Example B: Generate a report for a validation run

Use `SSPro_pipeline_report.py` when you want a saved HTML report and image folder.

Command:

```bash
cd "C:/Users/leongs/OneDrive - Takara Bio USA, Inc/2. Scripts/SSPro_Script"
python SSPro_pipeline_report.py
```

### Example C: Shasta-only analysis without CX comparator

Use `SSPro_pipeline_nocx.py` when you only have Shasta raw and Shasta downsampled outputs.

Command:

```bash
cd "C:/Users/leongs/OneDrive - Takara Bio USA, Inc/2. Scripts/SSPro_Script"
python SSPro_pipeline_nocx.py
```

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'SSPro_func'`

Run the script from the folder containing `SSPro_func.py`, or add the folder to `PYTHONPATH`.

```bash
cd "C:/Users/leongs/OneDrive - Takara Bio USA, Inc/2. Scripts/SSPro_Script"
python SSPro_pipeline.py
```

### `ModuleNotFoundError: No module named 'parse_CX_log'`

`SSPro_func.py` imports `parse_CX_log`. Make sure `parse_CX_log.py` is in the same folder, or comment out CX log parsing imports/sections if you are not using CX logs.

### `FileNotFoundError` for `analysis_stats.csv`

Check that `path0`, `path1`, `path2`, and `path3` point to folders that contain `analysis_stats.csv`.

### No well list found

`search_files_in_directory()` looks for files containing `WellList` in the filename. Make sure the relevant Mappa/Cogent output folder contains a matching well-list file.

### Heatmaps do not render correctly

Confirm the merged dataframe contains:

```text
Row
Col
Sample
```

and the metric you want to plot, such as:

```text
Barcoded_Reads
No_of_Genes
%Intergenic_Reads
%Ribosomal
%Mitochondrial
```

### Report images are missing

For `SSPro_pipeline_report.py`, confirm:

```python
masterpath = folder1 + "SSPro_Reports/"
image_folder = masterpath + title + "/images/"
```

and verify that the account running Python can create directories in `masterpath`.

---

## Suggested improvements for future maintainability

These scripts are useful for interactive analysis, but they would be easier to reuse if converted into command-line tools. Recommended next upgrades:

1. Add `argparse` support for input folders, labels, output folder, and gene-sensitivity threshold.
2. Move hard-coded paths into a YAML or JSON config file.
3. Add a `requirements.txt` file.
4. Add validation checks for expected input files and required columns.
5. Save all plots to a consistent output folder instead of relying on interactive matplotlib windows.
6. Split analysis sections into reusable functions such as `load_inputs()`, `calculate_specs()`, `make_plots()`, and `write_report()`.
7. Add a small test dataset or mock input folder to confirm that changes do not break the workflow.

---

## Minimal `requirements.txt`

```text
pandas
numpy
matplotlib
seaborn
plotly
mpld3
tabulate
```

Optional:

```text
dataframe_image
```
