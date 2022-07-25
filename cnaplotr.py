import argparse
from this import d
import pandas as pd
from matplotlib import pyplot as plt
from os import path
import seaborn as sns

COLOR_PALETTE = ['#db5e56', '#56db5e', '#5784DB', '#DB9457', '#57DB94', '#9457DB', '#57B9DB', '#DBC957', '#DB5784', '#1C1B1B']

def plot_cnv():
    args = parse_args()
    # load CNR file data
    raw_cnv_data = pd.read_csv(args.cnr_file, sep='\t', header=0)
    # filter out antitarget regions
    cnv_data = raw_cnv_data.loc[raw_cnv_data['gene'] != 'Antitarget']
    cnv_data.reset_index(inplace=True, drop=True)
    cnv_data.reset_index(inplace=True)

def get_chr_x_axis_ticks(cnv_data: pd.DataFrame):
    chromosomes = list(cnv_data['chromosome'].drop_duplicates())
    wg_x_tick_points = []
    chr_boundaries = []
    end_idx = 0
    for chromosome in chromosomes:
        chr_cnv_idx = cnv_data.loc[cnv_data['chromosome'] == chromosome]['index']
        start_idx = list(chr_cnv_idx[:1])[0]
        end_idx = list(chr_cnv_idx[-1:])[0]
        chr_boundaries.append(start_idx)
        idx_range = chr_cnv_idx.count()
        if idx_range > 0:
            mid_point = 0
            if idx_range % 2 == 0:
                mid_point = int(idx_range / 2)
            else:
                mid_point = int(idx_range / 2 + 0.5)
            mid_idx = start_idx + mid_point
            wg_x_tick_points.append(mid_idx)
    chr_boundaries.append(end_idx)
    return chromosomes, wg_x_tick_points, chr_boundaries    

def get_gene_x_axis_ticks(chr_cnv_data: pd.DataFrame):
    genes = list(chr_cnv_data['gene'].drop_duplicates())
    chr_x_tick_points = []
    gene_boundaries = []
    end_idx = 0
    for gene in genes:
        gene_cnv_idx = chr_cnv_data.loc[chr_cnv_data['gene'] == gene]['index']
        start_idx = list(gene_cnv_idx[:1])[0]
        end_idx = list(gene_cnv_idx[-1:])[0]
        gene_boundaries.append(start_idx)
        idx_range = gene_cnv_idx.count()
        if idx_range > 0:
            mid_point = 0
            if idx_range % 2 == 0:
                mid_point = int(idx_range / 2)
            else:
                mid_point = int(idx_range / 2 + 0.5)
            mid_idx = start_idx + mid_point
            chr_x_tick_points.append(mid_idx)
    gene_boundaries.append(end_idx)
    return genes, chr_x_tick_points, gene_boundaries

def create_chr_plot(cnv_data: pd.DataFrame, output_path: str, filetype: str):
    # set figure size (w, h)
    sns.set(rc={
        "figure.figsize": (25, 8)
    })
    # plot theme
    sns.set_style("white")
    # set color palette count to match unique chromosome numbers
    chrom_count = cnv_data["chromosome"].nunique()
    cpal_len = len(COLOR_PALETTE)
    chr_palette = COLOR_PALETTE * int(chrom_count / cpal_len)
    chr_palette.extend(COLOR_PALETTE[:chrom_count % cpal_len])
    # get the positions for the x axis ticks and corresponding labels
    chromosomes, x_tick_pos, vline_positions = get_chr_x_axis_ticks(cnv_data)
    # create chart axes
    ax = sns.scatterplot(x = 'index', y = 'log2', hue = 'chromosome', data = cnv_data, s = 5, legend = False, palette=chr_palette)
    # set custom x axis ticks data and rotate labels
    plt.xticks(x_tick_pos, chromosomes)
    ax.set_xticklabels(labels = chromosomes, rotation=90)
    # set vertical lines at the chromosome boundaries
    for pos in vline_positions:
        ax.axvline(x=pos, color='black', linewidth=0.6, alpha=0.1)
    # set y axis range dynamically
    ymin = -3
    ymax = cnv_data['log2'].max()
    if ymax < 4:
        ymax = 4
    else:
        ymax += 1
    ax.set(ylim=(ymin, ymax))
    output_file = path.join(output_path, "")
    plt.savefig(f'{output_path}', dpi=300)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("cnr_file",
                        type=str,
                        help="CNR file containing weighted log2 ratio info.")
    parser.add_argument("output",
                        type=str,
                        help="Output folder to save plot images.")
    parser.add_argument("--output-format",
                        type=str,
                        default="png",
                        help="Output file format. Supported types: PNG, JPG, TIFF, PDF, SVG. Default is PNG.")
    return parser.parse_args()


if __name__ == '__main__':
    plot_cnv()
    