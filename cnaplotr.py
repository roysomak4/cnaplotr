import argparse
import pandas as pd
from matplotlib import pyplot as plt
from os import path, mkdir
from PIL import Image
import seaborn as sns

COLOR_PALETTE = ['#db5e56', '#56db5e', '#5784DB', '#DB9457', '#57DB94', '#9457DB', '#57B9DB', '#DBC957', '#DB5784', '#1C1B1B']

def plot_cnv():
    args = parse_args()
    plots_dir = args.output_path
    # create directory to store all plots for the given sample
    sample_dir = path.join(plots_dir, args.sample_name)
    if not path.exists(sample_dir):
        mkdir(sample_dir)
    # load CNR file data
    print(f"Loading copy number data from {args.cnr_file}...")
    raw_cnv_data = pd.read_csv(args.cnr_file, sep='\t', header=0)
    # filter out antitarget regions
    cnv_data = raw_cnv_data.loc[raw_cnv_data['gene'] != 'Antitarget']
    cnv_data.reset_index(inplace=True, drop=True)
    cnv_data.reset_index(inplace=True)

    genome_view_file = create_genome_plot(cnv_data, sample_dir, args.output_format, args.sample_name)
    chr_plot_files = create_chrom_plots(cnv_data, sample_dir, args.output_format, args.sample_name)
    
    print("Merging plots into a single PDF file...")
    # Load genome view image and convert to RGB
    genome_image = Image.open(genome_view_file).convert("RGB")
    # Load per chromosome view images and convert to RGB
    chr_images = []
    for f in chr_plot_files:
        chr_images.append(Image.open(f).convert("RGB"))
    output_pdf_file = path.join(plots_dir, f"{args.sample_name}_cnv_plots.pdf")
    genome_image.save(output_pdf_file, save_all=True, append_images=chr_images)
    print(f"PDF file written to {output_pdf_file}")

def create_genome_plot(cnv_data: pd.DataFrame, output_path: str, file_format: str, sample_name: str):
    print("Generating genome wide plot...")
    # set figure size (w, h)
    sns.set(rc={
        "figure.figsize": (20, 8)
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
    # set chart title
    ax.set(title=f"{sample_name} - Genome View")
    output_file = path.join(output_path, f"{sample_name}_all_genome_view.{file_format}")
    plt.savefig(output_file, dpi=300)
    return output_file

def create_chrom_plots(cnv_data: pd.DataFrame, output_path: str, file_format: str, sample_name: str):
    print("Generating per chromosome plots...")
    chromosomes = list(cnv_data['chromosome'].drop_duplicates())
    chr_plot_files = []
    for chromosome in chromosomes:
        plt.clf() # clear prior figure
        # set figure size (w, h)
        sns.set(rc={
            "figure.figsize": (20, 8)
        })
        # plot theme
        sns.set_style("white")
        # get CNV data for specific chromosome
        chr_cnv = cnv_data.loc[cnv_data['chromosome'] == chromosome]
        gene_count = chr_cnv['gene'].nunique()
        # set color palatte for each gene sequenced in a chromosome
        cpal_len = len(COLOR_PALETTE)
        c_palette = COLOR_PALETTE * int(gene_count / cpal_len)
        c_palette.extend(COLOR_PALETTE[:gene_count % cpal_len])
        # get the plot axes
        ax = sns.scatterplot(x = 'index', y = 'log2', hue = 'gene', data = chr_cnv, s = 20, legend = False, palette=c_palette)
        # set the x-axis gene list and index positions
        genes, x_tick_pos, gene_vline_positions = get_gene_x_axis_ticks(chr_cnv)
        for pos in gene_vline_positions:
            ax.axvline(x=pos, color='black', linewidth=0.6, alpha=0.1)
        plt.xticks(x_tick_pos, genes)
        # rotate labels to 90 degree
        ax.set_xticklabels(labels = genes, rotation=90)
        # set y axis range
        ymin = -3
        ymax = chr_cnv['log2'].max()
        if ymax < 4:
            ymax = 4
        else:
            ymax += 1
        ax.set(ylim=(ymin, ymax))
        # set chart title
        ax.set(title=f"{sample_name} - {chromosome}")
        output_file = path.join(output_path, f"{sample_name}_{chromosome}.{file_format}")
        plt.savefig(output_file, dpi=300)
        chr_plot_files.append(output_file)
        print(f"Generated plot for {chromosome}")
    return chr_plot_files

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

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--cnr-file",
                        type=check_path,
                        required=True,
                        help="CNR file containing weighted log2 ratio info.")
    parser.add_argument("-o", "--output-path",
                        type=check_path,
                        required=True,
                        help="Output folder to save plot images. This folder must exist.")
    parser.add_argument("-f", "--output-format",
                        type=acceptable_formats,
                        default="png",
                        required=True,
                        help="Output file format. Supported types: png, jpg, tiff, pdf, svg. Default is png.")
    parser.add_argument("-s", "--sample-name",
                        type=str,
                        default="sample",
                        required=True,
                        help="Sample name to include in the chart title")
    return parser.parse_args()

def acceptable_formats(format: str):
    formats = ['png', 'jpg', 'tiff', 'svg', 'pdf']
    if format in formats:
        return format
    else:
        raise argparse.argparse.ArgumentTypeError(f'{format} is not a acceptable file format. Allowed types: png, jpg, tiff, pdf, svg.')

def check_path(file_path: str):
    if not path.exists(file_path):
        raise argparse.ArgumentTypeError("Error: Provided output directory does not exist.")
    else:
        return file_path

if __name__ == '__main__':
    plot_cnv()
    