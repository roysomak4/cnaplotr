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
    cnames = ['chromosome', 'pos', 'idx', 'af']
    cnv_data = pd.read_csv(args.cnr_file, delim_whitespace=True, names=cnames)
    print(cnv_data.keys())
    # transform bafs to heterozygosity
    # we want A_af - B_af but A_af + B_af = 1, so this is just (1 - B_af) - B_af)
    # = 1 - 2 * B_af
    cnv_data.reset_index(inplace=True)

    genome_view_file = create_genome_plot(cnv_data, sample_dir, args.output_format, args.sample_name)

    print("Merging plots into a single PDF file...")
    # Load genome view image and convert to RGB
    genome_image = Image.open(genome_view_file).convert("RGB")
    output_pdf_file = path.join(plots_dir, f"{args.sample_name}_cnv_plots.pdf")
    genome_image.save(output_pdf_file, save_all=True)
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
    ax = sns.scatterplot(x = 'index', y = 'af', hue = 'chromosome', data = cnv_data, s = 5, legend = False, palette=chr_palette)
    # set custom x axis ticks data and rotate labels
    plt.xticks(x_tick_pos, chromosomes)
    ax.set_xticklabels(labels = chromosomes, rotation=90)
    # set vertical lines at the chromosome boundaries
    for pos in vline_positions:
        ax.axvline(x=pos, color='black', linewidth=0.6, alpha=0.1)
    # afs are bounded in [0,1] so y (= A_af - B_af) is bounded within [-1,1]
    ymin = 0
    ymax = 1
    ax.set(ylim=(ymin, ymax))
    # set chart title
    ax.set(title=f"LOH: {sample_name} - Genome View")
    output_file = path.join(output_path, f"{sample_name}_all_genome_view.{file_format}")
    plt.savefig(output_file, dpi=300)
    return output_file

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

