import argparse
import gzip
from os import mkdir, path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from PIL import Image

from load_loh import extract_loh_data

COLOR_PALETTE = [
    "#db5e56",
    "#56db5e",
    "#5784DB",
    "#DB9457",
    "#57DB94",
    "#9457DB",
    "#57B9DB",
    "#DBC957",
    "#DB5784",
    "#1C1B1B",
]


def plot_cnv():
    args = parse_args()
    plots_dir = args.output_path
    # create directory to store all plots for the given sample
    sample_dir = path.join(plots_dir, args.sample_name)
    if not path.exists(sample_dir):
        mkdir(sample_dir)
    # load CNV data from CNR file data
    print(f"Loading copy number data from {args.cnr_file}...")
    raw_cnv_data = pd.read_csv(args.cnr_file, sep="\t", header=0)
    # filter out antitarget regions
    cnv_data = raw_cnv_data.loc[raw_cnv_data["gene"] != "Antitarget"]
    cnv_data.reset_index(inplace=True, drop=True)
    cnv_data.reset_index(inplace=True)

    output_filename_suffix = "_cnv_loh_plots"

    # Load LOH data
    if args.vcf_file != "no_vcf":
        loh_vars = extract_loh_data(args.vcf_file, args.known_snps_vcf)
        cnames = ["chromosome", "pos", "idx", "af"]
        loh_data = pd.DataFrame(loh_vars, columns=cnames)
        loh_data.reset_index(inplace=True)

        genome_view_file = create_genome_plot(
            sample_dir,
            args.output_format,
            args.sample_name,
            cnv_data=cnv_data,
            loh_data=loh_data,
        )
        output_filename_suffix = "cnv_loh_plots"
    else:
        genome_view_file = create_genome_plot(
            sample_dir,
            args.output_format,
            args.sample_name,
            cnv_data=cnv_data,
        )
        output_filename_suffix = "cnv_plots"
    chr_plot_files = create_chrom_plots(
        cnv_data, sample_dir, args.output_format, args.sample_name
    )

    print("Merging plots into a single PDF file...")
    # Load genome view image and convert to RGB
    genome_image = Image.open(genome_view_file).convert("RGB")
    # Load per chromosome view images and convert to RGB
    chr_images = []
    for f in chr_plot_files:
        chr_images.append(Image.open(f).convert("RGB"))
    output_pdf_file = path.join(
        plots_dir, f"{args.sample_name}_{output_filename_suffix}.pdf"
    )
    genome_image.save(output_pdf_file, save_all=True, append_images=chr_images)
    print(f"PDF file written to {output_pdf_file}")


def create_genome_plot(
    output_path: str,
    file_format: str,
    sample_name: str,
    **plot_datasets,
):
    if "loh_data" in plot_datasets:
        print("Generating genome wide plot for CNV & LOH...")

        plot_meta = get_plot_metadata(sample_name)

        # set seaborn plot style and plot figure dimensions
        sns.set_style("white")
        plt.figure(figsize=(20, 16))

        subplot_count = 0
        for k, v in plot_meta.items():
            subplot_count += 1
            # set color palette count to match unique chromosome numbers
            plot_data = plot_datasets[k]
            chrom_count = plot_data["chromosome"].nunique()
            cpal_len = len(COLOR_PALETTE)
            chr_palette = COLOR_PALETTE * int(chrom_count / cpal_len)
            chr_palette.extend(COLOR_PALETTE[: chrom_count % cpal_len])

            # get the positions for the x axis ticks and corresponding labels
            chromosomes, x_tick_pos, vline_positions = get_chr_x_axis_ticks(plot_data)

            # Chose the first subplot for CNV plot
            plt.subplot(2, 1, subplot_count)
            plot = sns.scatterplot(
                x="index",
                y=v["y_data"],
                hue="chromosome",
                data=plot_data,
                s=5,
                legend=False,
                palette=chr_palette,
            )

            # Set CNV plot title
            plt.title(v["plot_title"], fontsize=30)

            # set custom x axis ticks data, asix title, and rotate labels
            plt.xticks(x_tick_pos, chromosomes)
            plot.set_xticklabels(labels=chromosomes, rotation=0, fontdict={"size": 12})
            plot.set_xlabel("Chromosomes", fontdict={"size": 20})
            plot.set_ylabel(v["y_label"], fontdict={"size": 20})

            # set vertical lines at the chromosome boundaries
            for pos in vline_positions:
                plot.axvline(x=pos, color="black", linewidth=0.6, alpha=0.1)

            # set y axis range dynamically
            ymin = 0
            ymax = 0
            if k == "cnv_data":
                ymin, ymax = get_yaxis_limits(plot_data)
            else:
                # afs are bounded in [0,1] so y (= A_af - B_af) is bounded within [-1,1]
                ymin = 0
                ymax = 1
            plot.set(ylim=(ymin, ymax))

        # adjust horizontal space
        plt.subplots_adjust(hspace=0.4)
        output_file = path.join(
            output_path, f"{sample_name}_cnv_loh_genome_view.{file_format}"
        )
        plt.savefig(output_file, dpi=300)
        return output_file
    else:
        print("Generating genome wide plot for CNV only...")
        cnv_data = plot_datasets["cnv_data"]
        # set figure size (w, h)
        # sns.set(rc={"figure.figsize": (20, 8)})
        sns.set_theme(rc={"figure.figsize": (20, 8)})
        # plot theme
        sns.set_style("white")
        # set color palette count to match unique chromosome numbers
        chrom_count = cnv_data["chromosome"].nunique()
        cpal_len = len(COLOR_PALETTE)
        chr_palette = COLOR_PALETTE * int(chrom_count / cpal_len)
        chr_palette.extend(COLOR_PALETTE[: chrom_count % cpal_len])
        # get the positions for the x axis ticks and corresponding labels
        chromosomes, x_tick_pos, vline_positions = get_chr_x_axis_ticks(cnv_data)
        # create chart axes
        ax = sns.scatterplot(
            x="index",
            y="log2",
            hue="chromosome",
            data=cnv_data,
            s=5,
            legend=False,
            palette=chr_palette,
        )
        # set custom x axis ticks data and rotate labels
        plt.xticks(x_tick_pos, chromosomes)
        ax.set_xticklabels(labels=chromosomes, rotation=0)
        # set vertical lines at the chromosome boundaries
        for pos in vline_positions:
            ax.axvline(x=pos, color="black", linewidth=0.6, alpha=0.1)
        # set y axis range dynamically
        ymin, ymax = get_yaxis_limits(cnv_data)
        ax.set(ylim=(ymin, ymax))
        # set chart title
        ax.set(title=f"{sample_name} - Genome View")
        output_file = path.join(
            output_path, f"{sample_name}_all_genome_view.{file_format}"
        )
        plt.savefig(output_file, dpi=300)
        return output_file


def get_yaxis_limits(cnv_data: pd.DataFrame) -> tuple:
    ymin = cnv_data["log2"].min()
    if ymin > -3:
        ymin = -3
    else:
        ymin -= 1
    ymax = cnv_data["log2"].max()
    if ymax < 4:
        ymax = 4
    else:
        ymax += 1

    return ymin, ymax


def create_chrom_plots(
    cnv_data: pd.DataFrame, output_path: str, file_format: str, sample_name: str
):
    print("Generating CNV plots per chromosome...")
    chromosomes = list(cnv_data["chromosome"].drop_duplicates())
    chr_plot_files = []
    print("Loading gene chromosome cytoband info...")
    gene_cytobands = get_cytoband_transitions("gene_cytoband.txt.gz")
    for chromosome in chromosomes:
        plt.clf()  # clear prior figure
        # set figure size (w, h)
        # sns.set_theme(rc={"figure.figsize": (20, 8)})
        plt.figure(figsize=(20, 10))
        # plot theme
        sns.set_style("white")
        # get CNV data for specific chromosome
        chr_cnv = cnv_data.loc[cnv_data["chromosome"] == chromosome]
        gene_count = chr_cnv["gene"].nunique()
        # set color palatte for each gene sequenced in a chromosome
        cpal_len = len(COLOR_PALETTE)
        c_palette = COLOR_PALETTE * int(gene_count / cpal_len)
        c_palette.extend(COLOR_PALETTE[: gene_count % cpal_len])
        # get the plot axes
        ax = sns.scatterplot(
            x="index",
            y="log2",
            hue="gene",
            data=chr_cnv,
            s=20,
            legend=False,
            palette=c_palette,
        )
        # set y axis range
        ymin = chr_cnv["log2"].min()
        if ymin > -3:
            ymin = -3
        else:
            ymin -= 1
        ymax = chr_cnv["log2"].max()
        if ymax < 4:
            ymax = 4
        else:
            ymax += 1
        ax.set(ylim=(ymin, ymax))

        # set the x-axis gene list and index positions
        genes, x_tick_pos, gene_vline_positions, pq_trans_idx = get_gene_x_axis_ticks(
            chr_cnv, gene_cytobands
        )
        for pos in gene_vline_positions:
            if pos == pq_trans_idx:
                # set the vertical line indicating p arm to q arm transition
                ax.axvline(x=pos, color="blue", linewidth=0.6, alpha=0.5)
            else:
                ax.axvline(x=pos, color="black", linewidth=0.6, alpha=0.1)
        plt.xticks(x_tick_pos, genes)
        # set axis label title
        ax.set_xlabel("Gene", fontdict={"size": 20})
        ax.set_ylabel("Log2 ratio", fontdict={"size": 20})
        # rotate labels to 90 degree
        ax.set_xticklabels(labels=genes, rotation=90)

        # set labels for p and q arms of chromosomes
        if pq_trans_idx > gene_vline_positions[0]:
            ax.text(pq_trans_idx - 3, ymax - 0.3, "<- p arm", fontsize=14, ha="right")
        ax.text(pq_trans_idx + 3, ymax - 0.3, "q arm ->", fontsize=14, ha="left")
        # set chart title
        plt.title(f"{sample_name} - {chromosome}", fontsize=20)
        output_file = path.join(
            output_path, f"{sample_name}_{chromosome}.{file_format}"
        )
        plt.savefig(output_file, dpi=300)
        chr_plot_files.append(output_file)
        print(f"Generated plot for {chromosome}")
        plt.close()
    return chr_plot_files


def get_chr_x_axis_ticks(cnv_data: pd.DataFrame):
    chromosomes = list(cnv_data["chromosome"].drop_duplicates())
    wg_x_tick_points = []
    chr_boundaries = []
    end_idx = 0
    for chromosome in chromosomes:
        chr_cnv_idx = cnv_data.loc[cnv_data["chromosome"] == chromosome]["index"]
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
    # remove chr prefix for plotting
    chromosomes = [c.strip("chr") for c in chromosomes]
    return chromosomes, wg_x_tick_points, chr_boundaries


def get_gene_x_axis_ticks(chr_cnv_data: pd.DataFrame, gene_cytobands: dict):
    genes = list(chr_cnv_data["gene"].drop_duplicates())
    chr_x_tick_points = []
    gene_boundaries = []
    end_idx = 0
    pq_trans_idx = -1
    for gene in genes:
        gene_cnv_idx = chr_cnv_data.loc[chr_cnv_data["gene"] == gene]["index"]
        start_idx = list(gene_cnv_idx[:1])[0]
        end_idx = list(gene_cnv_idx[-1:])[0]
        gene_boundaries.append(start_idx)

        # find the gene index position where p to q transition occurs
        if pq_trans_idx == -1:
            cytoband = gene_cytobands.get(gene, "")
            if cytoband != "":
                if "q" in cytoband:
                    pq_trans_idx = start_idx

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
    return genes, chr_x_tick_points, gene_boundaries, pq_trans_idx


def get_cytoband_transitions(cytoband_file: str) -> dict:
    retInfo = {}
    with gzip.open(cytoband_file, "r") as f:
        for byteline in f:
            line = byteline.decode("UTF-8")
            gene, cytoband = line.strip().split("\t")
            retInfo[gene] = cytoband
    return retInfo


def get_plot_metadata(sample_name: str) -> dict:
    return {
        "cnv_data": {
            "y_data": "log2",
            "y_label": "Log2 ratio",
            "plot_title": f"{sample_name} - CNV plot",
        },
        "loh_data": {
            "y_data": "af",
            "y_label": "B-allele fraction",
            "plot_title": f"{sample_name} - LOH plot",
        },
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--cnr-file",
        type=check_path,
        required=True,
        help="CNR file containing weighted log2 ratio info.",
    )
    parser.add_argument(
        "-v",
        "--vcf-file",
        type=check_path,
        required=False,
        default="no_vcf",
        help="Input VCF file. Should be normalized if possible. Variant annotation or HGVS nomenclature is not required.",
    )
    parser.add_argument(
        "-k",
        "--known-snps-vcf",
        type=check_path,
        required=False,
        default="no_vcf",
        help="VCF file containing the known SNPs for plotting LOH. Should be normalized if possible. Variant annotation or HGVS nomenclature is not required.",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        type=check_path,
        required=True,
        help="Output folder to save plot images. This folder must exist.",
    )
    parser.add_argument(
        "-f",
        "--output-format",
        type=acceptable_formats,
        default="png",
        required=True,
        help="Output file format. Supported types: png, jpg, tiff, pdf, svg. Default is png.",
    )
    parser.add_argument(
        "-s",
        "--sample-name",
        type=str,
        default="sample",
        required=True,
        help="Sample name to include in the chart title",
    )
    return parser.parse_args()


def acceptable_formats(format: str):
    formats = ["png", "jpg", "tiff", "svg", "pdf"]
    if format in formats:
        return format
    else:
        raise argparse.ArgumentTypeError(
            f"{format} is not a acceptable file format. Allowed types: png, jpg, tiff, pdf, svg."
        )


def check_path(file_path: str):
    if file_path != "no_vcf":
        if not path.exists(file_path):
            raise argparse.ArgumentTypeError(
                "Error: Provided directory does not exist."
            )
    return file_path


if __name__ == "__main__":
    plot_cnv()
