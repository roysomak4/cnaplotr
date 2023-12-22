import argparse
from os import listdir, path, remove

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from PIL import Image


def covplot() -> None:
    args = parse_args()

    # set output filenames
    output_pdf_file = path.join(
        args.output_path, f"{args.sample_name}_gene_coverage_plot.pdf"
    )
    output_cov_file = path.join(
        args.output_path, f"{args.sample_name}_low_coverage_genes.txt"
    )
    # load coverage data into pandas dataframe
    print("Loading BED coverage data...")
    column_names = ["chrom", "start", "stop", "exon_info", "coverage"]
    cov_data = pd.read_csv(
        args.coverage_bedfile,
        compression="gzip",
        header=0,
        sep="\t",
        names=column_names,
    )
    # split the 4th column of bed file into individual columns
    cov_data[["gene", "tx", "exon"]] = cov_data["exon_info"].str.split(
        ";", n=3, expand=True
    )
    cov_data = cov_data.sort_values(by=["gene"], ascending=True)

    # get list of unique genes
    print("Listing unique genes for coverage calculation...")
    genes = list(cov_data["gene"].unique())

    # split coverage data in set of 60 genes
    cov_dfs = []
    chunk_size = 60
    gene_chunks = chunk_list(genes, chunk_size)
    for gene_list in gene_chunks:
        tmp_df = cov_data[cov_data["gene"].isin(gene_list)]
        cov_dfs.append([tmp_df, gene_list])

    # generate plots for each set of genes
    print("Generating intermediate plot files...")
    iteration = 0
    low_cov_genes = []
    for cov_df in cov_dfs:
        iteration += 1
        low_cov_genes.extend(generate_low_cov_gene_list(cov_df))
        generate_plot(
            cov_df, args.sample_name, args.output_path, args.output_format, iteration
        )

    print("Generate PDF plot file...")
    generate_pdf(output_pdf_file, args.sample_name, args.output_path)
    print(f"Plot saved at {output_pdf_file}")

    print("Writing low coverage gene list...")
    write_low_cov_genes(low_cov_genes, output_cov_file)
    print(f"Data saved at {output_cov_file}")


def generate_low_cov_gene_list(cov_df: list) -> list:
    low_cov_genes = []
    gene_cov_df = cov_df[0][["gene", "coverage"]]
    gene_cov_df = (
        gene_cov_df.groupby("gene")["coverage"].mean().reset_index(name="coverage")
    )
    gene_cov_df = gene_cov_df.astype({"coverage": "int"})
    low_cov_genes = gene_cov_df.loc[gene_cov_df["coverage"] < 300]["gene"].tolist()
    return low_cov_genes


def write_low_cov_genes(low_cov_genes: list, output_file: str) -> None:
    outbuffer = "None\n"
    if len(low_cov_genes) > 0:
        outbuffer = ", ".join(low_cov_genes) + "\n"
    with open(output_file, "w") as of:
        of.write(outbuffer)


def generate_pdf(output_pdf_file: str, sample_name: str, output_path: str) -> None:
    # find list of png files in the output folder
    png_files = [
        path.join(output_path, x)
        for x in listdir(output_path)
        if (x.endswith(".png") and sample_name in x)
    ]
    png_files.sort()
    output_pdf_hdl = Image.open(png_files[0]).convert("RGB")
    cov_images = []
    for image_file in png_files[1:]:
        cov_images.append(Image.open(image_file).convert("RGB"))
    output_pdf_hdl.save(output_pdf_file, save_all=True, append_images=cov_images)
    # delete the intermediate png files
    for image_file in png_files:
        remove(image_file)


def generate_plot(
    cov_df: list, sample_name: str, output_path: str, file_type: str, iteration: int
) -> None:
    plt.figure(figsize=(30, 17.14))
    plt.rcParams["axes.labelsize"] = 22
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16
    ax = sns.barplot(x="gene", y="coverage", data=cov_df[0], errorbar=("ci", False))
    ax.set_xticklabels(labels=cov_df[1], rotation=90)
    ax.axhline(300)
    plt.xlabel("Gene")
    plt.ylabel("Coverage")
    plt.title(f"Per gene coverage - {sample_name}", fontsize=40, pad=30)
    plt.subplots_adjust(bottom=0.2)
    output_file = path.join(
        output_path, f"{sample_name}_gene_coverage_set{iteration}.{file_type}"
    )
    plt.savefig(output_file, dpi=300)
    plt.close()


def chunk_list(input_list: list, chunk_size: int) -> list:
    retinfo = []
    counter = 0
    tmp_list = []
    for item in input_list:
        counter += 1
        if counter > chunk_size:
            counter = 0
            retinfo.append(tmp_list)
            tmp_list = []
        tmp_list.append(item)
    if len(tmp_list) > 0:
        retinfo.append(tmp_list)
    return retinfo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--coverage-bedfile",
        type=check_path,
        required=True,
        help="Coverage file generated by mosdepth.",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        type=check_path,
        required=True,
        help="Output folder to save coverage plot images. This folder must exist.",
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


def check_path(file_path: str) -> str:
    if not path.exists(file_path):
        raise argparse.ArgumentTypeError("Error: Provided directory does not exist.")
    else:
        return file_path


if __name__ == "__main__":
    covplot()
