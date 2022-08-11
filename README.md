# cnaplotr
Python package for generating copy number alteration (CNA) plots

<br>
<hr>

## Dependencies
The program requires python v3.8 or higher. The following python packages are required for this program
1. pandas - to read the cnr file
2. matplotlib and seaborn - generating the CNA plots.
The packages can be installed using pip
```
pip install pandas, matplotlib, seaborn
```
At the time of development the following version of the libraries were tested.
1. `pandas v1.4.3`
2. `matplotlib 3.5.2`
3. `seaborn 0.11.2`
<br>
<hr>

## Installation instruction
Installation is simple. Just download the python program to any directory of your chosing.
```
wget https://raw.githubusercontent.com/roysomak4/cnaplotr/main/cnaplotr.py
```
For convenience, a docker image of cnaplotr is available at <a href="https://gitlab.com/roysomak4/genbio_containers/container_registry/3271241" target="_blank">Gitlab container registry</a>. Use the following commands to download the image and run cnaplotr in a Docker container.
```
docker pull registry.gitlab.com/roysomak4/genbio_containers/cnaplotr:v0.2-bullseye

docker run --rm -v /path/to/cnvkit_cnr_files:/data registry.gitlab.com/roysomak4/genbio_containers/cnaplotr:v0.2-bullseye bash -c "python3 cnaplotr.py --cnr-file /data/sample1.cnr --output-path /data/output_folder --output-format png --sample-name sample1"
```
Please update the path to the data folder on the host machine and the appropriate sample name when running on your system.
<br>
<hr>

## How to run the program
and execute on the command line like such
```
python3 cnaplotr.py \
        --cnr-file sample1.cnr \
        --output-path /path/to/result_folder \
        --output-format png \
        --sample-name sample1
```

The details of all options can be viewed using the `-h` or `--help` flag
```
python3 cnaplotr.py -h                                                  

usage: cnaplotr.py [-h] -i CNR_FILE -o OUTPUT_PATH -f OUTPUT_FORMAT -s
                   SAMPLE_NAME

options:
  -h, --help            show this help message and exit
  -i CNR_FILE, --cnr-file CNR_FILE
                        CNR file containing weighted log2 ratio info.
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
                        Output folder to save plot images. This folder must
                        exist. A 'plots' folder will be created inside the
                        output path folder.
  -f OUTPUT_FORMAT, --output-format OUTPUT_FORMAT
                        Output file format. Supported types: png, jpg, tiff,
                        pdf, svg. Default is png.
  -s SAMPLE_NAME, --sample-name SAMPLE_NAME
                        Sample name to include in the chart title
```

## Output CNA Plot
The python program generates the following output
1. Image of a whole genome or all chromosome view plot. Depending on the targets sequenced, the plot will show only those chromosomes with sequence data. For exome sequencing or most comprehensive NGS panels, data for almost all or all chromosomes will be displayed.

Example of a whole genome CNA plot of a normal non-FFPE sample (NIST HG-003)
![Whole genome view normal non-FFPE sample](https://github.com/roysomak4/cnaplotr/raw/main/plots/Normal_sample_whole_genome_view.png)

Example of a whole genome CNA plot for a tumor sample (FFPE)
![Whole genome view tumor FFPE sample](https://github.com/roysomak4/cnaplotr/raw/main/plots/Tumor_WG_FFPE.png)

2. Per chromosome plot.

![Whole genome view tumor FFPE sample](https://github.com/roysomak4/cnaplotr/raw/main/plots/Tumor_FFPE_chr17.png)


## Questions / issues / feature requests
This is a beta release of this program. Feedback for bugs and feature requests are most welcome. Please use GitHub issues for putting in your questions. 