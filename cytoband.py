from sys import argv


def get_cytoband(args: list):
    cytoband_file = args[0]
    genome_pos = args[1]


def load_cytoband_data(cytoband_file: str):
    retinfo = {}
    starts = []
    ends = []
    cytobands = []
    with open(cytoband_file, "r") as f:
        old_chr = ""
        for line in f:
            arr = line.strip().split("\t")
            chrom = arr[0]
            if old_chr == "":
                old_chr = chrom
                starts.append(int(arr[1]))
                ends.append(int(arr[2]))
                cytobands.append(arr[3])
            elif chrom == old_chr:
                starts.append(int(arr[1]))
                ends.append(int(arr[2]))
                cytobands.append(arr[3])
            else:
                retinfo[chrom] = [starts, ends, cytobands]
                old_chr = chrom
                starts.clear()
                ends.clear()
                cytobands.clear()


if __name__ == "__main__":
    get_cytoband(argv[1:])
