from pysam import VariantFile
from sys import argv


def extract_loh_data(input_vcf: str, known_snps_vcf: str) -> list:
    # iterate through snps in the known snps vcf and query in the input vcf
    loh_vars = []
    var_idx = 0
    total_vars = 0
    soi = 0  # SNPs of interest
    print("Loading and searching SNP info for LOH analysis...")
    with VariantFile(known_snps_vcf) as ks:
        with VariantFile(input_vcf) as invcf:
            for snp in ks.fetch():
                total_vars += 1
                search_space = f"{snp.chrom}:{snp.pos}-{snp.pos}"
                for variant in invcf.fetch(region=search_space):
                    if variant:
                        tcov = int(variant.info["DP"])
                        if tcov >= 50:
                            soi += 1
                            var_info = []
                            var_info.append(variant.chrom)
                            var_info.append(variant.pos)
                            var_info.append(var_idx)
                            var_info.append(variant.info["AF"][0])
                            loh_vars.append(var_info)
                            var_idx += 1
    print("LOH SNP lookup complete.")
    print(f"Total SNPs queried: {total_vars}")
    print(f"SNPs identified:    {soi}")
    return loh_vars


if __name__ == "__main__":
    input_vcf = argv[1]
    known_snp_vcf = argv[2]
    loh_variants = extract_loh_data(input_vcf, known_snp_vcf)
    for x in range(6):
        print(loh_variants[x])
