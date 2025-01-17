import argparse
import re

parser = argparse.ArgumentParser(description="Extract PVE from GEMMA LMM logs.")
parser.add_argument("logs", nargs="+", help="One or more files containing concatenated logs.")
parser.add_argument("-o", "--output", help="Output file (TSV)")
# parser.add_argument("--vc_mode", action="store_true", default=False, help="Logs are from vc mode (PVE for GRM only).")
args = parser.parse_args()

with open(args.output, "w") as out:
    out.write("gene_id\tpve\tse_pve\n")
    for logfile in args.logs:
        with open(logfile, "r") as f:
            for line in f:
                m = re.search("-o (\S+) ", line)
                if m:
                    gene = m.group(1)
                # if args.vc_mode:
                #     m = re.search("pve estimates =\s+(\S+)\n$", line)
                # else:
                #     m = re.search("pve estimate in the null model = (.+)\n$", line)
                m = re.search("pve estimate.+\s(\S+)\n$", line)
                if m:
                    pve = m.group(1)
                # if args.vc_mode:
                #     m = re.search("se\(pve\) =\s+(\S+)\n$", line)
                # else:
                #     m = re.search("se\(pve\) in the null model = (.+)\n$", line)
                m = re.search("se\(pve\).+\s(\S+)\n$", line)
                if m:
                    se_pve = m.group(1)
                    out.write(f"{gene}\t{pve}\t{se_pve}\n")
                    # Clear variables in case of parsing issue later:
                    gene, pve, se_pve = None, None, None
