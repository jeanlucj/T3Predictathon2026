import gzip
import re

def get_samples(vcf_path):
    """
    Extract sample names from a VCF header.
    Handles:
      - tabs or space-delimited headers
      - multi-word sample names (e.g., 'SD ANDES')
      - missing FORMAT column (Tassel bug)
      - gzipped or plain VCFs
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt", errors="ignore") as f:
        for raw in f:
            if raw.startswith("#CHROM"):
                line = raw.rstrip("\n")

                # Case 1: tab-delimited (correct)
                if "\t" in line:
                    fields = line.split("\t")

                else:
                    # Case 2: space-delimited (broken)
                    # Split only on *2 or more* spaces, preserving single spaces inside sample names
                    fields = re.split(r" {2,}", line)

                # Fix Tassel bug: missing FORMAT column
                if len(fields) >= 8:
                    if fields[7].upper().startswith("INFO") and (
                        len(fields) < 9 or not fields[8].upper().startswith("FORMAT")
                    ):
                        fields.insert(8, "FORMAT")

                # Samples start at column 9
                return fields[9:]

    return []