import os
import csv
import urllib.request
import json
import time

# Directories
WORKDIR = "/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline"
READS_DIR = os.path.join(WORKDIR, "reads")
RESULTS_DIR = os.path.join(WORKDIR, "results")
os.makedirs(READS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

raw_file_path = os.path.join(READS_DIR, "23andme.txt")
raw_url = "https://9beeffb1f3ae3b4b0290c3a9935fb460-105.collections.ac2it.arvadosapi.com/_/genome_Lorena_Sandoval_v5_Full_20260429131650.txt"

# 1. Download raw 23andMe data file if it does not already exist
if not os.path.exists(raw_file_path):
    print(f"Downloading 23andMe raw data from {raw_url}...")
    req = urllib.request.Request(raw_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as r, open(raw_file_path, 'wb') as f:
        f.write(r.read())
    print("Download complete.")
else:
    print(f"Using cached raw 23andMe file: {raw_file_path}")

# 2. Parse 23andMe genotyping data into a fast-lookup dictionary: rsid -> genotype
print("Parsing 23andMe genotypes...")
genotype_map = {}
with open(raw_file_path, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        # Skip comments
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
            genotype_map[rsid] = genotype

print(f"Loaded {len(genotype_map)} genotyped markers.")

# 3. Define target SNPs, associated genes, and functional profiles for analysis
snp_targets = {
    "rs1801133": {
        "Gene": "MTHFR",
        "Trait": "Folate metabolism (C677T)",
        "Risk_Allele": "T",
        "Notes": "TT: ~16% folate conversion; CT: ~65%; CC: Typical"
    },
    "rs1801131": {
        "Gene": "MTHFR",
        "Trait": "Folate metabolism (A1298C)",
        "Risk_Allele": "C",
        "Notes": "CC: Reduced folate conversion; AC: Mildly reduced"
    },
    "rs429358": {
        "Gene": "APOE",
        "Trait": "Alzheimer's risk (C112R)",
        "Risk_Allele": "C",
        "Notes": "C allele determines APOE-epsilon4 (Increased Alzheimer's/lipid risk)"
    },
    "rs7412": {
        "Gene": "APOE",
        "Trait": "Alzheimer's risk (R158C)",
        "Risk_Allele": "T",
        "Notes": "T allele determines APOE-epsilon2 (Reduced Alzheimer's risk)"
    },
    "rs671": {
        "Gene": "ALDH2",
        "Trait": "Alcohol flush reaction (Glu504Lys)",
        "Risk_Allele": "A",
        "Notes": "AA: Severe flush/hangover; GA: Mild flush; GG: Normal tolerance"
    },
    "rs1815739": {
        "Gene": "ACTN3",
        "Trait": "Muscle performance (R577X)",
        "Risk_Allele": "T",
        "Notes": "CC: Sprinters/Power athletes; CT: Intermediate; TT: Endurance athletes"
    },
    "rs9939609": {
        "Gene": "FTO",
        "Trait": "Weight and obesity risk",
        "Risk_Allele": "A",
        "Notes": "AA: Increased obesity risk; TA: Intermediate; TT: Typical risk"
    },
    "rs4988235": {
        "Gene": "LCT",
        "Trait": "Lactase persistence (Lactose tolerance)",
        "Risk_Allele": "G",
        "Notes": "AA/AG: Lactose tolerant; GG: Lactose intolerant"
    },
    "rs1800562": {
        "Gene": "HFE",
        "Trait": "Hemochromatosis (C282Y)",
        "Risk_Allele": "A",
        "Notes": "AA: High risk of iron overload; GA: Carrier; GG: Typical"
    },
    "rs1805007": {
        "Gene": "MC1R",
        "Trait": "Red hair & fair skin",
        "Risk_Allele": "T",
        "Notes": "TT: High probability of red hair; CT: Carrier; CC: Typical"
    }
}

# 4. Query Ensembl REST API for dynamic metadata retrieval (consequence, coordinates)
def get_ensembl_metadata(rsid):
    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?content-type=application/json"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req) as r:
            res = json.loads(r.read().decode('utf-8'))
            if res and len(res) > 0:
                consequence = res[0].get("most_severe_consequence", "unknown")
                chrom = res[0].get("seq_region_name", "unknown")
                start = res[0].get("start", "unknown")
                return chrom, start, consequence
    except Exception as e:
        print(f"Ensembl query failed for {rsid}: {e}")
    return "unknown", "unknown", "unknown"

# 5. Build Gene-to-Feature matrix
rows = []
print("\nFetching Ensembl metadata and evaluating patient genotypes...")
for rsid, target in snp_targets.items():
    # Fetch Ensembl data
    chrom, pos_b38, consequence = get_ensembl_metadata(rsid)
    time.sleep(0.1)  # Respect Ensembl API rate limits
    
    # Retrieve genotype from 23andMe file
    genotype = genotype_map.get(rsid, "--")
    
    # Count risk alleles
    risk_allele = target["Risk_Allele"]
    if genotype == "--" or genotype == "null" or len(genotype) == 0:
        risk_count = 0
        status = "Not Genotyped"
    else:
        risk_count = sum(1 for base in genotype if base == risk_allele)
        if risk_count == 2:
            status = "Homozygous Risk"
        elif risk_count == 1:
            status = "Heterozygous Risk"
        else:
            status = "Homozygous Wildtype"
            
    # Structure row
    rows.append({
        "Gene_Symbol": target["Gene"],
        "rsID": rsid,
        "Chromosome": chrom,
        "Position_GRCh38": pos_b38,
        "Consequence": consequence,
        "Genotype": genotype,
        "Risk_Allele": risk_allele,
        "Risk_Allele_Count": risk_count,
        "Genotype_Status": status,
        "Trait_Association": target["Trait"],
        "Clinical_Implications": target["Notes"]
    })

# Save Gene-to-Feature Matrix CSV using built-in csv module
output_path = os.path.join(RESULTS_DIR, "gene_feature_matrix.csv")
fieldnames = [
    "Gene_Symbol", "rsID", "Chromosome", "Position_GRCh38", "Consequence", 
    "Genotype", "Risk_Allele", "Risk_Allele_Count", "Genotype_Status", 
    "Trait_Association", "Clinical_Implications"
]
with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"\nGene-to-Feature Matrix saved to {output_path}")

# Format console output as markdown table manually
print("\n==========================================================================")
print(" GENE-TO-FEATURE MATRIX FOR 23ANDME DATASET (LORENA SANDOVAL)")
print("==========================================================================")
print("| Gene | rsID | Chr | Pos (GRCh38) | Genotype | Risk Count | Genotype Status |")
print("| :--- | :--- | :-- | :----------- | :------- | :--------- | :-------------- |")
for r in rows:
    print(f"| {r['Gene_Symbol']} | {r['rsID']} | {r['Chromosome']} | {r['Position_GRCh38']} | {r['Genotype']} | {r['Risk_Allele_Count']} | {r['Genotype_Status']} |")
print("==========================================================================")
print("\nIndividual SNP Details:")
for r in rows:
    print(f"- {r['Gene_Symbol']} ({r['rsID']}): Genotype is {r['Genotype']}. {r['Trait_Association']} -> {r['Clinical_Implications']}")
print("==========================================================================")
