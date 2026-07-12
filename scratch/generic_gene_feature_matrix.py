import os
import csv
import json
import time
import urllib.request

# Configuration
WORKDIR = "/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline"
READS_DIR = os.path.join(WORKDIR, "reads")
RESULTS_DIR = os.path.join(WORKDIR, "results")
os.makedirs(READS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

raw_file_path = os.path.join(READS_DIR, "23andme.txt")
cache_file_path = os.path.join(RESULTS_DIR, "ensembl_metadata_cache.json")
snp_output_path = os.path.join(RESULTS_DIR, "generic_snp_features.csv")
gene_output_path = os.path.join(RESULTS_DIR, "generic_gene_features.csv")

# 1. Broad Set of Targeted Genetic Markers (27 SNPs across 21 Genes)
snp_database = {
    # Methylation & Brain Health
    "rs1801133": {"Gene": "MTHFR", "Category": "Methylation", "Trait": "Folate metabolism (C677T)", "Risk_Allele": "T", "Notes": "TT: ~16% folate conversion; CT: ~65%; CC: Typical"},
    "rs1801131": {"Gene": "MTHFR", "Category": "Methylation", "Trait": "Folate metabolism (A1298C)", "Risk_Allele": "C", "Notes": "CC: Reduced folate conversion; AC: Mildly reduced"},
    "rs4680": {"Gene": "COMT", "Category": "Methylation/Brain", "Trait": "Dopamine clearance (Val158Met)", "Risk_Allele": "A", "Notes": "AA (Met/Met): Slow dopamine breakdown ('Worrier'); GG (Val/Val): Fast breakdown ('Warrior'); AG: Intermediate"},
    "rs6265": {"Gene": "BDNF", "Category": "Brain Health", "Trait": "Neuroplasticity (Val66Met)", "Risk_Allele": "A", "Notes": "AA/AG: Met allele associated with lower BDNF secretion and memory/stress resilience variation"},
    "rs1800497": {"Gene": "DRD2", "Category": "Brain Health", "Trait": "Dopamine D2 Receptor (Taq1A)", "Risk_Allele": "A", "Notes": "AA/AG: Reduced receptor density; associated with reward-seeking behavior and addiction susceptibility"},
    "rs6318": {"Gene": "MAOA", "Category": "Brain Health", "Trait": "Monoamine oxidase A activity", "Risk_Allele": "G", "Notes": "G allele associated with altered serotonin/dopamine degradation rates"},

    # Cardiovascular, Lipids & Longevity
    "rs429358": {"Gene": "APOE", "Category": "Cardiovascular/Lipids", "Trait": "Alzheimer's & lipid risk (C112R)", "Risk_Allele": "C", "Notes": "C allele determines APOE-epsilon4 (Increased Alzheimer's & lipid risk)"},
    "rs7412": {"Gene": "APOE", "Category": "Cardiovascular/Lipids", "Trait": "Alzheimer's & lipid risk (R158C)", "Risk_Allele": "T", "Notes": "T allele determines APOE-epsilon2 (Reduced Alzheimer's risk)"},
    "rs1799983": {"Gene": "NOS3", "Category": "Cardiovascular", "Trait": "Nitric Oxide Synthase (vasodilation)", "Risk_Allele": "T", "Notes": "TT/GT: Reduced nitric oxide production, associated with increased hypertension risk"},
    "rs2802292": {"Gene": "FOXO3", "Category": "Longevity", "Trait": "Healthy aging & longevity", "Risk_Allele": "C", "Notes": "CC/CG: Associated with increased lifespan, protection against cardiovascular disease"},

    # Diet, Lactose & Metabolism
    "rs4988235": {"Gene": "LCT", "Category": "Dietary Tolerance", "Trait": "Lactase persistence (Lactose tolerance)", "Risk_Allele": "G", "Notes": "GG: Lactose intolerant; AA/AG: Lactose tolerant"},
    "rs9939609": {"Gene": "FTO", "Category": "Metabolism", "Trait": "Weight and obesity risk", "Risk_Allele": "A", "Notes": "AA: Increased obesity risk and fat mass; TA: Intermediate; TT: Typical risk"},
    "rs1801282": {"Gene": "PPARG", "Category": "Metabolism", "Trait": "Insulin sensitivity & fat storage", "Risk_Allele": "G", "Notes": "GG/CG: Typical sensitivity; CC (Pro/Pro): Associated with higher BMI but lower diabetes risk"},
    "rs1801260": {"Gene": "CLOCK", "Category": "Circadian Rhythm", "Trait": "Sleep duration & night owl tendency", "Risk_Allele": "G", "Notes": "GG/AG: Associated with evening chronotype ('night owl') and shorter sleep duration"},

    # Vitamins & Iron
    "rs1800562": {"Gene": "HFE", "Category": "Iron Metabolism", "Trait": "Hemochromatosis (C282Y)", "Risk_Allele": "A", "Notes": "AA: High risk of iron overload; GA: Carrier; GG: Typical"},
    "rs1799945": {"Gene": "HFE", "Category": "Iron Metabolism", "Trait": "Hemochromatosis (H63D)", "Risk_Allele": "G", "Notes": "GG: Moderate risk of iron overload; CG: Carrier; CC: Typical"},
    "rs731236": {"Gene": "VDR", "Category": "Vitamin Absorption", "Trait": "Vitamin D Receptor (TaqI)", "Risk_Allele": "G", "Notes": "GG (tt): Associated with reduced vitamin D absorption/bone density; AA: Normal"},
    "rs7501331": {"Gene": "BCMO1", "Category": "Vitamin Absorption", "Trait": "Beta-carotene to Vitamin A conversion", "Risk_Allele": "C", "Notes": "CC/TC: Up to 57% reduction in conversion of plant-based beta-carotene to active Vitamin A"},
    "rs602662": {"Gene": "FUT2", "Category": "Vitamin Absorption", "Trait": "Vitamin B12 absorption", "Risk_Allele": "A", "Notes": "AA: Associated with lower B12 absorption and plasma levels; GG: Normal"},
    "rs4654748": {"Gene": "NBPF3", "Category": "Vitamin Absorption", "Trait": "Vitamin B6 clearance rate", "Risk_Allele": "T", "Notes": "TT/CT: Associated with faster B6 clearance and lower plasma levels"},

    # Detoxification & Inflammation
    "rs1695": {"Gene": "GSTP1", "Category": "Detoxification", "Trait": "Glutathione S-transferase (detox speed)", "Risk_Allele": "G", "Notes": "GG: Reduced detoxification capacity for toxins/heavy metals; AG: Intermediate; AA: Normal"},
    "rs4880": {"Gene": "SOD2", "Category": "Antioxidant", "Trait": "Superoxide dismutase (antioxidant defense)", "Risk_Allele": "T", "Notes": "TT: Reduced transport of SOD2 to mitochondria, associated with higher oxidative stress; CC: Normal"},
    "rs671": {"Gene": "ALDH2", "Category": "Detoxification", "Trait": "Alcohol flush reaction (Glu504Lys)", "Risk_Allele": "A", "Notes": "AA: Severe flush/hangover; GA: Mild flush; GG: Normal tolerance"},
    "rs1800795": {"Gene": "IL6", "Category": "Inflammation", "Trait": "Interleukin-6 systemic inflammation", "Risk_Allele": "G", "Notes": "GG/CG: Associated with higher baseline IL-6 levels and systemic inflammation"},
    "rs1800629": {"Gene": "TNF", "Category": "Inflammation", "Trait": "Tumor Necrosis Factor (inflammatory response)", "Risk_Allele": "A", "Notes": "AA/GA: Higher baseline TNF-alpha, associated with increased inflammatory/autoimmune risk"},

    # Traits & Physical Characteristics
    "rs1815739": {"Gene": "ACTN3", "Category": "Physical Traits", "Trait": "Muscle performance (R577X)", "Risk_Allele": "T", "Notes": "CC: Sprinters/Power athletes; CT: Intermediate; TT: Endurance athletes"},
    "rs1805007": {"Gene": "MC1R", "Category": "Physical Traits", "Trait": "Red hair & fair skin", "Risk_Allele": "T", "Notes": "TT: High probability of red hair; CT: Carrier; CC: Typical"}
}

# 2. Load or Initialize Ensembl API Metadata Cache
metadata_cache = {}
if os.path.exists(cache_file_path):
    print(f"Loading Ensembl metadata cache from {cache_file_path}...")
    try:
        with open(cache_file_path, 'r') as f:
            metadata_cache = json.load(f)
        print(f"Loaded {len(metadata_cache)} cached SNPs.")
    except Exception as e:
        print(f"Failed to load cache: {e}")

# 3. Read 23andMe Genotyping File
print(f"Loading raw 23andMe genotypes from {raw_file_path}...")
genotype_map = {}
if os.path.exists(raw_file_path):
    with open(raw_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                genotype_map[parts[0]] = parts[3]
    print(f"Indexed {len(genotype_map)} genotypes.")
else:
    print(f"ERROR: Raw file {raw_file_path} not found. Please run the script inside your workspace.")
    exit(1)

# 4. Batch Resolve Missing SNPs from Ensembl VEP POST Endpoint
missing_snps = [rsid for rsid in snp_database if rsid not in metadata_cache]
if missing_snps:
    print(f"Querying Ensembl VEP API in batch for {len(missing_snps)} missing SNPs...")
    url = "https://rest.ensembl.org/vep/human/id"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0"
    }
    # Ensembl supports up to 1000 IDs per request
    chunk_size = 1000
    for i in range(0, len(missing_snps), chunk_size):
        chunk = missing_snps[i:i+chunk_size]
        body = json.dumps({"ids": chunk}).encode("utf-8")
        req = urllib.request.Request(url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req) as r:
                results = json.loads(r.read().decode("utf-8"))
                for item in results:
                    rsid = item.get("id")
                    consequence = item.get("most_severe_consequence", "unknown")
                    chrom = item.get("seq_region_name", "unknown")
                    start = item.get("start", "unknown")
                    
                    # Extract official gene symbol from transcripts
                    gene_symbol = "unknown"
                    tc = item.get("transcript_consequences", [])
                    if tc:
                        # Prefer transcripts that contain gene_symbol
                        for t in tc:
                            if "gene_symbol" in t:
                                gene_symbol = t["gene_symbol"]
                                break
                    
                    metadata_cache[rsid] = {
                        "Chromosome": chrom,
                        "Position": start,
                        "Consequence": consequence,
                        "Official_Gene": gene_symbol
                    }
            print(f"Successfully resolved chunk of {len(chunk)} SNPs.")
            time.sleep(0.1)  # Throttle limit safeguard
        except Exception as e:
            print(f"Ensembl batch query failed: {e}")
            
    # Save the updated cache back to disk
    with open(cache_file_path, 'w') as f:
        json.dump(metadata_cache, f, indent=2)
    print("Metadata cache updated.")

# 5. Process Individual SNPs and Compile Features
processed_snps = []
for rsid, target in snp_database.items():
    meta = metadata_cache.get(rsid, {"Chromosome": "unknown", "Position": "unknown", "Consequence": "unknown", "Official_Gene": target["Gene"]})
    genotype = genotype_map.get(rsid, "--")
    
    # Calculate Risk Allele Counts
    risk_allele = target["Risk_Allele"]
    if genotype == "--" or genotype == "null" or not genotype:
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
            
    processed_snps.append({
        "rsID": rsid,
        "Gene_Symbol": target["Gene"],
        "Category": target["Category"],
        "Chromosome": meta["Chromosome"],
        "Position_GRCh38": meta["Position"],
        "Consequence": meta["Consequence"],
        "Genotype": genotype,
        "Risk_Allele": risk_allele,
        "Risk_Allele_Count": risk_count,
        "Genotype_Status": status,
        "Trait_Association": target["Trait"],
        "Clinical_Implications": target["Notes"]
    })

# Save Individual SNP details
fieldnames_snp = [
    "rsID", "Gene_Symbol", "Category", "Chromosome", "Position_GRCh38", "Consequence",
    "Genotype", "Risk_Allele", "Risk_Allele_Count", "Genotype_Status", "Trait_Association", "Clinical_Implications"
]
with open(snp_output_path, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames_snp)
    writer.writeheader()
    for s in processed_snps:
        writer.writerow(s)
print(f"SNP-level features saved to {snp_output_path}")

# 6. Aggregate by Gene to Build the Gene-to-Feature Matrix
gene_aggregates = {}
for s in processed_snps:
    gene = s["Gene_Symbol"]
    if gene not in gene_aggregates:
        gene_aggregates[gene] = {
            "Gene_Symbol": gene,
            "Category": s["Category"],
            "Chromosome": s["Chromosome"],
            "Total_SNPs_Measured": 0,
            "Genotyped_SNPs": 0,
            "Risk_Alleles_Count": 0,
            "SNP_Genotypes": [],
            "Consequences": set()
        }
    
    ga = gene_aggregates[gene]
    ga["Total_SNPs_Measured"] += 1
    if s["Genotype"] != "--":
        ga["Genotyped_SNPs"] += 1
        ga["Risk_Alleles_Count"] += s["Risk_Allele_Count"]
        ga["SNP_Genotypes"].append(f"{s['rsID']}:{s['Genotype']}")
    if s["Consequence"] != "unknown":
        ga["Consequences"].add(s["Consequence"])

# Finalize Gene features and calculate Average Risk Score (genetic load)
gene_matrix_rows = []
for gene, ga in gene_aggregates.items():
    total_alleles = ga["Genotyped_SNPs"] * 2
    avg_risk_score = ga["Risk_Alleles_Count"] / total_alleles if total_alleles > 0 else 0.0
    
    # Format list fields as semicolon separated strings
    snp_string = "; ".join(ga["SNP_Genotypes"]) if ga["SNP_Genotypes"] else "None"
    consequence_string = "; ".join(ga["Consequences"]) if ga["Consequences"] else "unknown"
    
    gene_matrix_rows.append({
        "Gene_Symbol": ga["Gene_Symbol"],
        "Category": ga["Category"],
        "Chromosome": ga["Chromosome"],
        "Total_SNPs_Measured": ga["Total_SNPs_Measured"],
        "Genotyped_SNPs": ga["Genotyped_SNPs"],
        "Total_Risk_Alleles": ga["Risk_Alleles_Count"],
        "Average_Risk_Score": round(avg_risk_score, 3),
        "SNP_Genotypes": snp_string,
        "Most_Severe_Consequences": consequence_string
    })

# Save Gene-to-Feature Matrix
fieldnames_gene = [
    "Gene_Symbol", "Category", "Chromosome", "Total_SNPs_Measured", 
    "Genotyped_SNPs", "Total_Risk_Alleles", "Average_Risk_Score", 
    "SNP_Genotypes", "Most_Severe_Consequences"
]
with open(gene_output_path, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames_gene)
    writer.writeheader()
    for g in gene_matrix_rows:
        writer.writerow(g)
print(f"Gene-level feature matrix saved to {gene_output_path}")

# 7. Print aggregated markdown Gene-to-Feature Matrix
print("\n==========================================================================================")
print(" GENE-TO-FEATURE AGGREGATED MATRIX FOR LORENA SANDOVAL")
print("==========================================================================================")
print("| Gene | Category | Chromosome | SNPs | Genotyped | Risk Alleles | Risk Score | SNP Genotypes |")
print("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
for g in gene_matrix_rows:
    print(f"| {g['Gene_Symbol']} | {g['Category']} | {g['Chromosome']} | {g['Total_SNPs_Measured']} | {g['Genotyped_SNPs']} | {g['Total_Risk_Alleles']} | {g['Average_Risk_Score']} | {g['SNP_Genotypes']} |")
print("==========================================================================================")
