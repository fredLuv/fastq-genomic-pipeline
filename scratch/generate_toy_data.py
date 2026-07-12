import os
import random

# Set seed for reproducibility
random.seed(42)

# Directory setup
WORKDIR = "/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline"
REF_DIR = os.path.join(WORKDIR, "ref")
READS_DIR = os.path.join(WORKDIR, "reads")

os.makedirs(REF_DIR, exist_ok=True)
os.makedirs(READS_DIR, exist_ok=True)

# 1. Generate standard reference genome (50,000 bp E. coli mimic)
genome_length = 50000
bases = ['A', 'C', 'G', 'T']
ref_sequence = "".join(random.choices(bases, k=genome_length))

ref_fasta_path = os.path.join(REF_DIR, "ref.fa")
with open(ref_fasta_path, "w") as f:
    f.write(">ecoli_toy_ref\n")
    # Write 80 characters per line
    for i in range(0, len(ref_sequence), 80):
        f.write(ref_sequence[i:i+80] + "\n")

print(f"Generated reference genome FASTA at {ref_fasta_path}")

# 2. Introduce mutations to create the "sample" genome
# We will introduce a few distinct mutations:
# Mutation 1: SNP at coord 10,000
# Mutation 2: SNP at coord 25,000
# Mutation 3: 5bp Deletion at coord 35,000
mut_sequence = list(ref_sequence)

def mutate_base(seq, pos, new_base):
    original = seq[pos]
    seq[pos] = new_base
    print(f"Introduced SNP: pos={pos+1}, ref={original}, alt={new_base}")

# Introduced mutations (1-indexed for biological reference)
mutate_base(mut_sequence, 9999, 'C' if ref_sequence[9999] != 'C' else 'G')
mutate_base(mut_sequence, 24999, 'T' if ref_sequence[24999] != 'T' else 'A')

# Deletion of 5bp: replace 34999 to 35004 with empty string
deleted_seq = "".join(ref_sequence[34999:35004])
for idx in range(34999, 35004):
    mut_sequence[idx] = ""
print(f"Introduced Deletion: pos=35000-35004, sequence={deleted_seq}")

mutated_genome = "".join(mut_sequence)

# 3. Generate paired-end reads from the mutated genome
# Coverage: 40x
# Read length: 150bp
# Fragment size: mean=300bp, std=30bp
# Number of reads: 7000
num_reads = 7000
read_len = 150
fragment_mean = 300
fragment_std = 30

# Illumina TruSeq Adapter Sequence
adapter_r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

fastq_r1_path = os.path.join(READS_DIR, "read1.fastq")
fastq_r2_path = os.path.join(READS_DIR, "read2.fastq")

with open(fastq_r1_path, "w") as f1, open(fastq_r2_path, "w") as f2:
    for r_idx in range(num_reads):
        # Pick random fragment start position
        frag_len = int(random.normalvariate(fragment_mean, fragment_std))
        # Bounds check
        frag_len = max(read_len + 10, min(frag_len, 500))
        
        start_pos = random.randint(0, len(mutated_genome) - frag_len)
        end_pos = start_pos + frag_len
        
        fragment = mutated_genome[start_pos:end_pos]
        
        # Read 1 is prefix of fragment
        r1_seq = fragment[:read_len]
        # Read 2 is reverse complement of fragment suffix
        r2_raw = fragment[-read_len:]
        
        # Reverse complement function
        comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        r2_seq = "".join(comp.get(b, 'N') for b in reversed(r2_raw))
        
        # Add adapters to 10% of reads to let fastp show its trimming power
        if random.random() < 0.10:
            # Clip sequences and append adapter
            clip_pos = random.randint(100, 140)
            r1_seq = r1_seq[:clip_pos] + adapter_r1[:150-clip_pos]
            r2_seq = r2_seq[:clip_pos] + adapter_r2[:150-clip_pos]
            
        # Introduce random sequencing errors (0.5% rate) with low quality score
        def apply_errors_and_quality(seq):
            seq_list = list(seq)
            qual_list = []
            for i in range(len(seq_list)):
                if random.random() < 0.005:
                    seq_list[i] = random.choice([b for b in bases if b != seq_list[i]])
                    # Phred 10 (ASCII: '+' which is 10 + 33 = 43)
                    qual_list.append('+')
                else:
                    # Phred 35-40 (ASCII: 'H' to 'M')
                    qual_list.append(chr(random.randint(68, 73) + 33))
            return "".join(seq_list), "".join(qual_list)
        
        r1_seq, r1_qual = apply_errors_and_quality(r1_seq)
        r2_seq, r2_qual = apply_errors_and_quality(r2_seq)
        
        # Write FASTQ records
        header = f"@SRR1234567.{r_idx+1} 1/1"
        f1.write(f"{header}\n{r1_seq}\n+\n{r1_qual}\n")
        
        header2 = f"@SRR1234567.{r_idx+1} 2/2"
        f2.write(f"{header2}\n{r2_seq}\n+\n{r2_qual}\n")

print(f"Generated simulated reads in {READS_DIR}")
