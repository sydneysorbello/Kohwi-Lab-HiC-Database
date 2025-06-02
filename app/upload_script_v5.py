import subprocess
import os
import re
import pandas as pd
import mariadb
import getpass

# --- CONFIG --- 
HICEXPLORER_CMD = "/usr/local/Python-3.12/bin/hicConvertFormat"  # assume hicConvertFormat is installed
OUTPUT_DIR = "/var/www/html/students_25/Team3/app2/samples3"
os.makedirs(OUTPUT_DIR, exist_ok=True)

DB_USER = ###
DB_PASS = ###
DB_HOST = ###
DB_NAME = ###
DB_PORT = ###

# --- FUNCTIONS ---

def convert_h5_to_ginteractions(h5_file, output_dir="/var/www/html/students_25/Team3/app2/samples3"):
    """Convert .h5 file to .ginteractions using hicConvertFormat."""
    base_name = os.path.basename(h5_file).replace('.h5', '')
    output_file = os.path.join(output_dir, f"{base_name}.ginteractions")
    
    print(f"Converting {h5_file} to {output_file}...")
    os.environ["PATH"] += os.pathsep + "/usr/local/Python-3.12/bin/"    
    subprocess.run([
        HICEXPLORER_CMD,
        "-m", h5_file,
        "-o", output_file,
        "--inputFormat", "h5",
        "--outputFormat", "ginteractions"
    ], check=True)

    output_file_with_tsv = f"{output_file}.tsv"
    return output_file_with_tsv

def parse_ginteractions(filepath, experiment_id):
    """Parse ginteractions file into Experiment and Interaction tables."""
    bins = {}
    interactions = []
    bin_id = 0

    filename = os.path.basename(filepath).strip()
    print(f"Parsing file: {filename}")
    # same regex as before
    pattern = (
    r"^(?P<genotype>[^_]+)_"                  # [genotype]
    r"(?P<celltype>[^_]+)_"                   # [celltype]
    r"st(?P<stage>\d+)_"                      # st[stage]
    r"rep(?P<replicate>\d+)_"                 # rep[replicate]
    r"(?P<bin_size>\d+)bp_"                   # [bin_size]bp
    r"(?P<normalization_method>[^_]+)_"       # [normalization_method]
    r"(?P<contacts>[\d\.]+)M"                 # [contacts]M
    r"(?:_(?P<correction>KR))?"               # optional _KR
    r"(?:_(?P<transformation>OE))?"           # optional _OE
    r"\.ginteractions(\.tsv)?$"                  # .ginteractions.tsv
)
    print("Filename before match:", repr(filename))  # Check for hidden chars
    match = re.match(pattern, filename.strip())
    if not match:
        raise ValueError(f"HEY Filename format incorrect!!!: {filename}")

    data = match.groupdict()

    stage = int(re.search(r'\d+', data['stage']).group())
    replicate = int(data['replicate'])
    bin_size = int(data['bin_size'])
    hiC_contacts = float(data['contacts'].replace('M', '')) * 1_000_000
    normalization_method = data['normalization_method']
    KR = int(data['correction'] == 'KR')
    OE = int(data['transformation'] == 'OE')

    # Create DataFrame for experiment table
    print(f"Creating experiment table for Experiment {experiment_id}...")
    experiment_table = pd.DataFrame([{
        'expID': experiment_id,
        'genotype': data['genotype'],
        'cell_type': data['celltype'],
        'stage': stage,
        'replicate': replicate,
        'bin_size': bin_size,
        'HiC_contacts': int(hiC_contacts),
        'normalization_method': normalization_method,
        'KR': KR,
        'OE': OE
    }])

    # Read the ginteractions file
    with open(filepath, 'r') as f:
        for line in f:
            chrom1, start1, end1, chrom2, start2, end2, freq = line.strip().split('\t')
            start1, end1, start2, end2 = map(int, [start1, end1, start2, end2])
            freq = float(freq)
            chrom1 = chrom1.replace('chr', '')
            chrom2 = chrom2.replace('chr', '')

            bin1 = (chrom1, start1, end1)
            bin2 = (chrom2, start2, end2)

            if bin1 not in bins:
                bins[bin1] = bin_id
                bin_id += 1
            if bin2 not in bins:
                bins[bin2] = bin_id
                bin_id += 1

            interactions.append({
                'bin1': bins[bin1],
                'bin2': bins[bin2],
                'frequency': freq,
                'expID': experiment_id
            })
    # Create DataFrames for interactions
    print(f"Creating interactions table...")
    interact_table = pd.DataFrame(interactions)
    return experiment_table, interact_table

def upload_to_database(experiment_table, interact_table):
    """Upload parsed tables to MariaDB."""
    print("Connecting to database...")
    conn = mariadb.connect(
        user=DB_USER,
        port= DB_PORT,
        password=DB_PASS,
        host=DB_HOST,
        database=DB_NAME,
        local_infile=True
    )
    cur = conn.cursor()

    print("Uploading experiment table to database...")
    # upload experiment table
    for _, row in experiment_table.iterrows():
        cur.execute("""
            INSERT INTO Experiment 
            (expID, genotype, cell_type, stage, replicate, bin_size, HiC_contacts, normalization_method, KR, OE) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, tuple(row))

    print("Uploading interaction table to database this may take a while...")
    for _, row in interact_table.iterrows():
        cur.execute("""
            INSERT INTO Interaction
            (bin1, bin2, frequency, expID)
            VALUES (?, ?, ?, ?)
            """, tuple(row))

    conn.commit()
    cur.close()
    conn.close()

# --- MAIN ---

def process_h5_file(h5_file, experiment_id):
    """Process an uploaded .h5 file: convert, parse, and upload."""
    # 1. Convert .h5 → .ginteractions
    ginteractions_file = convert_h5_to_ginteractions(h5_file)
    
    # 2. Parse .ginteractions → experiment and interaction tables
    experiment_table, interact_table = parse_ginteractions(ginteractions_file, experiment_id)
    
    # 3. Upload parsed tables into database
    upload_to_database(experiment_table, interact_table)

