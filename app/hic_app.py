#!/usr/bin/env python3

from flask import Flask, request, render_template, jsonify, redirect, url_for, session, flash, send_from_directory
import mariadb
from functools import wraps
import os
import sys
import json
from werkzeug.utils import secure_filename
import re
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from upload_script_v5 import process_h5_file
import uuid
from datetime import datetime
import subprocess
import threading
import time
from upload_script_v5 import convert_h5_to_ginteractions, parse_ginteractions, upload_to_database
import traceback



#initialize the app
app = Flask(__name__, static_folder="static")

#Initialize a lock for thread safety
# This lock is used to ensure that only one thread can access the task status dictionary at a time
task_status_lock = threading.Lock()

# initialize notes file
NOTES_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'notes.txt')

# set a secret ket for session management (this is a placeholder, change it in production)
# This should be a strong random key in production
app.secret_key = ###

# set environment path - this is needed for the hicexplorer command line tool
os.environ["PATH"] += os.pathsep + "/usr/local/Python-3.12/bin"

# Global dictionary to track task statuses
task_status = {}

#database connection details
DB_CONFIG = {
    "host": ###,
    "user": ###,
    "password": ###,
    "database": ###,
    "port": ###
}

# Login credentials
USERNAME = ###
PASSWORD = ###

#DEFINE UPLOAD FOLDER 
UPLOAD_FOLDER = os.path.join(app.root_path, "samples3") 
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

#Define plot directory, track template and gene bed file paths
PLOT_DIR = os.path.join(app.root_path, "static/plots")
TRACK_TEMPLATE = os.path.join(app.root_path, "data/tracks_template.ini")
GENE_BED = os.path.join(app.root_path, "data/genes.bed")


# Regex pattern for filename validation
FILENAME_PATTERN = (
    r"^(?P<genotype>[^_]+)_"                      # genotype
    r"(?P<celltype>[^_]+)_"                       # cell type
    r"(?P<stage>st\d+)_"                          # stage
    r"rep(?P<replicate>\d+)_"                     # replicate
    r"(?P<bin_size>\d+)bp_"                       # bin size
    r"(?P<normalization_method>[^_]+)_"           # normalization method
    r"(?P<contacts>[\d\.]+M)"                     # contact number (e.g. 10.9M)
    r"(?:_(?P<correction>KR))?"                   # optional KR correction
    r"(?:_(?P<transformation>OE))?"               # optional OE transformation
    r"\.h5$"
)


#HELPER FUNCTIONS


#function to check if the upload file is allowed based on the filename regex pattern
def is_valid_file(filename):
    """
    Check if the filename matches the required pattern, including the .h5 extension.
    """
    return re.match(FILENAME_PATTERN, filename) is not None


#function for checking if the file is allowed
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

#database connection function
def get_db_connection():
    try:
        return mariadb.connect(**DB_CONFIG)
    except mariadb.Error as e:
        print(f"Error connecting to MariaDB Platform: {e}")
        return None

def get_hic_file_from_experiment(genotype, celltype, stage, resolution, norm, kr, oe):
    #MATCH THE USER INPUT TO THE HIC FILE NAME
    #List all .h5 files in the directory
    files = os.listdir(UPLOAD_FOLDER)
    matching_files=[]

    for file in files:
        if file.endswith(".h5"):
            match = re.match(FILENAME_PATTERN, file)
            if match:
                #extract the components from the filename
                file_genotype = match.group("genotype")
                file_celltype = match.group("celltype")
                file_stage = match.group("stage")
                file_resolution = match.group("bin_size")
                file_norm = match.group("normalization_method")
                file_kr = match.group("correction")
                file_oe = match.group("transformation")

                #check if the components match the user input
                if(
                    file_genotype == genotype and
                    file_celltype == celltype and
                    file_stage == f"st{stage}" and
                    file_resolution == str(resolution) and
                    file_norm == norm and
                    (kr and file_kr == "KR" or not kr and file_kr is None) and
                    (oe and file_oe == "OE" or not oe and file_oe is None)
                ):
                    matching_files.append(file)

    if not matching_files:
        raise ValueError(f"No matching Hi-C file found for the given parameters: {genotype}, {celltype}, st{stage}, {resolution}bp, {norm}, KR ={kr}, OE= {oe}")

    if len(matching_files) > 1:
        print(f"WARNING: Multiple matching files found. Using the first one: {matching_files[0]}")

    #Return the first matching file
    return os.path.join(UPLOAD_FOLDER, matching_files[0])

def get_region(input_str):
    if ":" in input_str:
        return input_str  # already formatted as chrom:start-end
    else:
        return lookup_gene_location(input_str)

def lookup_gene_location(gene_name):
    with open(GENE_BED, "r") as bed:
        for line in bed:
            chrom, start, end, name = line.strip().split()[:4]
            if name.lower() == gene_name.lower():
                start, end = int(start), int(end)
                midpoint = (start + end) // 2
                start_region = max(0, midpoint - 500000)
                end_region = midpoint + 500000
                return f"{chrom}:{start_region}-{end_region}"
    raise ValueError(f"Gene '{gene_name}' not found in annotation file.")


def create_ini_file(hic_file):
    new_ini = os.path.join(PLOT_DIR, "temp_tracks.ini")
    print(f"Creating .ini file with Hi-C file: {hic_file}")
    prev_section = ""
    with open(TRACK_TEMPLATE, "r") as fin, open(new_ini, "w") as fout:
        for line in fin:
            print(f"Processing line: {line.strip()}")  # Debug log
            if "file =" in line.strip() and "hic_matrix" in prev_section:
                full_path = os.path.abspath(hic_file)
                print(f"Replacing file line with: file = {full_path}")  # Debug log
                fout.write(f"file = {full_path}\n")
            else:
                fout.write(line)
            if line.strip().startswith("[") and "hic_matrix" in line.lower():
                prev_section = "hic_matrix"
                print(f"Entered hic_matrix section: {line.strip()}")  # Debug log
            elif line.strip().startswith("[") and "genes" in line.lower():
                prev_section = "genes"
                print(f"Entered genes section: {line.strip()}")  # Debug log
    print(f".ini file created at: {new_ini}")
    return new_ini


def run_pygenometracks(region, ini_path):
    output_path = os.path.join(PLOT_DIR, "output_plot.png")
    command = [
        "pyGenomeTracks",
        "--tracks", ini_path,
        "--region", region,
        "--out", output_path,
        "--width", "30",
        "--height", "15"
    ]

    result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")  # Set encoding to utf-8

    # NEW: Print stderr and stdout to logs
    if result.returncode != 0:
        print(f"pyGenomeTracks command failed! Region: {region}, INI Path: {ini_path}")
        print("Command:", " ".join(command))
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)  # ← This will contain the real error!
        raise RuntimeError("pyGenomeTracks failed.")

    return output_path

def get_experiments():
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM Experiment")
        experiments = cursor.fetchall()
        cursor.close()
        conn.close()
        return experiments
    except mariadb.Error as e:
        print(f"Error querying database: {e}")
        return []

#get the experiment ID from the database
def get_next_experiment_id():
    """
    Query the database to get the highest experiment ID and return the next ID.
    """
    try:
        conn = get_db_connection()
        if conn is None:
            raise RuntimeError("Database connection failed")

        cursor = conn.cursor()
        cursor.execute("SELECT MAX(expID) FROM Experiment")
        result = cursor.fetchone()
        cursor.close()
        conn.close()

        # If there are no records, start with ID 1
        if result[0] is None:
            return 1
        return result[0] + 1
    except mariadb.Error as e:
        print(f"Error querying database for next experiment ID: {e}")
        raise


def delete_experiment_task(experiment_id):
    try:
        conn = get_db_connection()
        if conn is None:
            print("Database connection failed.")
            return

        cursor = conn.cursor()

        # Delete related interactions first
        cursor.execute("DELETE FROM Interaction WHERE expID = %s", (experiment_id,))

        # Delete the experiment
        cursor.execute("DELETE FROM Experiment WHERE expID = %s", (experiment_id,))

        conn.commit()
        cursor.close()
        conn.close()

        print(f"Experiment {experiment_id} and its interactions successfully deleted.")
    except mariadb.Error as e:
        print(f"Error deleting experiment {experiment_id}: {e}")

# base url
@app.route('/<tab>')
def tab_view(tab):
    return render_template("kohwi_webpage.html", current_tab=tab, error=None)

#define the home and welcome routes 
@app.route("/")
@app.route("/welcome")
def index():
    return render_template("kohwi_webpage.html", current_tab="welcome", error = None)

# hi-c base tab
@app.route("/hic", methods=["POST"])
def hic_query():
    data = request.get_json()

    visualize = data.get("visualize", False)

    # Extract query parameters

    genotype1 = data["genotype1"]
    genotype2 = data["genotype2"]
    celltype1 = data["cellType1"]
    celltype2 = data["cellType2"]
    stage1 = data["stages1"]
    stage2 = data["stages2"]
    genename = data["geneName"]
    chrom = data["chrom"]
    start = data["start"]
    end = data["end"]
    no_norm = data["no_norm"]
    cooltools_norm = data["cooltools_norm"]
    hic_norm = data["hic_norm"]
    kr_correction = data["kr_correction"]
    oe_correction = data["oe_correction"]
    resolution = data["resolution"]

    # Type checks
    
    stage1 = int(stage1)
    resolution = int(resolution)

    # Map to DB values
    genotype_map = {
        "wildType": "WT",
        "UAS_Dan_FL": "UASDan",
        "UAS_Dan_DeltaPSQ": "UASDanDeltaPSQ",
        "deltaLARKS": "UASDanDeltaLARKS"
    }
    celltype_map = {
        "neuroblast": "NB",
        "wholeEmbryo": "whole"
    }

    genotype1 = genotype_map.get(genotype1, genotype1)
    celltype1 = celltype_map.get(celltype1, celltype1)
    genotype2 = genotype_map.get(genotype2, genotype2)
    celltype2 = celltype_map.get(celltype2, celltype2)

    #TEAM PLS CHECK MY LOGIC HERE 

    if (no_norm):
        norm = "adjusted"
    elif (cooltools_norm):
        norm = "mcool"
    elif (hic_norm):
        norm = "hic"

    if (kr_correction):
        kr = int(1)
    else:
        kr = int(0)

    if (oe_correction):
        oe = int(1)
    else:
        oe = int(0)

    #do not allow the user to input a gene name if they are using positional arguments
    if (chrom and start and end) and genename:
        return jsonify({"error": "Please provide either a gene name or positional arguments, both are not allowed."}), 400

    #visualization logic 
    if visualize:
        try:
            #Get the Hi-C file based on the user input
            hic_file = get_hic_file_from_experiment(genotype1, celltype1, stage1, resolution, norm, kr, oe)

            print(f"DEBUG: Hi-C file path: {hic_file}")  # Debug log

            #Get the region for visualization
            region = get_region(genename if genename else f"{chrom}:{start}-{end}")

            #create ini file
            ini_path = create_ini_file(hic_file)

            #run pyGenomeTracks
            output_plot = run_pygenometracks(region, ini_path)
            print(f"DEBUG: Output plot path: {output_plot}")  # Debug log

            #return the path to the generated plot 
            return jsonify({"success": True, "plot_path": url_for('static', filename=f'plots/{os.path.basename(output_plot)}')}), 200
        
        except ValueError as e:
            return jsonify({"success": False, "error": str(e)}), 400
        except RuntimeError as e:
            return jsonify({"success": False, "error": str(e)}), 500
        except Exception as e:
            return jsonify({"success": False, "error": f"Unexpected error: {str(e)}"}), 500


    # Query type 1: One experiment query, no gene name, no positional arguments
    if(genotype2 == "N/A" or celltype2 == "N/A" or stage2 == "N/A") and not (chrom and start and end) and not genename:
        # Query to get the interaction data
        query1 = f"""
            SELECT g1.name AS gene1_name, 
            b1.chrom AS bin1_chrom, 
            b1.start AS bin1_start, 
            b1.end AS bin1_end,                 
            g2.name AS gene2_name,
            b2.chrom as bin2_chrom , 
            b2.start as bin2_start, 
            b2.end as bin2_end, 
            I.frequency as frequency,
            E.genotype as genotype,
            E.stage as stage,
            E.cell_type as cell_type,
            E.bin_size as bin_size
            FROM Interaction I 
            JOIN Experiment E ON I.expID = E.expID
            JOIN Bin b1 ON I.bin1 = b1.binID
            JOIN Bin b2 ON I.bin2 = b2.binID
            JOIN Map m1 on m1.binID = b1.binID
            JOIN Gene g1 on m1.gID = g1.gid
            JOIN Map m2 on m2.binID = b2.binID
            JOIN Gene g2 on m2.gID = g2.gid
            WHERE E.expID IN (select E1.expID
                            from Experiment E1
                            WHERE E1.genotype =  %s
                            AND E1.stage = %s
                            AND E1.cell_type = %s
                            AND E1.bin_size = %s
                            AND E1.KR = %s
                            AND E1.OE = %s
                            AND E1.normalization_method = %s)
            LIMIT 5000;
            """
        params = [genotype1, stage1, celltype1, resolution, kr, oe, norm]

    #Query type 2: User has output experimment 1 and gene name 
    if (genotype2 == "N/A" or celltype2 == "N/A" or stage2 == "N/A") and not (chrom and start and end) and genename:
        query1 = f"""
            SELECT g1.name AS gene1_name, 
            b1.chrom AS bin1_chrom, 
            b1.start AS bin1_start, 
            b1.end AS bin1_end, 
            g2.name AS gene2_name,
            b2.chrom as bin2_chrom , 
            b2.start as bin2_start, 
            b2.end as bin2_end, 
            I.frequency as frequency,
            E.genotype as genotype,
            E.stage as stage,
            E.cell_type as cell_type,
            E.bin_size as bin_size
            FROM Interaction I 
            JOIN Experiment E ON I.expID = E.expID
            JOIN Bin b1 ON I.bin1 = b1.binID
            JOIN Bin b2 ON I.bin2 = b2.binID
            JOIN Map m1 on m1.binID = b1.binID
            JOIN Gene g1 on m1.gID = g1.gid
            JOIN Map m2 on m2.binID = b2.binID
            JOIN Gene g2 on m2.gID = g2.gid
            WHERE E.expID IN (select E1.expID
                            from Experiment E1
                            WHERE E1.genotype =  %s
                            AND E1.stage = %s
                            AND E1.cell_type = %s
                            AND E1.bin_size = %s
                            AND E1.KR = %s
                            AND E1.OE = %s
                            AND E1.normalization_method = %s)
            AND g1.name = %s
            LIMIT 100;
            """
        params = [genotype1, stage1, celltype1, resolution, kr, oe, norm, genename]

    #Query type 3: User has input one experiment and positional arguments, no gene name 
    if (genotype2 == "N/A" or celltype2 == "N/A" or stage2 == "N/A") and (chrom and start and end) and not genename:
        query1 = f"""
            SELECT g1.name AS gene1_name, 
            b1.chrom AS bin1_chrom, 
            b1.start AS bin1_start, 
            b1.end AS bin1_end, 
            g2.name AS gene2_name,
            b2.chrom as bin2_chrom , 
            b2.start as bin2_start, 
            b2.end as bin2_end, 
            I.frequency as frequency,
            E.genotype as genotype,
            E.stage as stage,
            E.cell_type as cell_type,
            E.bin_size as bin_size
            FROM Interaction I 
            JOIN Experiment E ON I.expID = E.expID
            JOIN Bin b1 ON I.bin1 = b1.binID
            JOIN Bin b2 ON I.bin2 = b2.binID
            JOIN Map m1 on m1.binID = b1.binID
            JOIN Gene g1 on m1.gID = g1.gid
            JOIN Map m2 on m2.binID = b2.binID
            JOIN Gene g2 on m2.gID = g2.gid
            WHERE E.expID IN (select E1.expID
                            from Experiment E1
                            WHERE E1.genotype =  %s
                            AND E1.stage = %s
                            AND E1.cell_type = %s
                            AND E1.bin_size = %s
                            AND E1.KR = %s
                            AND E1.OE = %s
                            AND E1.normalization_method = %s)
            AND b1.binID NOT IN (SELECT b3.binID
                            FROM Bin b3
                            WHERE b3.chrom = %s
                            AND b3.start >= %s 
                            AND b3.end <= %s)
            LIMIT 5;
        """
        params = [genotype1, stage1, celltype1, resolution, kr, oe, norm, chrom, int(end), int(start)]
    
    #Query type 4: User has input 2 experiments, no gene name, no positional arguments
    if (genotype2 != "N/A" and celltype2 != "N/A" and stage2 != "N/A") and not (chrom and start and end) and not genename:
        query1 = f"""
            SELECT g1.name AS gene1_name, 
            b1.chrom AS bin1_chrom, 
            b1.start AS bin1_start, 
            b1.end AS bin1_end, 
            g2.name AS gene2_name,
            b2.chrom as bin2_chrom , 
            b2.start as bin2_start, 
            b2.end as bin2_end, 
            I.frequency as frequency,
            E.genotype as genotype,
            E.stage as stage,
            E.cell_type as cell_type,
            E.bin_size as bin_size
            FROM Interaction I 
            JOIN Experiment E ON I.expID = E.expID
            JOIN Bin b1 ON I.bin1 = b1.binID
            JOIN Bin b2 ON I.bin2 = b2.binID
            JOIN Map m1 on m1.binID = b1.binID
            JOIN Gene g1 on m1.gID = g1.gid
            JOIN Map m2 on m2.binID = b2.binID
            JOIN Gene g2 on m2.gID = g2.gid
            WHERE E.expID IN (select E1.expID
                            from Experiment E1
                            WHERE E1.genotype =  %s 
                            OR E1.genotype = %s
                            AND E1.stage = %s
                            OR E1.stage = %s
                            AND E1.cell_type = %s
                            OR E1.cell_type = %s
                            AND E1.bin_size = %s
                            AND E1.KR = %s
                            AND E1.OE = %s
                            AND E1.normalization_method = %s)
                LIMIT 50000;
        """
        params = [genotype1, genotype2, stage1, stage2, celltype1, celltype2, resolution, kr, oe, norm]
    
    #Query type 5: User has input 2 experiments, gene name, no positional arguments
    if (genotype2 != "N/A" and celltype2 != "N/A" and stage2 != "N/A") and not (chrom and start and end) and genename:
        query1 = f"""
            SELECT g1.name AS gene1_name, 
            b1.chrom AS bin1_chrom, 
            b1.start AS bin1_start, 
            b1.end AS bin1_end, 
            g2.name AS gene2_name,
            b2.chrom as bin2_chrom , 
            b2.start as bin2_start, 
            b2.end as bin2_end, 
            I.frequency as frequency,
            E.genotype as genotype,
            E.stage as stage,
            E.cell_type as cell_type,
            E.bin_size as bin_size
            FROM Interaction I 
            JOIN Experiment E ON I.expID = E.expID
            JOIN Bin b1 ON I.bin1 = b1.binID
            JOIN Bin b2 ON I.bin2 = b2.binID
            JOIN Map m1 on m1.binID = b1.binID
            JOIN Gene g1 on m1.gID = g1.gid
            JOIN Map m2 on m2.binID = b2.binID
            JOIN Gene g2 on m2.gID = g2.gid
            WHERE E.expID IN (select E1.expID
                            from Experiment E1
                            WHERE E1.genotype =  %s 
                            OR E1.genotype = %s
                            AND E1.stage = %s
                            OR E1.stage = %s
                            AND E1.cell_type = %s
                            OR E1.cell_type = %s
                            AND E1.bin_size = %s
                            AND E1.KR = %s
                            AND E1.OE = %s
                            AND E1.normalization_method = %s)
                AND g1.name = %s
                LIMIT 5000;
        """
        params = [genotype1, genotype2, stage1, stage2, celltype1, celltype2, resolution, kr, oe, norm, genename]

    #Query type 6: User has input 2 experiments, positional arguments, no gene name

    if (genotype2 != "N/A" and celltype2 != "N/A" and stage2 != "N/A") and (chrom and start and end) and not genename:
        query1 = f"""
            SELECT g1.name AS gene1_name, 
            b1.chrom AS bin1_chrom, 
            b1.start AS bin1_start, 
            b1.end AS bin1_end, 
            g2.name AS gene2_name,
            b2.chrom as bin2_chrom , 
            b2.start as bin2_start, 
            b2.end as bin2_end, 
            I.frequency as frequency,
            E.genotype as genotype,
            E.stage as stage,
            E.cell_type as cell_type,
            E.bin_size as bin_size
            FROM Interaction I 
            JOIN Experiment E ON I.expID = E.expID
            JOIN Bin b1 ON I.bin1 = b1.binID
            JOIN Bin b2 ON I.bin2 = b2.binID
            JOIN Map m1 on m1.binID = b1.binID
            JOIN Gene g1 on m1.gID = g1.gid
            JOIN Map m2 on m2.binID = b2.binID
            JOIN Gene g2 on m2.GID = g2.gid
            WHERE E.expID IN (select E1.expID
                            from Experiment E1
                            WHERE E1.genotype =  %s 
                            OR E1.genotype = %s
                            AND E1.stage = %s
                            OR E1.stage = %s
                            AND E1.cell_type = %s
                            OR E1.cell_type = %s
                            AND E1.bin_size = %s
                            AND E1.KR = %s
                            AND E1.OE = %s
                            AND E1.normalization_method = %s)
                AND b1.binID NOT IN (SELECT b3.binID
                            FROM Bin b3
                            WHERE b3.chrom = %s
                            AND b3.start >= %s 
                            AND b3.end <= %s)
                LIMIT 5;
        """
        params = [genotype1, genotype2, stage1, stage2, celltype1, celltype2, resolution, kr, oe, norm, chrom, int(end), int(start)]

    try:
        conn = get_db_connection()
        if conn is None:
            return jsonify({"error": "Database connection failed"}), 500

        cursor = conn.cursor()
        cursor.execute(query1, params)
        results = cursor.fetchall()
        cursor.close()
        conn.close()

        if results:
            return (json.dumps(results))
        else:
            return jsonify({"no results returned"})  # Return empty list to prevent JS crash

    except mariadb.Error as e:
            return jsonify({"error": f"MariaDB error: {str(e)}"}), 501

    return jsonify({"No queries ran, input errors"})  # Also fallback to empty list if no query matched


@app.route('/genome', methods=["POST"])
def genome_query():
    data = request.get_json()
    browser_gene_name = data["browser_gene_name"]
    
    try:
        if not data:
            return jsonify({"error": "Missing gene name"}), 400
    
        

        query = """
           SELECT g.chrom, g.start, g.end 
           FROM Gene g
           WHERE g.name = %s
        """
		# Connect to the database and execute the query
        connection = get_db_connection()
        cursor = connection.cursor()
        cursor.execute(query, [browser_gene_name])
        result = cursor.fetchone()
        cursor.close()
        connection.close()

        if result:
            return(json.dumps(result))
        else:
            return jsonify({"error": "Gene not found"}), 404

    except Exception as e:
        import traceback
        traceback.print_exc()  # This will print the full error to the terminal
        return jsonify({"error: ": str(e)}), 500


# Help Tab
@app.route('/help')
def help():
    return render_template('kohwi_webpage.html', current_tab='help', error=None)


# Login Tab
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if 'logged_in' not in session:
            if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
                return jsonify({"error": "Authentication required"}), 401
            return redirect(url_for('login'))
        return f(*args, **kwargs)
    return decorated_function

# Create the login url
@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        if request.form['username'] == USERNAME and request.form['password'] == PASSWORD:
            session['logged_in'] = True
            # Here we redirect to the upload tab
            return redirect(url_for('upload'))
        else:
            error = 'Invalid credentials. Please try again.'
    return render_template('kohwi_webpage.html', current_tab='login', error=error)

# Logout redirect
@app.route('/logout')
def logout():
    session.pop('logged_in', None)
    return redirect(url_for('login'))

@app.route('/upload', methods=['GET', 'POST'])
@login_required
def upload():
    if request.method == 'POST':
        # Handle file uploads
        if 'file' in request.files:
            file = request.files['file']
            if file.filename == '':
                flash('No selected file')
                return redirect(request.url)
            if file and is_valid_file(file.filename):
                filename = file.filename
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(filepath)

                # Generate a unique task ID
                task_id = str(uuid.uuid4())
                with task_status_lock:
                    task_status[task_id] = {"status": "Processing", "progress": 0}

                # Start the background task
                thread = threading.Thread(target=process_file_task, args=(filepath, task_id))
                thread.start()

                flash('File successfully uploaded and is being processed. This will take a while depending on the size of your file. Please check back later.')
                return redirect(url_for('upload'))

        # Handle lab note submissions
        if 'note' in request.form:
            note = request.form.get('note')
            if note:
                timestamp = datetime.now().strftime("%m/%d %H:%M")
                entry = f"{timestamp}: {note}\n"
                with open(NOTES_FILE, 'a') as f:
                    f.write(entry)
                flash('Note successfully saved.')
            else:
                flash('Note cannot be empty.')

    # Fetch experiments from the database
    experiments = get_experiments()

    # Load lab notes for display
    notes = []
    notes_error = None
    try:
        if not os.path.exists(NOTES_FILE):
            open(NOTES_FILE, 'w').close()
        with open(NOTES_FILE, 'r', encoding='utf-8') as f:
            notes = f.readlines()
            notes = notes[::-1]  # Show the latest notes first
    except Exception as e:
        notes_error = f'Error loading notes: {str(e)}'

    return render_template('kohwi_webpage.html', current_tab='upload', notes=notes, error=notes_error, experiments=experiments)


def process_file_task(filepath, task_id):
    try:
        with task_status_lock:
            if task_id not in task_status:
                raise KeyError(f"Task ID {task_id} not found in task_status dictionary.")
            task_status[task_id]["status"] = "Converting file"
            task_status[task_id]["progress"] = 10
            print(f"DEBUG: Task status after setting to 'Converting file': {task_status}")

        # Step 1: Convert .h5 to .ginteractions
        ginteractions_file = convert_h5_to_ginteractions(filepath)

       # Step 2: Generate experiment ID
        experiment_id = get_next_experiment_id()
        print(f"DEBUG: Generated experiment ID: {experiment_id}")

        # Step 3: Parse .ginteractions into tables
        experiment_table, interact_table = parse_ginteractions(ginteractions_file, experiment_id)

        # Step 4: Upload tables to database
        upload_to_database(experiment_table, interact_table)

    except Exception as e:
        with task_status_lock:
            if task_id in task_status:
                task_status[task_id]["status"] = "Failed"
                task_status[task_id]["error"] = str(e)
            print(f"DEBUG: Task status after failure: {task_status}")
        # Log the full exception traceback
        print(f"Error processing file with task_id {task_id}: {str(e)}")
        traceback.print_exc()  # Log the full traceback for debugging


@app.route('/delete_experiment/<int:experiment_id>', methods=['POST'])
@login_required
def delete_experiment(experiment_id):
    try:
        # Start the deletion task in a background thread
        thread = threading.Thread(target=delete_experiment_task, args=(experiment_id,))
        thread.start()

        flash(f"Deletion of experiment {experiment_id} has started. This may take a while, depending on the size of the experiment. Come back later!", "info")
    except Exception as e:
        print(f"Error starting deletion task: {e}")
        flash(f"Error starting deletion task: {e}", "error")

    return redirect(url_for('upload'))


if __name__ == "__main__":
    app.run(debug=True)
