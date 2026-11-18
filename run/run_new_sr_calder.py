#!/usr/bin/env python3

################################################################################
#
# run_new_sr_calder.py
#
# CALDER-ONLY launcher for PlasmaDissociation-python_updated
#
# Workflow:
#
#       SETUP
#         |
#         | afterok
#         v
#   BOMD ARRAY - STAGE 1
#         |
#         | afterany
#         v
#     VALIDATOR
#         |
#         | incomplete
#         v
#   BOMD ARRAY - STAGE 2
#         |
#         | afterany
#         v
#     VALIDATOR
#         |
#         | incomplete
#         v
#       ...
#
# The validator decides whether another BOMD stage is required.
#
# Restart point is determined by the BOMD trajectory files:
#
#   output/xyz.all
#   output/momenta.all
#   output/forces.all
#   output/time.all
#
# The latest COMMON timestep in these four files is used as the
# restart point by code/init.py.
#
################################################################################


import sys
import socket
import os
import subprocess
import getpass
import shutil
import json
from pathlib import Path


################################################################################
# BASIC SETUP
################################################################################

if __name__ != "__main__":
    raise SystemExit


################################################################################
# PATHS
################################################################################

SCRIPT_DIR = Path(__file__).resolve().parent

# Repository root:
#
# PlasmaDissociation-python_updated/
# ├── code/
# ├── run/
# └── Initial_Geometry/
#
PROJECT_DIR = SCRIPT_DIR.parent

CODE_DIR = PROJECT_DIR / "code"
INITIAL_GEOMETRY_DIR = PROJECT_DIR / "Initial_Geometry"

SETUP_PY = CODE_DIR / "setup.py"
MAIN_PY = CODE_DIR / "main.py"


################################################################################
# INPUT FILE
################################################################################

INPUT_FILE = SCRIPT_DIR / "inputs.json"

if not INPUT_FILE.is_file():

    sys.exit(
        "ERROR: inputs.json was not found at:\n"
        + str(INPUT_FILE)
    )


with open(INPUT_FILE, "r") as f:

    inputs = json.load(f)


################################################################################
# PYTHON INTERPRETER
################################################################################

# Use exactly the Python interpreter running this launcher.
#
# On Calder this should be the plasma environment:
#
# /users/lgjs798/.conda/envs/plasma/bin/python
#
PYTHON_EXE = sys.executable


################################################################################
# CHECK PROJECT FILES
################################################################################

if not SETUP_PY.is_file():

    sys.exit(
        "ERROR: setup.py was not found at:\n"
        + str(SETUP_PY)
    )


if not MAIN_PY.is_file():

    sys.exit(
        "ERROR: main.py was not found at:\n"
        + str(MAIN_PY)
    )


################################################################################
# CHECK INPUTS
################################################################################

try:

    repeats = inputs["HPC"]["repeats"]
    cores = inputs["HPC"]["cores"]
    hours = inputs["HPC"]["hours"]
    memory = inputs["HPC"]["memory"]
    account = inputs["HPC"]["account"]
    runfolder = inputs["HPC"]["Runfolder"]

except KeyError as e:

    sys.exit(
        "ERROR: Missing HPC input:\n"
        + str(e)
    )


if not isinstance(repeats, int):

    sys.exit(
        "ERROR: HPC repeats must be an integer."
    )


if not isinstance(cores, int):

    sys.exit(
        "ERROR: HPC cores must be an integer."
    )


if repeats < 1:

    sys.exit(
        "ERROR: Number of repeats must be at least 1."
    )


if cores < 1:

    sys.exit(
        "ERROR: Number of cores must be at least 1."
    )


if hours < 1:

    sys.exit(
        "ERROR: HPC hours must be at least 1."
    )


################################################################################
# CALDER SETTINGS
################################################################################

# Calder CPU partition.
#
# Your BOMD calculation is CPU based.
#
CALDER_PARTITION = "calder-cpu-turin"


# Calder scratch.
#
# Calder documentation specifies:
#
# /scratch-calder/users/<username>
#
CALDER_SCRATCH = Path(
    "/scratch-calder/users"
) / getpass.getuser()


################################################################################
# ELECTRONIC STRUCTURE
################################################################################

try:

    METHOD = inputs["Elec_structure"]["method"]

except KeyError:

    sys.exit(
        "ERROR: inputs.json does not contain:\n"
        "\"Elec_structure\": {\"method\": ...}"
    )


if METHOD != "QChem":

    sys.exit(
        "ERROR: This Calder launcher is currently designed "
        "for QChem BOMD.\n"
        "Current method = "
        + str(METHOD)
    )


################################################################################
# INITIAL GEOMETRY SETTINGS
################################################################################

try:

    GEOM_FLAG = inputs["Molecule_data"]["Geom_flg"]

except KeyError:

    sys.exit(
        "ERROR: Molecule_data.Geom_flg is missing from inputs.json."
    )


################################################################################
# PRINT INFORMATION
################################################################################

print()
print("=" * 70)
print("CALDER BOMD RESTART WORKFLOW")
print("=" * 70)

print(
    "Hostname              :",
    socket.gethostname()
)

print(
    "Python                :",
    PYTHON_EXE
)

print(
    "Repository            :",
    PROJECT_DIR
)

print(
    "Code directory        :",
    CODE_DIR
)

print(
    "Calder partition      :",
    CALDER_PARTITION
)

print(
    "Account               :",
    account
)

print(
    "Runfolder             :",
    runfolder
)

print(
    "Repeats               :",
    repeats
)

print(
    "CPUs                  :",
    cores
)

print(
    "Memory                :",
    memory,
    "GB"
)

print(
    "Walltime              :",
    hours,
    "hours"
)

print(
    "Electronic structure  :",
    METHOD
)

print(
    "Geometry flag         :",
    GEOM_FLAG
)

print(
    "Calder scratch        :",
    CALDER_SCRATCH
)

print("=" * 70)


################################################################################
# CREATE CALDER SCRATCH
################################################################################

CALDER_SCRATCH.mkdir(
    parents=True,
    exist_ok=True
)


################################################################################
# RUN DIRECTORY
################################################################################

RUN_DIR = CALDER_SCRATCH / runfolder


################################################################################
# SAFETY CHECK
################################################################################

if RUN_DIR.exists():

    print()
    print(
        "WARNING: Runfolder already exists:"
    )

    print(
        RUN_DIR
    )

    print()

    answer = input(
        "Do you want to DELETE this existing runfolder? "
        "Type y to continue: "
    ).strip().lower()

    if answer != "y":

        sys.exit(
            "Existing runfolder preserved. "
            "No job was submitted."
        )

    print(
        "Deleting existing runfolder..."
    )

    shutil.rmtree(RUN_DIR)


################################################################################
# CREATE DIRECTORY STRUCTURE
################################################################################

RUN_DIR.mkdir(
    parents=True
)

(RUN_DIR / "results").mkdir()
(RUN_DIR / "results" / "bonds").mkdir()
(RUN_DIR / "results" / "specifics").mkdir()
(RUN_DIR / "results" / "graphs").mkdir()

(RUN_DIR / "repetitions").mkdir()

(RUN_DIR / "setup").mkdir()
(RUN_DIR / "setup" / "tmp").mkdir()

(RUN_DIR / "logs").mkdir()


################################################################################
# COPY INPUT FILE
################################################################################

shutil.copy2(
    INPUT_FILE,
    RUN_DIR / "inputs.json"
)


################################################################################
# INITIAL GEOMETRY
################################################################################

if GEOM_FLAG == "Initial":

    molecule = inputs["Molecule_data"]["Molecule"]

    source = (
        INITIAL_GEOMETRY_DIR
        / (molecule + ".txt")
    )

    destination = (
        RUN_DIR
        / "setup"
        / "initial_guess.txt"
    )

    if not source.is_file():

        sys.exit(
            "ERROR: Initial geometry not found:\n"
            + str(source)
        )

    shutil.copy2(
        source,
        destination
    )

    print(
        "Initial geometry copied:"
    )

    print(
        source
    )


################################################################################
# PUBCHEM
################################################################################

if GEOM_FLAG == "PubChem":

    print(
        "Geometry mode = PubChem"
    )

    print(
        "setup.py will obtain the required geometry."
    )


################################################################################
# CREATE REPETITION DIRECTORIES
################################################################################

for i in range(1, repeats + 1):

    repdir = (
        RUN_DIR
        / "repetitions"
        / ("rep-" + str(i))
    )

    repdir.mkdir()

    (repdir / "output").mkdir()
    (repdir / "tmp").mkdir()
    (repdir / "checks").mkdir()


################################################################################
# FULL GEOMETRY
################################################################################

if GEOM_FLAG == "Full":

    molecule = inputs["Molecule_data"]["Molecule"]

    bond_source = (
        PROJECT_DIR
        / molecule
        / "bondarr.txt"
    )

    if not bond_source.is_file():

        sys.exit(
            "ERROR: bondarr.txt not found:\n"
            + str(bond_source)
        )

    shutil.copy2(
        bond_source,
        RUN_DIR / "results"
    )

    geom_start = inputs[
        "Molecule_data"
    ]["Geom_file_start"]

    for i in range(repeats):

        geometry_number = (
            i + geom_start
        )

        source = (
            PROJECT_DIR
            / molecule
            / "Geom"
            / ("Geometry." + str(geometry_number))
        )

        destination = (
            RUN_DIR
            / "repetitions"
            / ("rep-" + str(i + 1))
            / "Geometry"
        )

        if not source.is_file():

            sys.exit(
                "ERROR: Geometry file not found:\n"
                + str(source)
            )

        shutil.copy2(
            source,
            destination
        )


################################################################################
# WRITE SETUP SLURM SCRIPT
################################################################################

SETUP_SCRIPT = (
    RUN_DIR
    / ("Setup_" + runfolder + ".sh")
)


with open(SETUP_SCRIPT, "w") as f:

    f.write("#!/bin/bash\n")

    f.write(
        "#SBATCH --job-name=Setup_"
        + runfolder
        + "\n"
    )

    f.write(
        "#SBATCH --partition="
        + CALDER_PARTITION
        + "\n"
    )

    f.write(
        "#SBATCH --account="
        + account
        + "\n"
    )

    f.write(
        "#SBATCH --nodes=1\n"
    )

    f.write(
        "#SBATCH --ntasks=1\n"
    )

    f.write(
        "#SBATCH --cpus-per-task="
        + str(cores)
        + "\n"
    )

    f.write(
        "#SBATCH --mem="
        + str(memory)
        + "G\n"
    )

    f.write(
        "#SBATCH --time=02:00:00\n"
    )

    f.write(
        "#SBATCH --output="
        + str(RUN_DIR / "logs" / "setup_%j.out")
        + "\n"
    )

    f.write(
        "#SBATCH --error="
        + str(RUN_DIR / "logs" / "setup_%j.err")
        + "\n"
    )

    f.write("\n")

    f.write(
        "set -euo pipefail\n"
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Calder module environment
    # ------------------------------------------------------------

    f.write(
        "export MODULEPATH=/uol/apps/calder/modules\n"
    )

    f.write(
        "module load calder/cpu\n"
    )

    f.write(
        "module load miniforge3/26.1.1-3\n"
    )

    f.write(
        'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
    )

    f.write(
        "conda activate plasma\n"
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------

    f.write(
        'echo "HOSTNAME=$(hostname)"\n'
    )

    f.write(
        'echo "PYTHON=$(which python)"\n'
    )

    f.write(
        'python --version\n'
    )

    f.write(
        'python -c "import pyqchem; print(\'PyQChem:\', pyqchem.__version__)"\n'
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Q-Chem
    #
    # IMPORTANT:
    # The Q-Chem module name should be confirmed on Calder with:
    #
    # module spider qchem
    #
    # We use the same module name currently used in your workflow.
    # ------------------------------------------------------------

    f.write(
        "module load qchem/6.0.2\n"
    )

    f.write("\n")

    f.write(
        'echo "QCHEM_HOME=$QCHEM_HOME"\n'
    )

    f.write(
        "which qchem\n"
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Q-Chem environment
    # ------------------------------------------------------------

    f.write(
        'export QC="$QCHEM_HOME"\n'
    )

    f.write(
        'export QCAUX="$QC/qcaux"\n'
    )

    f.write(
        'export QCPROG="$QC/exe/qcprog.exe"\n'
    )

    f.write(
        'export QCPROG_S="$QC/exe/qcprog.exe_s"\n'
    )

    f.write(
        'export PATH="$PATH:$QC/exe:$QC/bin"\n'
    )

    f.write(
        'export LM_LICENSE_FILE="27000@uol-lnx-lic01.leeds.ac.uk"\n'
    )

    f.write(
        'export QCSCRATCH="'
        + str(RUN_DIR / "setup" / "tmp")
        + '"\n'
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Run setup.py
    #
    # setup.py expects to run from:
    #
    # RUN_DIR/setup
    #
    # and reads:
    #
    # ../inputs.json
    #
    # ------------------------------------------------------------

    f.write(
        "cd "
        + str(RUN_DIR / "setup")
        + "\n"
    )

    f.write(
        "echo \"Running setup.py\"\n"
    )

    f.write(
        '"'
        + PYTHON_EXE
        + '" "'
        + str(SETUP_PY)
        + '"\n'
    )


os.chmod(
    SETUP_SCRIPT,
    0o755
)


################################################################################
# SUBMIT SETUP
################################################################################

setup_job_id = subprocess.check_output(
    [
        "sbatch",
        "--parsable",
        str(SETUP_SCRIPT)
    ],
    text=True
).strip().split(";")[0]


print()
print(
    "Setup job submitted:",
    setup_job_id
)


################################################################################
# WRITE BOMD ARRAY SCRIPT
################################################################################

PLASMA_SCRIPT = (
    RUN_DIR
    / ("Plasma_" + runfolder + "_1.sh")
)


with open(PLASMA_SCRIPT, "w") as f:

    f.write("#!/bin/bash\n")

    f.write(
        "#SBATCH --job-name=MACE_" + runfolder + "\n"
    )

    f.write(
        "#SBATCH --partition="
        + CALDER_PARTITION
        + "\n"
    )

    f.write(
        "#SBATCH --account="
        + account
        + "\n"
    )

    f.write(
        "#SBATCH --nodes=1\n"
    )

    f.write(
        "#SBATCH --ntasks=1\n"
    )

    f.write(
        "#SBATCH --cpus-per-task="
        + str(cores)
        + "\n"
    )

    f.write(
        "#SBATCH --mem="
        + str(memory)
        + "G\n"
    )

    f.write(
        "#SBATCH --time="
        + str(hours)
        + ":00:00\n"
    )

    f.write(
        "#SBATCH --array=1-"
        + str(repeats)
        + "\n"
    )

    f.write(
        "#SBATCH --output="
        + str(
            RUN_DIR
            / "logs"
            / "plasma_%A_%a.out"
        )
        + "\n"
    )

    f.write(
        "#SBATCH --error="
        + str(
            RUN_DIR
            / "logs"
            / "plasma_%A_%a.err"
        )
        + "\n"
    )

    f.write("\n")

    f.write(
        "set -euo pipefail\n"
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Calder environment
    # ------------------------------------------------------------

    f.write(
        "export MODULEPATH=/uol/apps/calder/modules\n"
    )

    f.write(
        "module load calder/cpu\n"
    )

    f.write(
        "module load miniforge3/26.1.1-3\n"
    )

    f.write(
        'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
    )

    f.write(
        "conda activate plasma\n"
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------

    f.write(
        'echo "HOSTNAME=$(hostname)"\n'
    )

    f.write(
        'echo "PYTHON=$(which python)"\n'
    )

    f.write(
        'python --version\n'
    )

    f.write(
        'python -c "import pyqchem; print(\'PyQChem:\', pyqchem.__version__)"\n'
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Q-Chem
    # ------------------------------------------------------------

    f.write(
        "module load qchem/6.0.2\n"
    )

    f.write(
        'export QC="$QCHEM_HOME"\n'
    )

    f.write(
        'export QCAUX="$QC/qcaux"\n'
    )

    f.write(
        'export QCPROG="$QC/exe/qcprog.exe"\n'
    )

    f.write(
        'export QCPROG_S="$QC/exe/qcprog.exe_s"\n'
    )

    f.write(
        'export PATH="$PATH:$QC/exe:$QC/bin"\n'
    )

    f.write(
        'export LM_LICENSE_FILE="27000@uol-lnx-lic01.leeds.ac.uk"\n'
    )

    f.write(
        'export QCSCRATCH="'
        + str(
            RUN_DIR
            / "repetitions"
            / "rep-$SLURM_ARRAY_TASK_ID"
            / "tmp"
        )
        + '"\n'
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Select repetition directory
    # ------------------------------------------------------------

    f.write(
        'REP_DIR="'
        + str(
            RUN_DIR
            / "repetitions"
        )
        + '/rep-$SLURM_ARRAY_TASK_ID"\n'
    )

    f.write(
        'cd "$REP_DIR"\n'
    )

    f.write("\n")

    f.write(
        'echo "========================================"\n'
    )

    f.write(
        'echo "BOMD ARRAY WORKER"\n'
    )

    f.write(
        'echo "JOB ID       : $SLURM_JOB_ID"\n'
    )

    f.write(
        'echo "ARRAY JOB ID : $SLURM_ARRAY_JOB_ID"\n'
    )

    f.write(
        'echo "ARRAY TASK   : $SLURM_ARRAY_TASK_ID"\n'
    )

    f.write(
        'echo "REPETITION   : $REP_DIR"\n'
    )

    f.write(
        'echo "PYTHON       : $(which python)"\n'
    )

    f.write(
        'echo "========================================"\n'
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Run main.py
    #
    # main.py will decide whether this is:
    #
    #   new trajectory
    #
    # or:
    #
    #   restart from latest common timestep
    #
    # ------------------------------------------------------------

    f.write(
        '"'
        + PYTHON_EXE
        + '" "'
        + str(MAIN_PY)
        + '"\n'
    )


os.chmod(
    PLASMA_SCRIPT,
    0o755
)


################################################################################
# SUBMIT INITIAL BOMD ARRAY
################################################################################

plasma_job_id = subprocess.check_output(
    [
        "sbatch",
        "--parsable",
        "--dependency=afterok:" + setup_job_id,
        str(PLASMA_SCRIPT)
    ],
    text=True
).strip().split(";")[0]


print(
    "Initial BOMD array submitted:",
    plasma_job_id
)


################################################################################
# WRITE VALIDATION SCRIPT
################################################################################

VALIDATION_SCRIPT = (
    RUN_DIR
    / "restart_validation_job.sh"
)


with open(VALIDATION_SCRIPT, "w") as f:

    f.write("#!/bin/bash\n")

    f.write(
        "#SBATCH --job-name=Validate_"
        + runfolder
        + "\n"
    )

    f.write(
        "#SBATCH --partition="
        + CALDER_PARTITION
        + "\n"
    )

    f.write(
        "#SBATCH --account="
        + account
        + "\n"
    )

    f.write(
        "#SBATCH --nodes=1\n"
    )

    f.write(
        "#SBATCH --ntasks=1\n"
    )

    f.write(
        "#SBATCH --cpus-per-task=1\n"
    )

    f.write(
        "#SBATCH --mem=2G\n"
    )

    f.write(
        "#SBATCH --time=00:30:00\n"
    )

    f.write(
        "#SBATCH --output="
        + str(
            RUN_DIR
            / "logs"
            / "restart_validation_%j.out"
        )
        + "\n"
    )

    f.write(
        "#SBATCH --error="
        + str(
            RUN_DIR
            / "logs"
            / "restart_validation_%j.err"
        )
        + "\n"
    )

    f.write("\n")

    f.write(
        "set -euo pipefail\n"
    )

    f.write("\n")

    # ------------------------------------------------------------
    # Calder environment
    # ------------------------------------------------------------

    f.write(
        "export MODULEPATH=/uol/apps/calder/modules\n"
    )

    f.write(
        "module load calder/cpu\n"
    )

    f.write(
        "module load miniforge3/26.1.1-3\n"
    )

    f.write(
        'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
    )

    f.write(
        "conda activate plasma\n"
    )

    f.write("\n")

    f.write(
        'echo "HOSTNAME=$(hostname)"\n'
    )

    f.write(
        'echo "PYTHON=$(which python)"\n'
    )

    f.write(
        'echo "PREVIOUS_PLASMA_JOB_ID=$PREVIOUS_PLASMA_JOB_ID"\n'
    )

    f.write("\n")

    f.write(
        "cd "
        + str(RUN_DIR)
        + "\n"
    )

    f.write("\n")

    f.write(
        '"'
        + PYTHON_EXE
        + '" "'
        + str(
            RUN_DIR
            / "restart_validation.py"
        )
        + '"\n'
    )


os.chmod(
    VALIDATION_SCRIPT,
    0o755
)


################################################################################
# SUBMIT FIRST VALIDATOR
################################################################################

validation_job_id = subprocess.check_output(
    [
        "sbatch",
        "--parsable",
        "--export=ALL,PREVIOUS_PLASMA_JOB_ID=" + plasma_job_id,
        "--dependency=afterany:" + plasma_job_id,
        str(VALIDATION_SCRIPT)
    ],
    text=True
).strip().split(";")[0]


################################################################################
# FINAL INFORMATION
################################################################################

print()
print("=" * 70)
print("CALDER WORKFLOW SUBMITTED")
print("=" * 70)

print(
    "Run directory       :",
    RUN_DIR
)

print(
    "Setup job            :",
    setup_job_id
)

print(
    "Initial BOMD array   :",
    plasma_job_id
)

print(
    "Initial validator    :",
    validation_job_id
)

print(
    "Repetitions          :",
    repeats
)

print(
    "CPUs per repetition  :",
    cores
)

print(
    "Walltime             :",
    hours,
    "hours"
)

print()
print(
    "Workflow:"
)

print(
    "  Setup"
)

print(
    "    |"
)

print(
    "    | afterok"
)

print(
    "    v"
)

print(
    "  BOMD array"
)

print(
    "    |"
)

print(
    "    | afterany"
)

print(
    "    v"
)

print(
    "  Validator"
)

print(
    "    |"
)

print(
    "    | incomplete"
)

print(
    "    v"
)

print(
    "  Restart BOMD array"
)

print(
    "    |"
)

print(
    "    | afterany"
)

print(
    "    v"
)

print(
    "  Validator ..."
)

print()
print(
    "Restart point is determined from the latest common"
)

print(
    "trajectory timestep in xyz.all, momenta.all,"
)

print(
    "forces.all and time.all."
)

print("=" * 70)
