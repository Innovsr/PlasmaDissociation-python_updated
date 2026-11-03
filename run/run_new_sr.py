#########################################################################################
#
# Python run script for PlasmaDissociation
#
# Dynamic Slurm restart/validation system for AIRES
#
# Workflow:
#
#   Setup
#      |
#      | afterok
#      v
#   Plasma stage 1
#      |
#      | afterany
#      v
#   restart_validation.py
#      |
#      | if incomplete
#      v
#   Plasma stage 2
#      |
#      | afterany
#      v
#   restart_validation.py
#      |
#      | if incomplete
#      v
#   Plasma stage 3
#      |
#      ...
#
# The validator decides whether another restart is required.
#
#########################################################################################

import sys
import socket
import os
import subprocess
import getpass
import shutil
import json
from pathlib import Path


#########################################################################################
# BASIC SETUP
#########################################################################################

if __name__ != "__main__":
    raise SystemExit


# Read inputs from the directory from which run_new_sr.py is executed
with open("inputs.json", "r") as f:
    inputs = json.load(f)


# Use exactly the Python interpreter running this script
PYTHON_EXE = sys.executable

# Directory containing run_new_sr.py
SCRIPT_DIR = Path(__file__).resolve().parent


#########################################################################################
# CHECK INPUTS
#########################################################################################

if not isinstance(inputs["HPC"]["repeats"], int):
    sys.exit("Number of repeats must be an integer")

if not isinstance(inputs["HPC"]["cores"], int):
    sys.exit("Number of parallel cores must be an integer")

if inputs["HPC"]["repeats"] < 1:
    sys.exit("Number of repeats must be at least 1")

if inputs["HPC"]["cores"] < 1:
    sys.exit("Number of cores must be at least 1")

if inputs["HPC"]["cores"] > 16:
    sys.exit("Maximum of 16 cores available")


print("Arguments checked")


#########################################################################################
# DETECT HPC
#########################################################################################

hostname = socket.gethostname()

HPCFLG = 0

if "arc4" in hostname or "arc3" in hostname:
    HPCFLG = "arc"

elif "aire" in hostname:
    HPCFLG = "aire"

print("Hostname:", hostname)
print("HPCFLG:", HPCFLG)


#########################################################################################
# DETERMINE EXECUTION DIRECTORY
#########################################################################################

if HPCFLG == 0:

    EXDIR = "../EXEC"

    if not os.path.exists(EXDIR):
        os.mkdir(EXDIR)

elif HPCFLG == "arc":

    EXDIR = "/nobackup/" + getpass.getuser()

    if not os.path.exists(EXDIR):
        os.mkdir(EXDIR)

elif HPCFLG == "aire":

    EXDIR = "/mnt/scratch/" + getpass.getuser()

    if not os.path.exists(EXDIR):
        os.mkdir(EXDIR)

else:
    sys.exit("Unknown HPC system")


#########################################################################################
# RUN DIRECTORY
#########################################################################################

runfolder_path = os.path.join(
    EXDIR,
    inputs["HPC"]["Runfolder"]
)


if os.path.exists(runfolder_path):

    value = input(
        "Runfolder already exists. Do you want to delete it? y/n\n"
    ).strip().lower()

    if value == "y":
        shutil.rmtree(runfolder_path)

    else:
        sys.exit(
            "Runfolder already exists. "
            "Change Runfolder name or delete/move the existing run."
        )


os.mkdir(runfolder_path)

EXDIR1 = os.path.abspath(runfolder_path)

print("Run directory:", EXDIR1)


#########################################################################################
# CREATE DIRECTORY STRUCTURE
#########################################################################################

os.mkdir(os.path.join(EXDIR1, "results"))
os.mkdir(os.path.join(EXDIR1, "repetitions"))
os.mkdir(os.path.join(EXDIR1, "results", "bonds"))
os.mkdir(os.path.join(EXDIR1, "results", "specifics"))
os.mkdir(os.path.join(EXDIR1, "results", "graphs"))
os.mkdir(os.path.join(EXDIR1, "setup"))
os.mkdir(os.path.join(EXDIR1, "setup", "tmp"))
os.mkdir(os.path.join(EXDIR1, "logs"))


with open(
    os.path.join(EXDIR1, "results", "bonds", "allbonds.out"),
    "w"
):
    pass


#########################################################################################
# INITIAL GEOMETRY
#########################################################################################

if inputs["Molecule_data"]["Geom_flg"] == "Initial":

    source = os.path.join(
        SCRIPT_DIR,
        "..",
        "Initial_Geometry",
        inputs["Molecule_data"]["Molecule"] + ".txt"
    )

    destination = os.path.join(
        EXDIR1,
        "setup",
        "initial_guess.txt"
    )

    shutil.copy2(source, destination)


#########################################################################################
# COPY IMPORTANT FILES INTO RUN DIRECTORY
#########################################################################################

# Copy inputs.json
shutil.copy2(
    os.path.join(SCRIPT_DIR, "inputs.json"),
    os.path.join(EXDIR1, "inputs.json")
)


# Copy restart.py if it exists
restart_py_source = os.path.join(SCRIPT_DIR, "restart.py")

if os.path.isfile(restart_py_source):
    shutil.copy2(
        restart_py_source,
        os.path.join(EXDIR1, "restart.py")
    )


# Copy restart_validation.py
validation_source = os.path.join(
    SCRIPT_DIR,
    "restart_validation.py"
)

if not os.path.isfile(validation_source):
    sys.exit(
        "ERROR: restart_validation.py was not found in "
        + str(SCRIPT_DIR)
    )

shutil.copy2(
    validation_source,
    os.path.join(EXDIR1, "restart_validation.py")
)


#########################################################################################
# CREATE REPETITION DIRECTORIES
#########################################################################################

for i in range(inputs["HPC"]["repeats"]):

    repdir = os.path.join(
        EXDIR1,
        "repetitions",
        "rep-" + str(i + 1)
    )

    os.mkdir(repdir)
    os.mkdir(os.path.join(repdir, "output"))
    os.mkdir(os.path.join(repdir, "tmp"))
    os.mkdir(os.path.join(repdir, "checks"))


    if inputs["Molecule_data"]["Geom_flg"] == "Full":

        source = os.path.join(
            SCRIPT_DIR,
            "..",
            inputs["Molecule_data"]["Molecule"],
            "Geom",
            "Geometry."
            + str(
                i + inputs["Molecule_data"]["Geom_file_start"]
            )
        )

        destination = os.path.join(
            repdir,
            "Geometry"
        )

        shutil.copy2(source, destination)


#########################################################################################
# LOCAL RUN: COPY CODE
#########################################################################################

if HPCFLG == 0:

    shutil.copytree(
        os.path.join(SCRIPT_DIR, "..", "code"),
        os.path.join(EXDIR1, "code")
    )


#########################################################################################
# FULL GEOMETRY: COPY BOND FILE
#########################################################################################

if inputs["Molecule_data"]["Geom_flg"] == "Full":

    source = os.path.join(
        SCRIPT_DIR,
        "..",
        inputs["Molecule_data"]["Molecule"],
        "bondarr.txt"
    )

    destination = os.path.join(
        EXDIR1,
        "results"
    )

    shutil.copy2(source, destination)


#########################################################################################
# CHANGE TO RUN DIRECTORY
#########################################################################################

os.chdir(EXDIR1)

EXDIR1 = os.getcwd()

print("Working directory:", EXDIR1)


#########################################################################################
# ARC
#########################################################################################

if HPCFLG == "arc":

    HOME_DIR = str(SCRIPT_DIR)

    #############################################################################
    # SETUP
    #############################################################################

    if inputs["Molecule_data"]["Geom_flg"] in [
        "PubChem",
        "Initial"
    ]:

        setup_script = (
            "Setup_"
            + inputs["HPC"]["Runfolder"]
            + ".sh"
        )

        with open(setup_script, "w") as f:

            f.write("#$ -cwd -V\n")
            f.write("#$ -l h_vmem=4G,h_rt=06:00:00\n")

            f.write(
                "#$ -N Setup_"
                + inputs["HPC"]["Runfolder"]
                + "\n"
            )

            f.write(
                "#$ -pe smp "
                + str(inputs["HPC"]["cores"])
                + "\n"
            )

            if inputs["Elec_structure"]["method"] == "QChem":

                f.write("module load test qchem\n")
                f.write("mkdir -p $TMPDIR/qchemlocal\n")

                f.write(
                    "tar -xzvf /nobackup/"
                    + getpass.getuser()
                    + "/qchem.tar.gz "
                    "-C $TMPDIR/qchemlocal\n"
                )

                f.write(
                    "qchemlocal=$TMPDIR/qchemlocal\n"
                )

                f.write(
                    'export QCHEM_HOME="$qchemlocal"\n'
                )

                f.write(
                    'export QC="$qchemlocal"\n'
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
                    "export QCSCRATCH="
                    + EXDIR1
                    + "/setup/tmp\n"
                )

            f.write(
                "cd "
                + EXDIR1
                + "/setup\n"
            )

            f.write(
                PYTHON_EXE
                + " "
                + HOME_DIR
                + "/../code/setup.py\n"
            )

        subprocess.call(
            [
                "qsub",
                "-N",
                "Setup_" + inputs["HPC"]["Runfolder"],
                setup_script
            ]
        )


    #############################################################################
    # INITIAL PLASMA
    #############################################################################

    plasma_script = (
        "Plasma_"
        + inputs["HPC"]["Runfolder"]
        + "_1.sh"
    )

    with open(plasma_script, "w") as f:

        f.write("#$ -cwd -V\n")

        f.write(
            "#$ -l h_vmem=8G,h_rt="
            + str(inputs["HPC"]["hours"])
            + ":00:00\n"
        )

        f.write(
            "#$ -N Plasma_"
            + inputs["HPC"]["Runfolder"]
            + "_1\n"
        )

        f.write(
            "#$ -pe smp "
            + str(inputs["HPC"]["cores"])
            + "\n"
        )

        f.write(
            "#$ -t 1-"
            + str(inputs["HPC"]["repeats"])
            + "\n"
        )

        f.write(
            "cd "
            + EXDIR1
            + "/repetitions/rep-$SGE_TASK_ID\n"
        )

        f.write(
            PYTHON_EXE
            + " "
            + HOME_DIR
            + "/../code/main.py\n"
        )


    if inputs["Molecule_data"]["Geom_flg"] in [
        "PubChem",
        "Initial"
    ]:

        subprocess.call(
            [
                "qsub",
                "-N",
                "Plasma_" + inputs["HPC"]["Runfolder"] + "_1",
                "-hold_jid",
                "Setup_" + inputs["HPC"]["Runfolder"],
                plasma_script
            ]
        )

    else:

        subprocess.call(
            [
                "qsub",
                "-N",
                "Plasma_" + inputs["HPC"]["Runfolder"] + "_1",
                plasma_script
            ]
        )


#########################################################################################
# AIRES
#########################################################################################

elif HPCFLG == "aire":

    #############################################################################
    # PATHS
    #############################################################################

    # Project code is discovered from the location of run_new_sr.py.
    PROJECT_DIR = SCRIPT_DIR.parent
    CODE_DIR = PROJECT_DIR / "code"

    SETUP_PY = CODE_DIR / "setup.py"
    MAIN_PY = CODE_DIR / "main.py"

    if not SETUP_PY.is_file():
        sys.exit(
            "Cannot find setup.py at "
            + str(SETUP_PY)
        )

    if not MAIN_PY.is_file():
        sys.exit(
            "Cannot find main.py at "
            + str(MAIN_PY)
        )


    #############################################################################
    # SETUP JOB
    #############################################################################

    setup_job_id = None

    if inputs["Molecule_data"]["Geom_flg"] in [
        "PubChem",
        "Initial"
    ]:

        setup_script = (
            "Setup_"
            + inputs["HPC"]["Runfolder"]
            + ".sh"
        )

        with open(setup_script, "w") as f:

            f.write("#!/bin/bash\n")

            f.write(
                "#SBATCH --job-name=Setup_"
                + inputs["HPC"]["Runfolder"]
                + "\n"
            )

            f.write("#SBATCH --time=02:00:00\n")

            f.write(
                "#SBATCH --mem="
                + str(inputs["HPC"]["memory"])
                + "G\n"
            )

            f.write(
                "#SBATCH --cpus-per-task="
                + str(inputs["HPC"]["cores"])
                + "\n"
            )

            f.write(
                "#SBATCH --output="
                + EXDIR1
                + "/logs/setup_%j.out\n"
            )

            f.write(
                "#SBATCH --error="
                + EXDIR1
                + "/logs/setup_%j.err\n"
            )

            f.write("set -euo pipefail\n")

            f.write(
                'echo "HOSTNAME=$(hostname)"\n'
            )

            f.write(
                'echo "PYTHON='
                + PYTHON_EXE
                + '"\n'
            )

            f.write(
                'echo "WORKDIR=$(pwd)"\n'
            )

            f.write(
                'export PLASMA_INPUTS="'
                + EXDIR1
                + '/inputs.json"\n'
            )

            f.write(
                'export TMPDIR="${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}"\n'
            )

            f.write(
                'mkdir -p "$TMPDIR"\n'
            )

            if inputs["Elec_structure"]["method"] == "QChem":

                f.write(
                    "module unload qchem "
                    "2>/dev/null || true\n"
                )

                f.write(
                    "module load qchem/6.0.2\n"
                )

                f.write(
                    "export LM_LICENSE_FILE="
                    "27000@uol-lnx-lic01.leeds.ac.uk\n"
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
                    'echo "QCHEM_HOME=$QCHEM_HOME"\n'
                )

                f.write(
                    'ls "$QCHEM_HOME/exe/qcprog.exe"\n'
                )

                f.write(
                    'export QCSCRATCH="'
                    + EXDIR1
                    + '/setup/tmp"\n'
                )

            f.write(
                'cd "'
                + EXDIR1
                + '/setup"\n'
            )

            f.write(
                'which python || true\n'
            )

            f.write(
                '"'
                + PYTHON_EXE
                + '" "'
                + str(SETUP_PY)
                + '"\n'
            )


        setup_job_id = subprocess.check_output(
            [
                "sbatch",
                "--parsable",
                setup_script
            ],
            text=True
        ).strip().split(";")[0]

        print(
            "Setup job submitted:",
            setup_job_id
        )


    #############################################################################
    # FUNCTION: WRITE PLASMA SCRIPT
    #############################################################################

    def write_aire_plasma_script(
        script_name,
        stage_number,
        array_spec
    ):

        with open(script_name, "w") as f:

            f.write("#!/bin/bash\n")

            f.write(
                "#SBATCH --job-name=Plasma_"
                + inputs["HPC"]["Runfolder"]
                + "_"
                + str(stage_number)
                + "\n"
            )

            f.write(
                "#SBATCH --mem="
                + str(inputs["HPC"]["memory"])
                + "G\n"
            )

            f.write(
                "#SBATCH --time="
                + str(inputs["HPC"]["hours"])
                + ":00:00\n"
            )

            f.write(
                "#SBATCH --cpus-per-task="
                + str(inputs["HPC"]["cores"])
                + "\n"
            )

            f.write(
                "#SBATCH --array="
                + str(array_spec)
                + "\n"
            )

            f.write(
                "#SBATCH --output="
                + EXDIR1
                + "/logs/plasma_"
                + str(stage_number)
                + "_%A_%a.out\n"
            )

            f.write(
                "#SBATCH --error="
                + EXDIR1
                + "/logs/plasma_"
                + str(stage_number)
                + "_%A_%a.err\n"
            )

            f.write(
                "set -euo pipefail\n"
            )

            f.write(
                'echo "HOSTNAME=$(hostname)"\n'
            )

            f.write(
                'echo "PYTHON='
                + PYTHON_EXE
                + '"\n'
            )

            f.write(
                'echo "TASK_ID=$SLURM_ARRAY_TASK_ID"\n'
            )

            f.write(
                'echo "PLASMA_STAGE='
                + str(stage_number)
                + '"\n'
            )

            f.write(
                'export PLASMA_INPUTS="'
                + EXDIR1
                + '/inputs.json"\n'
            )

            f.write(
                'export TMPDIR="${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}"\n'
            )

            f.write(
                'mkdir -p "$TMPDIR"\n'
            )

            f.write(
                'cd "'
                + EXDIR1
                + '/repetitions/rep-$SLURM_ARRAY_TASK_ID"\n'
            )


            #####################################################################
            # IMPORTANT:
            #
            # DO NOT check .bomd_complete.
            #
            # The validator has already inspected the trajectory files and
            # decided that this repetition is incomplete.
            #####################################################################


            if inputs["Elec_structure"]["method"] == "QChem":

                f.write(
                    "module unload qchem "
                    "2>/dev/null || true\n"
                )

                f.write(
                    "module load qchem/6.0.2\n"
                )

                f.write(
                    "export LM_LICENSE_FILE="
                    "27000@uol-lnx-lic01.leeds.ac.uk\n"
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
                    'echo "QCHEM_HOME=$QCHEM_HOME"\n'
                )

                f.write(
                    'ls "$QCHEM_HOME/exe/qcprog.exe"\n'
                )

                f.write(
                    'export QCSCRATCH="'
                    + EXDIR1
                    + '/repetitions/rep-$SLURM_ARRAY_TASK_ID/tmp"\n'
                )


            f.write(
                "which python || true\n"
            )

            f.write(
                '"'
                + PYTHON_EXE
                + '" "'
                + str(MAIN_PY)
                + '"\n'
            )

            f.write(
                'echo "Plasma stage '
                + str(stage_number)
                + ' finished."\n'
            )


        os.chmod(script_name, 0o755)


    #############################################################################
    # INITIAL PLASMA ARRAY
    #############################################################################

    initial_array = (
        "1-"
        + str(inputs["HPC"]["repeats"])
    )

    plasma_script = (
        "Plasma_"
        + inputs["HPC"]["Runfolder"]
        + "_1.sh"
    )

    write_aire_plasma_script(
        plasma_script,
        1,
        initial_array
    )


    if setup_job_id is not None:

        plasma_command = [
            "sbatch",
            "--parsable",
            "--dependency=afterok:" + setup_job_id,
            plasma_script
        ]

    else:

        plasma_command = [
            "sbatch",
            "--parsable",
            plasma_script
        ]


    plasma_job_id = subprocess.check_output(
        plasma_command,
        text=True
    ).strip().split(";")[0]


    print(
        "Initial Plasma job:",
        plasma_job_id
    )


    #############################################################################
    # VALIDATION WRAPPER
    #
    # This wrapper is used for EVERY validation cycle.
    #
    # It is deliberately a Slurm script rather than submitting the Python
    # file directly.
    #############################################################################

    validation_wrapper = os.path.join(
        EXDIR1,
        "restart_validation_job.sh"
    )

    with open(validation_wrapper, "w") as f:

        f.write("#!/bin/bash\n")

        f.write(
            "#SBATCH --job-name=Validate_"
            + inputs["HPC"]["Runfolder"]
            + "\n"
        )

        f.write(
            "#SBATCH --time=00:30:00\n"
        )

        f.write(
            "#SBATCH --mem=2G\n"
        )

        f.write(
            "#SBATCH --cpus-per-task=1\n"
        )

        f.write(
            "#SBATCH --output="
            + EXDIR1
            + "/logs/restart_validation_%j.out\n"
        )

        f.write(
            "#SBATCH --error="
            + EXDIR1
            + "/logs/restart_validation_%j.err\n"
        )

        f.write(
            "set -euo pipefail\n"
        )

        f.write(
            'export PLASMA_INPUTS="'
            + EXDIR1
            + '/inputs.json"\n'
        )

        f.write(
            'cd "'
            + EXDIR1
            + '"\n'
        )

        f.write(
            'echo "PREVIOUS_PLASMA_JOB_ID=$PREVIOUS_PLASMA_JOB_ID"\n'
        )

        f.write(
            '"'
            + PYTHON_EXE
            + '" "'
            + EXDIR1
            + '/restart_validation.py"\n'
        )


    os.chmod(
        validation_wrapper,
        0o755
    )


    #############################################################################
    # SUBMIT FIRST VALIDATOR
    #
    # IMPORTANT:
    #
    # afterany means validator runs whether Plasma:
    #
    #   COMPLETED
    #   TIMEOUT
    #   FAILED
    #   CANCELLED
    #
    #############################################################################

    validation_job_id = subprocess.check_output(
        [
            "sbatch",
            "--parsable",
            "--export=ALL,PREVIOUS_PLASMA_JOB_ID=" + plasma_job_id,
            "--dependency=afterany:" + plasma_job_id,
            validation_wrapper
        ],
        text=True
    ).strip().split(";")[0]


    print(
        "Initial validation job:",
        validation_job_id
    )


    print()
    print("============================================================")
    print("Dynamic Plasma restart system submitted.")
    print("============================================================")
    print("Setup job:       ", setup_job_id)
    print("Initial Plasma:  ", plasma_job_id)
    print("Validator:       ", validation_job_id)
    print("============================================================")


#########################################################################################
# LOCAL EXECUTION
#########################################################################################

else:

    for i in range(
        inputs["HPC"]["repeats"]
    ):

        repdir = os.path.join(
            EXDIR1,
            "repetitions",
            "rep-" + str(i + 1)
        )

        os.chdir(repdir)

        if inputs["Molecule_data"]["Geom_flg"] in [
            "PubChem",
            "Initial"
        ]:

            subprocess.call(
                [
                    PYTHON_EXE,
                    "../../../code/setup.py"
                ]
            )

        subprocess.call(
            [
                PYTHON_EXE,
                "../../../code/main.py"
            ]
        )
