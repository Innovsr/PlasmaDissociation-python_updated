#!/usr/bin/env python3

"""
restart_validation.py

Validation-driven Slurm restart controller for PlasmaDissociation.

The script:

1. Reads inputs.json automatically from the current run directory.
2. Checks every repetition.
3. Finds the latest timestep common to:
       output/xyz.all
       output/momenta.all
       output/forces.all
       output/time.all
4. Compares it with:
       Propagation.Tot_timesteps * Propagation.Timestep
5. Restarts ONLY incomplete repetitions.
6. Submits another validator using an afterany dependency.
7. Stops automatically when all repetitions reach the requested
   final timestep.

.bomd_complete is NOT used as the completion criterion.

The project code path is discovered automatically.
No hard-coded PlasmaDissociation-python path is required.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path


# ============================================================
# AUTOMATIC PATH DISCOVERY
# ============================================================

# Directory containing this script.
SCRIPT_DIR = Path(__file__).resolve().parent

# The run directory is where this validator is located.
RUN_DIR = SCRIPT_DIR

# inputs.json is in the run directory.
INPUTS_FILE = RUN_DIR / "inputs.json"

# Repetition directories.
REPETITIONS_DIR = RUN_DIR / "repetitions"

# Log directory.
LOG_DIR = RUN_DIR / "logs"


def find_code_directory():
    """
    Automatically find the project code directory.

    We look for a directory containing BOTH:

        code/main.py
        code/setup.py

    This means the code can be located in:

        PlasmaDissociation-python/
        PlasmaDissociation-python_updated/
        or another project directory.

    Nothing is hard-coded.
    """

    # --------------------------------------------------------
    # First search upwards from this script.
    # --------------------------------------------------------

    current = SCRIPT_DIR

    for directory in [current] + list(current.parents):

        candidate = directory / "code"

        if (
            candidate.is_dir()
            and (candidate / "main.py").is_file()
            and (candidate / "setup.py").is_file()
        ):
            return candidate

    # --------------------------------------------------------
    # If the validator is in the scratch run directory,
    # search under the user's home directory.
    # --------------------------------------------------------

    home = Path.home()

    search_roots = [
        home,
        home / "PlasmaDissociation-python",
        home / "PlasmaDissociation-python_updated",
    ]

    checked = set()

    for root in search_roots:

        if not root.exists():
            continue

        try:
            root = root.resolve()
        except OSError:
            continue

        if root in checked:
            continue

        checked.add(root)

        candidate = root / "code"

        if (
            candidate.is_dir()
            and (candidate / "main.py").is_file()
            and (candidate / "setup.py").is_file()
        ):
            return candidate

    # --------------------------------------------------------
    # Last automatic search under HOME.
    # --------------------------------------------------------

    try:

        for candidate in home.glob("**/code"):

            if (
                candidate.is_dir()
                and (candidate / "main.py").is_file()
                and (candidate / "setup.py").is_file()
            ):
                return candidate

    except Exception:
        pass

    raise FileNotFoundError(
        "Could not automatically find the project code directory "
        "containing code/main.py and code/setup.py."
    )


# Automatically determine project code location.
CODE_DIR = find_code_directory()

MAIN_PY = CODE_DIR / "main.py"
SETUP_PY = CODE_DIR / "setup.py"


# ============================================================
# TRAJECTORY FILES
# ============================================================

TRAJECTORY_FILES = (
    "xyz.all",
    "momenta.all",
    "forces.all",
    "time.all",
)

# Expected line:
#
# Timestep: 123
#
TIMESTEP_RE = re.compile(
    r"^\s*Timestep:\s*(-?\d+)\s*$"
)


# ============================================================
# READ INPUT FILE
# ============================================================

def load_inputs():

    if not INPUTS_FILE.exists():

        raise FileNotFoundError(
            f"Cannot find inputs.json at {INPUTS_FILE}"
        )

    with INPUTS_FILE.open("r") as f:

        return json.load(f)


# ============================================================
# READ TIMESTEPS FROM A FILE
# ============================================================

def timestep_set(path):

    """
    Return all timestep values found in one trajectory file.
    """

    if not path.exists():

        return set()

    timesteps = set()

    with path.open(
        "r",
        errors="replace"
    ) as f:

        for line in f:

            match = TIMESTEP_RE.match(line)

            if match:

                timesteps.add(
                    int(match.group(1))
                )

    return timesteps


# ============================================================
# FIND COMMON TIMESTEP
# ============================================================

def latest_common_timestep(rep_dir):

    """
    Find the latest timestep which exists in ALL required
    trajectory files.

    Required:

        xyz.all
        momenta.all
        forces.all
        time.all
    """

    output_dir = rep_dir / "output"

    timestep_sets = {}

    for filename in TRAJECTORY_FILES:

        timestep_sets[filename] = timestep_set(
            output_dir / filename
        )

    # --------------------------------------------------------
    # Intersection of all four timestep sets.
    # --------------------------------------------------------

    common = set.intersection(
        *(
            timestep_sets[name]
            for name in TRAJECTORY_FILES
        )
    )

    if common:

        latest_common = max(common)

    else:

        latest_common = None

    # --------------------------------------------------------
    # Latest timestep in each individual file.
    # Useful for debugging.
    # --------------------------------------------------------

    latest_each = {}

    for name, values in timestep_sets.items():

        if values:

            latest_each[name] = max(values)

        else:

            latest_each[name] = None

    return latest_common, latest_each


# ============================================================
# FIND REQUESTED FINAL TIMESTEP
# ============================================================

def requested_final_timestep(inputs):

    """
    Read:

        Propagation.Timestep
        Propagation.Tot_timesteps

    and calculate:

        final timestep =
            Timestep * Tot_timesteps
    """

    propagation = inputs.get(
        "Propagation"
    )

    if not isinstance(propagation, dict):

        raise KeyError(
            "Propagation section is missing from inputs.json."
        )

    if (
        "Tot_timesteps" not in propagation
        or "Timestep" not in propagation
    ):

        raise KeyError(
            "Propagation.Tot_timesteps or "
            "Propagation.Timestep is missing."
        )

    timestep = int(
        propagation["Timestep"]
    )

    total_steps = int(
        propagation["Tot_timesteps"]
    )

    final_timestep = (
        timestep
        *
        total_steps
    )

    return final_timestep


# ============================================================
# MAXIMUM RESTARTS
# ============================================================

def get_max_restarts(inputs):

    """
    Read HPC.initre.

    This controls how many automatic restart attempts
    are allowed.
    """

    value = inputs["HPC"].get(
        "initre",
        0
    )

    try:

        return int(value)

    except (TypeError, ValueError):

        raise ValueError(
            "HPC.initre must be an integer."
        )


# ============================================================
# NUMBER OF REPETITIONS
# ============================================================

def get_repeats(inputs):

    try:

        repeats = int(
            inputs["HPC"]["repeats"]
        )

    except (
        KeyError,
        TypeError,
        ValueError
    ):

        raise ValueError(
            "HPC.repeats must be a valid integer."
        )

    if repeats < 1:

        raise ValueError(
            "HPC.repeats must be at least 1."
        )

    return repeats


# ============================================================
# CREATE RESTART ARRAY SCRIPT
# ============================================================

def write_restart_script(
    inputs,
    current_stage,
    repetition_ids
):

    """
    Create the next Plasma array.

    Example:

        repetition_ids = [2, 5]

    gives:

        #SBATCH --array=2,5

    so ONLY rep-2 and rep-5 are restarted.
    """

    next_stage = current_stage + 1

    runfolder = inputs["HPC"]["Runfolder"]

    script_name = (
        RUN_DIR
        /
        f"Plasma_{runfolder}_{next_stage}.sh"
    )

    array_spec = ",".join(
        str(rep)
        for rep in repetition_ids
    )

    python_exe = sys.executable

    with script_name.open("w") as f:

        f.write("#!/bin/bash\n")

        f.write(
            f"#SBATCH --job-name="
            f"Plasma_{runfolder}_{next_stage}\n"
        )

        f.write(
            f"#SBATCH --mem="
            f"{inputs['HPC']['memory']}G\n"
        )

        f.write(
            f"#SBATCH --time="
            f"{inputs['HPC']['hours']}:00:00\n"
        )

        f.write(
            f"#SBATCH --cpus-per-task="
            f"{inputs['HPC']['cores']}\n"
        )

        f.write(
            f"#SBATCH --array={array_spec}\n"
        )

        f.write(
            f"#SBATCH --output="
            f"{LOG_DIR}/"
            f"plasma_{next_stage}_%A_%a.out\n"
        )

        f.write(
            f"#SBATCH --error="
            f"{LOG_DIR}/"
            f"plasma_{next_stage}_%A_%a.err\n"
        )

        f.write(
            "set -euo pipefail\n"
        )

        f.write(
            'echo "HOSTNAME=$(hostname)"\n'
        )

        f.write(
            f'echo "PYTHON={python_exe}"\n'
        )

        f.write(
            'echo "TASK_ID=$SLURM_ARRAY_TASK_ID"\n'
        )

        f.write(
            f'echo "RESTART_STAGE={next_stage}"\n'
        )

        f.write(
            'echo "WORKDIR=$(pwd)"\n'
        )

        f.write(
            'export TMPDIR="${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}"\n'
        )

        f.write(
            'mkdir -p "$TMPDIR"\n'
        )

        # ----------------------------------------------------
        # Go to the correct repetition.
        # ----------------------------------------------------

        f.write(
            f'cd "{RUN_DIR}/repetitions/'
            '$SLURM_ARRAY_TASK_ID"\n'
        )

        # ----------------------------------------------------
        # Q-Chem environment.
        # ----------------------------------------------------

        if (
            inputs
            .get("Elec_structure", {})
            .get("method")
            == "QChem"
        ):

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

        # ----------------------------------------------------
        # IMPORTANT:
        #
        # Do NOT check .bomd_complete here.
        #
        # The validator already decided this repetition
        # is incomplete.
        # ----------------------------------------------------

        f.write(
            "which python || true\n"
        )

        f.write(
            f'"{python_exe}" "{MAIN_PY}"\n'
        )

        f.write(
            'echo "Restart stage completed."\n'
        )

    script_name.chmod(
        0o755
    )

    return script_name


# ============================================================
# SUBMIT SLURM COMMAND
# ============================================================

def submit(command):

    result = subprocess.check_output(
        command,
        text=True
    ).strip()

    # sbatch --parsable normally gives:
    #
    # 123456
    #
    # or:
    #
    # 123456;cluster

    return result.split(";")[0]


# ============================================================
# SUBMIT NEXT VALIDATOR
# ============================================================

def submit_next_validator(
    restart_job_id,
    next_stage
):

    """
    Submit another validator after the restart array.

    afterany is intentional.

    The validator must run after:

        COMPLETED
        TIMEOUT
        FAILED
        CANCELLED
    """

    validator = (
        RUN_DIR
        /
        "restart_validation.py"
    )

    LOG_DIR.mkdir(
        exist_ok=True
    )

    command = [

        "sbatch",

        "--parsable",

        f"--dependency=afterany:{restart_job_id}",

        (
            "--export=ALL,"
            f"VALIDATION_STAGE={next_stage}"
        ),

        f"--job-name=Validate_{next_stage}",

        "--time=00:30:00",

        "--mem=1G",

        "--cpus-per-task=1",

        (
            f"--output={LOG_DIR}/"
            "restart_validation_%j.out"
        ),

        (
            f"--error={LOG_DIR}/"
            "restart_validation_%j.err"
        ),

        str(validator),
    ]

    job_id = submit(
        command
    )

    print(
        "Next validator submitted:",
        job_id
    )

    return job_id


# ============================================================
# MAIN VALIDATION
# ============================================================

def main():

    inputs = load_inputs()

    repeats = get_repeats(
        inputs
    )

    max_restarts = get_max_restarts(
        inputs
    )

    target = requested_final_timestep(
        inputs
    )

    # --------------------------------------------------------
    # Validation stage.
    #
    # Initial Plasma stage = 1
    #
    # Validator after stage 1:
    #
    # VALIDATION_STAGE=1
    #
    # If incomplete:
    #     submit Plasma stage 2
    #
    # Then stage-2 validator is submitted.
    # --------------------------------------------------------

    stage = int(
        os.environ.get(
            "VALIDATION_STAGE",
            "1"
        )
    )

    print("=" * 72)

    print(
        "PLASMA TRAJECTORY RESTART VALIDATION"
    )

    print("=" * 72)

    print(
        "Validator:",
        Path(__file__).resolve()
    )

    print(
        "Run directory:",
        RUN_DIR
    )

    print(
        "Code directory:",
        CODE_DIR
    )

    print(
        "main.py:",
        MAIN_PY
    )

    print(
        "Validation stage:",
        stage
    )

    print(
        "Requested final timestep:",
        target
    )

    print(
        "Number of repetitions:",
        repeats
    )

    print(
        "Maximum restart attempts:",
        max_restarts
    )

    print("=" * 72)

    incomplete = []

    complete = []

    # ========================================================
    # CHECK EVERY REPETITION
    # ========================================================

    for rep in range(
        1,
        repeats + 1
    ):

        rep_dir = (
            REPETITIONS_DIR
            /
            f"rep-{rep}"
        )

        latest, latest_each = (
            latest_common_timestep(
                rep_dir
            )
        )

        # ----------------------------------------------------
        # COMPLETE
        # ----------------------------------------------------

        if (
            latest is not None
            and latest >= target
        ):

            complete.append(
                rep
            )

            status = "COMPLETE"

        # ----------------------------------------------------
        # INCOMPLETE
        # ----------------------------------------------------

        else:

            incomplete.append(
                rep
            )

            status = "INCOMPLETE"

        print(
            f"rep-{rep}: "
            f"{status:<10} "
            f"common={latest} "
            f"target={target}"
        )

        print(
            "          "
            f"xyz={latest_each['xyz.all']} "
            f"momenta={latest_each['momenta.all']} "
            f"forces={latest_each['forces.all']} "
            f"time={latest_each['time.all']}"
        )

    print("-" * 72)

    print(
        "Complete:",
        complete
    )

    print(
        "Incomplete:",
        incomplete
    )

    print("-" * 72)

    # ========================================================
    # EVERYTHING FINISHED
    # ========================================================

    if not incomplete:

        print(
            "ALL TRAJECTORIES REACHED "
            "THE REQUESTED FINAL TIMESTEP."
        )

        print(
            "NO RESTART ARRAY WILL BE SUBMITTED."
        )

        print(
            "VALIDATION LOOP FINISHED SUCCESSFULLY."
        )

        return 0

    # ========================================================
    # MAXIMUM RESTARTS
    # ========================================================

    if stage > max_restarts:

        print(
            "Maximum automatic restart attempts "
            f"({max_restarts}) reached."
        )

        print(
            "No further restart will be submitted."
        )

        return 2

    # ========================================================
    # RESTART ONLY INCOMPLETE REPETITIONS
    # ========================================================

    print(
        f"Submitting restart attempt {stage}"
    )

    print(
        "Only these repetitions will restart:",
        incomplete
    )

    restart_script = (
        write_restart_script(
            inputs,
            stage,
            incomplete
        )
    )

    print(
        "Restart script:",
        restart_script
    )

    restart_job = submit(
        [
            "sbatch",
            "--parsable",
            str(restart_script)
        ]
    )

    print(
        "Restart array submitted:",
        restart_job
    )

    # ========================================================
    # CHAIN ANOTHER VALIDATOR
    # ========================================================

    submit_next_validator(
        restart_job,
        stage + 1
    )

    print(
        "Restart/validation chain "
        "successfully extended."
    )

    return 0


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":

    try:

        sys.exit(
            main()
        )

    except Exception as exc:

        print(
            "ERROR:",
            exc,
            file=sys.stderr
        )

        sys.exit(1)
