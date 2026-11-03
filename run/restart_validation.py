#!/usr/bin/env python3

"""
restart_validation.py

Validation-driven Slurm restart controller for PlasmaDissociation.

Workflow
--------

1. Read inputs.json automatically from the current run directory.

2. Check every repetition.

3. Find the latest timestep common to:

       output/xyz.all
       output/momenta.all
       output/forces.all
       output/time.all

4. Compare the common timestep with:

       Propagation.Tot_timesteps * Propagation.Timestep

5. Restart ONLY incomplete repetitions.

6. Submit another validator using an afterany dependency.

7. Continue the validation/restart loop until all repetitions reach
   the requested final timestep.

8. If the previous Plasma job was manually CANCELLED, terminate the
   restart chain.

9. If the previous Plasma job FAILED, terminate the restart chain.

10. TIMEOUT and OUT_OF_MEMORY are treated as restartable conditions.

.bomd_complete is NOT used as the completion criterion.

No VALIDATION_STAGE is used.
No restart number is used to decide whether to continue.
No HPC.initre limit is used.

The actual trajectory timestep is the completion criterion.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path


################################################################################
# AUTOMATIC PATH DISCOVERY
################################################################################

SCRIPT_DIR = Path(__file__).resolve().parent

RUN_DIR = SCRIPT_DIR

INPUTS_FILE = RUN_DIR / "inputs.json"

REPETITIONS_DIR = RUN_DIR / "repetitions"

LOG_DIR = RUN_DIR / "logs"


################################################################################
# FIND PROJECT CODE AUTOMATICALLY
################################################################################

def find_code_directory():

    """
    Automatically find a directory containing:

        code/main.py
        code/setup.py

    No hard-coded PlasmaDissociation path is required.
    """

    current = SCRIPT_DIR

    # First search upwards from the current run directory.
    for directory in [current] + list(current.parents):

        candidate = directory / "code"

        if (
            candidate.is_dir()
            and (candidate / "main.py").is_file()
            and (candidate / "setup.py").is_file()
        ):
            return candidate

    # Then search common locations in the user's home directory.
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

    # Last search under HOME.
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
        "Could not automatically find project code directory "
        "containing code/main.py and code/setup.py."
    )


CODE_DIR = find_code_directory()

MAIN_PY = CODE_DIR / "main.py"

SETUP_PY = CODE_DIR / "setup.py"


################################################################################
# TRAJECTORY FILES
################################################################################

TRAJECTORY_FILES = (
    "xyz.all",
    "momenta.all",
    "forces.all",
    "time.all",
)


################################################################################
# TIMESTEP REGEX
################################################################################

TIMESTEP_RE = re.compile(
    r"^\s*\*\*Timestep:\s*(-?\d+)\s*$"
)


################################################################################
# LOAD INPUTS
################################################################################

def load_inputs():

    if not INPUTS_FILE.exists():

        raise FileNotFoundError(
            f"Cannot find inputs.json at {INPUTS_FILE}"
        )

    with INPUTS_FILE.open("r") as f:

        return json.load(f)


################################################################################
# READ TIMESTEPS
################################################################################

def timestep_set(path):

    """
    Return all timestep values found in a trajectory file.
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


################################################################################
# FIND COMMON TIMESTEP
################################################################################

def latest_common_timestep(rep_dir):

    """
    Find the latest timestep which exists in ALL required
    trajectory files.
    """

    output_dir = rep_dir / "output"

    timestep_sets = {}

    for filename in TRAJECTORY_FILES:

        timestep_sets[filename] = timestep_set(
            output_dir / filename
        )

    # If any required file has no timestep, there can be no
    # common timestep.
    if any(
        not timestep_sets[name]
        for name in TRAJECTORY_FILES
    ):

        common = set()

    else:

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

    latest_each = {}

    for name, values in timestep_sets.items():

        if values:
            latest_each[name] = max(values)

        else:
            latest_each[name] = None

    return latest_common, latest_each


################################################################################
# REQUESTED FINAL TIMESTEP
################################################################################

def requested_final_timestep(inputs):

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

    return timestep * total_steps


################################################################################
# NUMBER OF REPETITIONS
################################################################################

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


################################################################################
# CREATE RESTART ARRAY SCRIPT
################################################################################

def write_restart_script(
    inputs,
    repetition_ids
):

    """
    Create the next Plasma array.

    There is NO stage number.

    The restart array contains only repetitions which have not
    reached the requested final timestep.
    """

    runfolder = inputs["HPC"]["Runfolder"]

    # Use the current validator job ID to make the generated
    # script filename unique.
    validator_id = os.environ.get(
        "SLURM_JOB_ID",
        "manual"
    )

    script_name = (
        RUN_DIR
        /
        f"Plasma_{runfolder}_restart_{validator_id}.sh"
    )

    array_spec = ",".join(
        str(rep)
        for rep in repetition_ids
    )

    python_exe = sys.executable

    with script_name.open("w") as f:

        f.write("#!/bin/bash\n")

        f.write(
            f"#SBATCH --job-name=Plasma_{runfolder}_restart\n"
        )

        f.write(
            f"#SBATCH --mem={inputs['HPC']['memory']}G\n"
        )

        f.write(
            f"#SBATCH --time={inputs['HPC']['hours']}:00:00\n"
        )

        f.write(
            f"#SBATCH --cpus-per-task={inputs['HPC']['cores']}\n"
        )

        f.write(
            f"#SBATCH --array={array_spec}\n"
        )

        f.write(
            f"#SBATCH --output={LOG_DIR}/"
            f"plasma_restart_%A_%a.out\n"
        )

        f.write(
            f"#SBATCH --error={LOG_DIR}/"
            f"plasma_restart_%A_%a.err\n"
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
            'echo "WORKDIR=$(pwd)"\n'
        )

        f.write(
            f'export PLASMA_INPUTS="{RUN_DIR / "inputs.json"}"\n'
        )

        f.write(
            'export TMPDIR="${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}"\n'
        )

        f.write(
            'mkdir -p "$TMPDIR"\n'
        )

        #######################################################################
        # Correct repetition directory
        #######################################################################

        f.write(
            f'cd "{RUN_DIR}/repetitions/rep-$SLURM_ARRAY_TASK_ID"\n'
        )

        #######################################################################
        # Q-Chem
        #######################################################################

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

        #######################################################################
        # IMPORTANT:
        #
        # Do NOT check .bomd_complete here.
        #
        # The validator has already determined which trajectories are
        # incomplete.
        #######################################################################

        f.write(
            "which python || true\n"
        )

        f.write(
            f'"{python_exe}" "{MAIN_PY}"\n'
        )

        f.write(
            'echo "Restart Plasma job finished."\n'
        )

    script_name.chmod(0o755)

    return script_name


################################################################################
# SUBMIT COMMAND
################################################################################

def submit(command):

    result = subprocess.check_output(
        command,
        text=True
    ).strip()

    return result.split(";")[0]


################################################################################
# CREATE VALIDATOR WRAPPER
################################################################################

def create_validator_wrapper():

    """
    Create the Slurm wrapper which runs restart_validation.py.

    No stage number is used.

    The previous Plasma job ID is supplied by SLURM through:

        PREVIOUS_PLASMA_JOB_ID
    """

    wrapper = (
        RUN_DIR
        /
        "restart_validation_job.sh"
    )

    python_exe = sys.executable

    with wrapper.open("w") as f:

        f.write("#!/bin/bash\n")

        f.write(
            "#SBATCH --job-name=RestartValidator\n"
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
            f"#SBATCH --output={LOG_DIR}/"
            "restart_validation_%j.out\n"
        )

        f.write(
            f"#SBATCH --error={LOG_DIR}/"
            "restart_validation_%j.err\n"
        )

        f.write(
            "set -euo pipefail\n"
        )

        f.write(
            f'export PLASMA_INPUTS="{RUN_DIR}/inputs.json"\n'
        )

        f.write(
            'echo "PREVIOUS_PLASMA_JOB_ID='
            '$PREVIOUS_PLASMA_JOB_ID"\n'
        )

        f.write(
            f'cd "{RUN_DIR}"\n'
        )

        f.write(
            f'"{python_exe}" '
            f'"{RUN_DIR / "restart_validation.py"}"\n'
        )

    wrapper.chmod(0o755)

    return wrapper


################################################################################
# SUBMIT NEXT VALIDATOR
################################################################################

def submit_next_validator(
    restart_job_id
):

    """
    Submit another validator after a Plasma restart job.

    The validator depends on the Plasma job using afterany.

    Therefore the validator runs after:

        COMPLETED
        TIMEOUT
        FAILED
        CANCELLED
        OUT_OF_MEMORY

    The validator then decides what should happen.

    No stage number is used.
    """

    validator_wrapper = create_validator_wrapper()

    command = [
        "sbatch",
        "--parsable",

        # Pass the actual Plasma job ID to the validator.
        "--export=ALL,PREVIOUS_PLASMA_JOB_ID="
        + str(restart_job_id),

        # Validator must run regardless of Plasma termination state.
        "--dependency=afterany:"
        + str(restart_job_id),

        str(validator_wrapper)
    ]

    job_id = submit(command)

    print(
        "Next validator submitted:",
        job_id,
        "after Plasma job:",
        restart_job_id
    )

    return job_id


################################################################################
# GET PREVIOUS SLURM JOB STATE
################################################################################

def get_previous_job_state(
    job_id
):

    """
    Return the overall state of the previous Plasma array.

    Important:

    This function does NOT use a stage number.

    For an array job there may be several states simultaneously.

    Examples:

        COMPLETED
        TIMEOUT
        CANCELLED
        FAILED
        OUT_OF_MEMORY

    Priority:

        CANCELLED
            -> stop

        TIMEOUT
            -> validate/restart

        OUT_OF_MEMORY
            -> validate/restart

        COMPLETED
            -> validate

        FAILED
            -> stop
    """

    if not job_id:
        return None

    try:

        result = subprocess.check_output(
            [
                "sacct",
                "-n",
                "-P",
                "-j",
                str(job_id),
                "--format=State"
            ],
            text=True,
            stderr=subprocess.DEVNULL
        )

        states = []

        for line in result.splitlines():

            state = line.strip().split("|")[0]

            if not state:
                continue

            # Remove Slurm suffixes such as:
            #
            # CANCELLED+
            # TIMEOUT+
            #
            state = state.split("+")[0]

            if state not in states:
                states.append(state)

        print(
            "Previous Plasma job states:",
            states
        )

        if not states:
            return None

        #######################################################################
        # MANUAL CANCELLATION HAS HIGHEST PRIORITY
        #######################################################################

        if "CANCELLED" in states:

            return "CANCELLED"

        #######################################################################
        # TIMEOUT ALLOWS AUTOMATIC VALIDATION/RESTART
        #######################################################################

        if "TIMEOUT" in states:

            return "TIMEOUT"

        #######################################################################
        # OOM ALLOWS AUTOMATIC VALIDATION/RESTART
        #######################################################################

        if "OUT_OF_MEMORY" in states:

            return "OUT_OF_MEMORY"

        #######################################################################
        # ALL COMPLETED
        #######################################################################

        if all(
            state == "COMPLETED"
            for state in states
        ):

            return "COMPLETED"

        #######################################################################
        # GENUINE FAILURE
        #######################################################################

        if "FAILED" in states:

            return "FAILED"

        #######################################################################
        # Otherwise return the first state.
        #######################################################################

        return states[0]

    except Exception as exc:

        print(
            "WARNING: Could not determine previous "
            "Plasma job state:",
            exc
        )

    return None


################################################################################
# MAIN VALIDATION
################################################################################

def main():

    ###########################################################################
    # LOAD INPUTS
    ###########################################################################

    inputs = load_inputs()

    repeats = get_repeats(inputs)

    target = requested_final_timestep(inputs)

    ###########################################################################
    # PREVIOUS PLASMA JOB
    ###########################################################################

    previous_plasma_job_id = os.environ.get(
        "PREVIOUS_PLASMA_JOB_ID"
    )

    previous_job_state = get_previous_job_state(
        previous_plasma_job_id
    )

    print(
        "Previous Plasma job ID:",
        previous_plasma_job_id
    )

    print(
        "Previous Plasma job state:",
        previous_job_state
    )

    ###########################################################################
    # MANUAL CANCELLATION
    ###########################################################################

    if previous_job_state == "CANCELLED":

        print("=" * 72)

        print(
            "PREVIOUS PLASMA JOB WAS CANCELLED."
        )

        print(
            "NO RESTART WILL BE SUBMITTED."
        )

        print(
            "VALIDATION LOOP TERMINATED."
        )

        print("=" * 72)

        return 3

    ###########################################################################
    # GENUINE FAILURE
    ###########################################################################

    if previous_job_state == "FAILED":

        print("=" * 72)

        print(
            "PREVIOUS PLASMA JOB FAILED."
        )

        print(
            "NO AUTOMATIC RESTART WILL BE SUBMITTED."
        )

        print(
            "VALIDATION LOOP TERMINATED."
        )

        print("=" * 72)

        return 4

    ###########################################################################
    # HEADER
    ###########################################################################

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
        "Previous Plasma job:",
        previous_plasma_job_id
    )

    print(
        "Previous Plasma state:",
        previous_job_state
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
        "Restart limit: NONE"
    )

    print(
        "Completion criterion: trajectory timestep"
    )

    print("=" * 72)

    ###########################################################################
    # CHECK EVERY REPETITION
    ###########################################################################

    incomplete = []

    complete = []

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

        #######################################################################
        # COMPLETE
        #######################################################################

        if (
            latest is not None
            and latest >= target
        ):

            complete.append(rep)

            status = "COMPLETE"

        #######################################################################
        # INCOMPLETE
        #######################################################################

        else:

            incomplete.append(rep)

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

    ###########################################################################
    # SUMMARY
    ###########################################################################

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

    ###########################################################################
    # EVERYTHING FINISHED
    ###########################################################################

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

    ###########################################################################
    # RESTART ONLY INCOMPLETE REPETITIONS
    ###########################################################################

    print(
        "Trajectory data are incomplete."
    )

    print(
        "Only these repetitions will restart:",
        incomplete
    )

    ###########################################################################
    # CREATE RESTART SCRIPT
    ###########################################################################

    restart_script = write_restart_script(
        inputs,
        incomplete
    )

    print(
        "Restart script:",
        restart_script
    )

    ###########################################################################
    # SUBMIT RESTART PLASMA ARRAY
    ###########################################################################

    restart_job = submit(
        [
            "sbatch",
            "--parsable",
            str(restart_script)
        ]
    )

    print(
        "Restart Plasma array submitted:",
        restart_job
    )

    ###########################################################################
    # CRITICAL:
    #
    # IMMEDIATELY SUBMIT ANOTHER VALIDATOR.
    #
    # The new validator receives:
    #
    #     PREVIOUS_PLASMA_JOB_ID=<restart_job>
    #
    # and depends on the restart Plasma job using afterany.
    #
    # This creates the continuous loop:
    #
    #     Plasma
    #        |
    #        v
    #     Validator
    #        |
    #        v
    #     Plasma
    #        |
    #        v
    #     Validator
    #        |
    #        ...
    #
    # The loop ends ONLY when:
    #
    #     all trajectories >= target
    #
    # or when:
    #
    #     Plasma job is CANCELLED
    #
    # or:
    #
    #     Plasma job genuinely FAILED
    ###########################################################################

    submit_next_validator(
        restart_job
    )

    print(
        "Restart/validation chain successfully extended."
    )

    return 0


################################################################################
# ENTRY POINT
################################################################################

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
