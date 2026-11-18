import numpy as np
import json
import pubchempy as pcp
import init
from scipy.spatial.distance import pdist, squareform
import os
import csv


# =====================================================================
# EXISTING MOMENTUM SAMPLING FUNCTION
# =====================================================================
#
# IMPORTANT:
# This function is intentionally kept the same as your existing code.
#
# It generates the random thermal momentum samples using:
#
#   1. Random Gaussian numbers
#   2. Normal-mode velocities
#   3. Transformation to Cartesian velocities
#   4. Conversion of velocity to momentum
#
# Therefore, enabling trajectory selection DOES NOT introduce a new
# momentum-sampling method.
#
# =====================================================================

def create_geom(n, nmod, T, modes, m, mom_num):

    Ax = modes[:, 0]
    Ay = modes[:, 1]
    Az = modes[:, 2]

    Ax = Ax.reshape(n, nmod, order='F')
    Ay = Ay.reshape(n, nmod, order='F')
    Az = Az.reshape(n, nmod, order='F')

    # ---------------------------------------------------------------
    # Generate independent Gaussian random numbers for every
    # normal mode and every momentum sample.
    #
    # Shape:
    #     nmod x mom_num
    #
    # Example:
    #     15 modes x 1,000,000 samples
    # ---------------------------------------------------------------

    rn = np.random.randn(nmod, mom_num)

    # Convert temperature from Kelvin to atomic units.
    T = T * 0.0000031668

    # ---------------------------------------------------------------
    # Calculate effective mass and normal-mode velocities.
    # ---------------------------------------------------------------

    Meff = np.zeros(nmod)

    rv = np.zeros((nmod, mom_num))

    for i in range(nmod):

        for j in range(n):

            Meff[i] = Meff[i] + np.sum(
                (
                    (Ax[j, i] ** 2)
                    + (Ay[j, i] ** 2)
                    + (Az[j, i] ** 2)
                ) * m[j]
            )

        rv[i, :] = (
            rn[i, :]
            * np.sqrt(2 * T / Meff[i])
        )

    # ---------------------------------------------------------------
    # Transform normal-mode velocities into Cartesian velocities.
    # ---------------------------------------------------------------

    Vx = np.dot(Ax, rv)
    Vy = np.dot(Ay, rv)
    Vz = np.dot(Az, rv)

    # ---------------------------------------------------------------
    # Convert Cartesian velocities into momenta:
    #
    #       p = m v
    #
    # Px, Py and Pz have dimensions:
    #
    #       number_of_atoms x number_of_samples
    # ---------------------------------------------------------------

    Px = np.zeros((n, mom_num))
    Py = np.zeros((n, mom_num))
    Pz = np.zeros((n, mom_num))

    for i in range(n):

        Px[i, :] = Vx[i, :] * m[i]
        Py[i, :] = Vy[i, :] * m[i]
        Pz[i, :] = Vz[i, :] * m[i]

    return Px, Py, Pz


# =====================================================================
# READ INITIAL GEOMETRY
# =====================================================================

def get_geometry_file():

    symbols = []
    coordinates = []

    with open("initial_guess.txt", "r") as file:

        for line in file:

            parts = line.split()

            if len(parts) == 4:

                symbols.append(parts[0])

                x, y, z = map(
                    float,
                    parts[1:]
                )

                # Convert Angstrom to Bohr.
                coordinates.append([
                    x * 1.88973,
                    y * 1.88973,
                    z * 1.88973
                ])

    return symbols, coordinates


# =====================================================================
# MOMENTUM CHECKS
# =====================================================================

def momenta_checks(
        Px,
        Py,
        Pz,
        symbols,
        temperature,
        output_file="momentum_checks.txt"):

    """
    Calculate kinetic-energy checks for the generated momentum samples.

    This is the original checking procedure.

    It is retained for the normal BOMD workflow.
    """

    natoms = Px.shape[0]
    repeats = Px.shape[1]

    masses = init.setup_masses(symbols)

    kinetic_energies = np.zeros(repeats)

    for j in range(repeats):

        file_kinetic_energy = 0.0

        for i in range(natoms):

            magnitude = np.sqrt(
                Px[i, j] ** 2
                + Py[i, j] ** 2
                + Pz[i, j] ** 2
            )

            kinetic_energy = (
                0.5
                * (magnitude ** 2 / masses[i])
            )

            file_kinetic_energy += kinetic_energy

        kinetic_energies[j] = file_kinetic_energy

    average_kinetic_energy = np.mean(
        kinetic_energies
    )

    expected_ke = (
        (3 * natoms - 6)
        * 3.16881e-6
        * temperature
    )

    avg_momentum_magnitude = np.zeros(
        natoms
    )

    for i in range(natoms):

        total_momentum_magnitude = 0.0

        for j in range(repeats):

            total_momentum_magnitude += np.sqrt(
                Px[i, j] ** 2
                + Py[i, j] ** 2
                + Pz[i, j] ** 2
            )

        avg_momentum_magnitude[i] = (
            total_momentum_magnitude / repeats
        )

    with open(
        output_file,
        "w"
    ) as file:

        file.write(
            "Kinetic Energy Checks:\n"
        )

        file.write(
            f"Expected Kinetic Energy: "
            f"{expected_ke:.6f}\n"
        )

        file.write(
            f"Average Kinetic Energy: "
            f"{average_kinetic_energy:.6f}\n"
        )

        file.write(
            "Average Momentum Magnitude per Atom:\n"
        )

        for i in range(natoms):

            file.write(
                f"Atom {i+1} ({symbols[i]}): "
                f"{avg_momentum_magnitude[i]:.6f}\n"
            )

        for j in range(repeats):

            file.write(
                f"Repeat {j+1}:\n"
            )

            file.write(
                f"  Kinetic Energy: "
                f"{kinetic_energies[j]:.6f}\n"
            )

            if np.isclose(
                kinetic_energies[j],
                expected_ke,
                rtol=0.1 * expected_ke
            ):

                file.write(
                    "  Kinetic Energy check passed.\n"
                )

            else:

                file.write(
                    "  Kinetic Energy check failed.\n"
                )

        file.write("\n")


# =====================================================================
# GENERATE INITIAL BOND INFORMATION
# =====================================================================

def generate_bondarr(symbols, coordinates):

    BOND_THRESHOLDS = {

        ("C", "C"): [
            (2.7401, "C=C"),
            (3.25, "C-C")
        ],

        ("C", "F"): [
            (2.8346, "C-F")
        ],

        ("C", "H"): [
            (2.1732, "C-H")
        ],

        ("C", "O"): [
            (2.75, "C-O")
        ],

        ("O", "H"): [
            (2.0, "O-H")
        ],

        ("H", "O"): [
            (2.0, "O-H")
        ]
    }

    natoms = len(symbols)

    distances = squareform(
        pdist(coordinates)
    )

    with open(
        "../results/bondarr.txt",
        "w"
    ) as file:

        for i in range(natoms):

            for j in range(i + 1, natoms):

                elem1 = symbols[i]
                elem2 = symbols[j]

                d = distances[i, j]

                key = tuple(
                    sorted(
                        (elem1, elem2)
                    )
                )

                if key in BOND_THRESHOLDS:

                    for max_dist, bond_type in (
                        BOND_THRESHOLDS[key]
                    ):

                        if d <= max_dist:

                            file.write(
                                f"{i+1}-{j+1}:"
                                f"{bond_type}\n"
                            )

    print(
        "bondarr.txt successfully written!"
    )


# =====================================================================
# GET GEOMETRY FROM PUBCHEM
# =====================================================================

def get_geometry_pubchem(molecule_name):

    search = pcp.get_compounds(
        molecule_name,
        'name',
        record_type='3d'
    )

    molecule = search[0]

    bonds = molecule.bonds
    atoms = molecule.atoms

    atom_id_to_symbol = {
        atom.aid: atom.element
        for atom in atoms
    }

    symbols = [
        atom.element
        for atom in atoms
    ]

    coordinates = [
        [
            atom.x * 1.88973,
            atom.y * 1.88973,
            atom.z * 1.88973
        ]
        for atom in atoms
    ]

    with open(
        "../results/bondarr.txt",
        "w"
    ) as file:

        for bond in bonds:

            atom1 = bond.aid1
            atom2 = bond.aid2

            order = bond.order

            if order == 1:
                bond_type = "-"

            elif order == 2:
                bond_type = "="

            elif order == 3:
                bond_type = "#"

            else:
                bond_type = "-"

            atom1_symbol = (
                atom_id_to_symbol.get(
                    atom1,
                    "Unknown"
                )
            )

            atom2_symbol = (
                atom_id_to_symbol.get(
                    atom2,
                    "Unknown"
                )
            )

            file.write(
                f"{atom1}-{atom2}:"
                f"{atom1_symbol}{bond_type}"
                f"{atom2_symbol}\n"
            )

    return symbols, coordinates


# =====================================================================
# ORGANISE NORMAL MODES
# =====================================================================

def organise_modes(modes, atoms):

    numeric_modes = np.zeros(
        (len(modes) * atoms, 3)
    )

    cnt = 0

    for i in range(len(modes)):

        for j in range(atoms):

            numeric_modes[cnt, 0] = float(
                modes[i]['displacement'][j][0]
            )

            numeric_modes[cnt, 1] = float(
                modes[i]['displacement'][j][1]
            )

            numeric_modes[cnt, 2] = float(
                modes[i]['displacement'][j][2]
            )

            cnt += 1

    return numeric_modes


# =====================================================================
# NEW TRAJECTORY-SELECTION FUNCTION
# =====================================================================

def get_selection_bonds(category):

    """
    Determine which bonds should be examined for one category.

    Example:

        {
            "name": "C-C",
            "bonds": ["2-4", "4-7"],
            "select": 10
        }

    means only C2-C4 and C4-C7 are considered.

    Example:

        {
            "name": "C-F",
            "bonds": [],
            "select": 10
        }

    means all C-F bonds in bondarr.txt are considered.

    NOTE:
    The bond labels are taken directly from bondarr.txt.
    """

    name = category["name"]

    requested_bonds = category["bonds"]

    # ---------------------------------------------------------------
    # Specific bonds supplied by the user.
    # ---------------------------------------------------------------

    if len(requested_bonds) > 0:

        bonds = []

        for bond in requested_bonds:

            parts = bond.split("-")

            if len(parts) != 2:
                continue

            i = int(parts[0])
            j = int(parts[1])

            bonds.append(
                (i, j)
            )

        return bonds

    # ---------------------------------------------------------------
    # Empty [] means ALL bonds of this type.
    # ---------------------------------------------------------------

    bonds = []

    bond_file = (
        "../results/bondarr.txt"
    )

    if not os.path.exists(
        bond_file
    ):

        return bonds

    with open(
        bond_file,
        "r"
    ) as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            parts = line.split(":")

            if len(parts) != 2:
                continue

            atom_part = parts[0]

            bond_type = parts[1]

            if bond_type != name:
                continue

            atom_numbers = (
                atom_part.split("-")
            )

            if len(atom_numbers) != 2:
                continue

            i = int(
                atom_numbers[0]
            )

            j = int(
                atom_numbers[1]
            )

            bonds.append(
                (i, j)
            )

    return bonds


# =====================================================================
# NEW TRAJECTORY-SELECTION FUNCTION
# =====================================================================

def calculate_vrad(
        Px,
        Py,
        Pz,
        coordinates,
        masses,
        i,
        j):

    """
    Calculate the initial radial velocity for ONE bond.

    Equations:

        v_i = p_i / m_i

        r_hat = (r_i-r_j) / |r_i-r_j|

        v_ij = v_i-v_j

        Vrad = v_ij . r_hat

    Positive Vrad:
        atoms are initially moving apart along the bond.

    Negative Vrad:
        atoms are initially moving towards each other.

    The function is vectorised over ALL momentum samples.

    Therefore, Px/Py/Pz contain:

        atoms x momentum_samples

    and the returned Vrad contains:

        momentum_samples
    """

    # ---------------------------------------------------------------
    # Convert atom numbering from:
    #
    #     user numbering: 1, 2, 3, ...
    #
    # to Python:
    #
    #     0, 1, 2, ...
    # ---------------------------------------------------------------

    i = i - 1
    j = j - 1

    # ---------------------------------------------------------------
    # Convert momentum to velocity.
    #
    # p = m v
    #
    # Therefore:
    #
    # v = p/m
    # ---------------------------------------------------------------

    vi_x = Px[i, :] / masses[i]
    vi_y = Py[i, :] / masses[i]
    vi_z = Pz[i, :] / masses[i]

    vj_x = Px[j, :] / masses[j]
    vj_y = Py[j, :] / masses[j]
    vj_z = Pz[j, :] / masses[j]

    # ---------------------------------------------------------------
    # Relative velocity.
    # ---------------------------------------------------------------

    vij_x = vi_x - vj_x
    vij_y = vi_y - vj_y
    vij_z = vi_z - vj_z

    # ---------------------------------------------------------------
    # Initial bond direction.
    # ---------------------------------------------------------------

    rij = (
        coordinates[i]
        - coordinates[j]
    )

    distance = np.linalg.norm(
        rij
    )

    if distance == 0:

        return np.zeros(
            Px.shape[1]
        )

    rhat = rij / distance

    # ---------------------------------------------------------------
    # Radial velocity:
    #
    # Vrad = (vi-vj) . r_hat
    # ---------------------------------------------------------------

    vrad = (
        vij_x * rhat[0]
        + vij_y * rhat[1]
        + vij_z * rhat[2]
    )

    return vrad


# =====================================================================
# CALCULATE SELECTION SCORES
# =====================================================================

def calculate_selection_scores(
        Px,
        Py,
        Pz,
        coordinates,
        masses,
        categories):

    """
    Calculate the best positive Vrad for every momentum sample.

    For each category:

        1. Find the requested bonds.
        2. Calculate Vrad for every bond.
        3. For each momentum sample, keep the largest Vrad.
        4. Record the bond producing that value.

    Therefore:

        C-C:
            sample 25
            C2-C4 = 0.0012
            C4-C7 = 0.0041

        gives:

            C-C Vrad = 0.0041
            bond = 4-7
    """

    number_of_samples = (
        Px.shape[1]
    )

    scores = {}

    for category in categories:

        name = category["name"]

        bonds = get_selection_bonds(
            category
        )

        # -----------------------------------------------------------
        # Initially every sample has no positive Vrad.
        # -----------------------------------------------------------

        category_scores = np.full(
            number_of_samples,
            -np.inf
        )

        category_bonds = np.full(
            number_of_samples,
            "",
            dtype=object
        )

        # -----------------------------------------------------------
        # Examine every bond in this category.
        # -----------------------------------------------------------

        for i, j in bonds:

            vrad = calculate_vrad(
                Px,
                Py,
                Pz,
                coordinates,
                masses,
                i,
                j
            )

            # -------------------------------------------------------
            # For every momentum sample keep the largest Vrad.
            # -------------------------------------------------------

            better = (
                vrad > category_scores
            )

            category_scores[better] = (
                vrad[better]
            )

            category_bonds[better] = (
                f"{i}-{j}"
            )

        scores[name] = {
            "vrad": category_scores,
            "bond": category_bonds
        }

    return scores


# =====================================================================
# SELECT MOMENTUM SAMPLES
# =====================================================================

def select_momentum_samples(
        scores,
        categories):

    """
    Select the requested number of momentum samples.

    Selection rule:

        - Vrad must be positive.
        - Larger Vrad has higher priority.
        - A momentum sample is used only once.

    If a sample is a strong candidate for more than one category,
    the category appearing first in inputs.json gets priority.

    This guarantees that, for example:

        10 C-C
        10 C-F
        10 C-H

    gives up to 30 DIFFERENT BOMD trajectories.
    """

    selected = []

    used_samples = set()

    for category in categories:

        name = category["name"]

        number_required = int(
            category["select"]
        )

        vrad_values = (
            scores[name]["vrad"]
        )

        bond_values = (
            scores[name]["bond"]
        )

        # -----------------------------------------------------------
        # Sort from largest Vrad to smallest.
        # -----------------------------------------------------------

        candidate_indices = np.argsort(
            vrad_values
        )[::-1]

        count = 0

        for sample in candidate_indices:

            vrad = vrad_values[sample]

            # -------------------------------------------------------
            # Only positive Vrad is allowed.
            # -------------------------------------------------------

            if vrad <= 0:
                continue

            # -------------------------------------------------------
            # Do not use the same momentum sample twice.
            # -------------------------------------------------------

            if int(sample) in used_samples:
                continue

            selected.append({
                "sample": int(sample),
                "category": name,
                "bond": str(
                    bond_values[sample]
                ),
                "vrad": float(vrad)
            })

            used_samples.add(
                int(sample)
            )

            count += 1

            if count >= number_required:
                break

    return selected


# =====================================================================
# CHECK HPC.REPEATS
# =====================================================================

def check_selection_repeats(
        inputs,
        categories):

    """
    Compare HPC.repeats with the total number requested by the
    trajectory-selection categories.

    IMPORTANT:

    A mismatch DOES NOT stop the calculation.

    A warning is written to:

        ../results/trajectory_selection_error.txt

    The values in Trajectory_selection are treated as authoritative.
    """

    requested_total = sum(
        int(category["select"])
        for category in categories
    )

    user_repeats = int(
        inputs["HPC"]["repeats"]
    )

    if user_repeats != requested_total:

        error_file = (
            "../results/"
            "trajectory_selection_error.txt"
        )

        with open(
            error_file,
            "w"
        ) as file:

            file.write(
                "=" * 78 + "\n"
            )

            file.write(
                "TRAJECTORY SELECTION WARNING\n"
            )

            file.write(
                "=" * 78 + "\n\n"
            )

            file.write(
                "HPC.repeats does not match "
                "Trajectory_selection.\n\n"
            )

            file.write(
                f"HPC.repeats        : "
                f"{user_repeats}\n"
            )

            file.write(
                f"Selection requested: "
                f"{requested_total}\n\n"
            )

            file.write(
                "The calculation will NOT be stopped.\n"
            )

            file.write(
                "The Trajectory_selection values will "
                "be used as the authoritative selection.\n"
            )

            file.write(
                f"Actual requested trajectory count: "
                f"{requested_total}\n"
            )

            file.write(
                "\n" + "=" * 78 + "\n"
            )

        print(
            "WARNING: HPC.repeats does not match "
            "Trajectory_selection."
        )

        print(
            f"HPC.repeats = {user_repeats}"
        )

        print(
            f"Selection requested = {requested_total}"
        )

        print(
            "Continuing using the selection criteria."
        )

    return requested_total


# =====================================================================
# WRITE HUMAN-READABLE SELECTION FILE
# =====================================================================

def write_selection_txt(
        selected,
        scores,
        categories,
        momentum_samples,
        temperature,
        output_file):

    """
    Write a human-readable record of the entire screening process.

    Contains:

        1. Temperature
        2. Number of momentum samples
        3. Selection criteria
        4. Selected trajectories
        5. Every sampled momentum ID
        6. Best Vrad for each category
        7. Bond producing that Vrad
        8. Whether the sample was selected
    """

    selected_lookup = {}

    for rep, item in enumerate(
        selected,
        start=1
    ):

        selected_lookup[
            item["sample"]
        ] = {
            "rep": rep,
            "category": item["category"],
            "bond": item["bond"],
            "vrad": item["vrad"]
        }

    with open(
        output_file,
        "w"
    ) as file:

        # -----------------------------------------------------------
        # HEADER
        # -----------------------------------------------------------

        file.write(
            "=" * 120 + "\n"
        )

        file.write(
            "TRAJECTORY MOMENTUM SELECTION RESULTS\n"
        )

        file.write(
            "=" * 120 + "\n\n"
        )

        file.write(
            f"Temperature           : "
            f"{temperature} K\n"
        )

        file.write(
            f"Momentum samples      : "
            f"{momentum_samples}\n"
        )

        file.write(
            f"Selected trajectories : "
            f"{len(selected)}\n\n"
        )

        # -----------------------------------------------------------
        # SELECTION CRITERIA
        # -----------------------------------------------------------

        file.write(
            "=" * 120 + "\n"
        )

        file.write(
            "SELECTION CRITERIA\n"
        )

        file.write(
            "=" * 120 + "\n\n"
        )

        for category in categories:

            name = category["name"]

            bonds = category["bonds"]

            number = category["select"]

            if len(bonds) == 0:

                bond_text = "ALL"

            else:

                bond_text = (
                    ", ".join(bonds)
                )

            file.write(
                f"{name:<10} "
                f"Bonds: {bond_text:<35} "
                f"Select: {number}\n"
            )

        # -----------------------------------------------------------
        # SELECTED TRAJECTORIES
        # -----------------------------------------------------------

        file.write("\n")

        file.write(
            "=" * 120 + "\n"
        )

        file.write(
            "SELECTED TRAJECTORIES\n"
        )

        file.write(
            "=" * 120 + "\n\n"
        )

        file.write(
            f"{'BOMD rep':<12}"
            f"{'Sample':<12}"
            f"{'Category':<12}"
            f"{'Bond':<12}"
            f"{'Positive Vrad':>20}\n"
        )

        file.write(
            "-" * 120 + "\n"
        )

        for rep, item in enumerate(
            selected,
            start=1
        ):

            file.write(
                f"{'rep-' + str(rep):<12}"
                f"{item['sample'] + 1:<12}"
                f"{item['category']:<12}"
                f"{item['bond']:<12}"
                f"{item['vrad']:>20.10e}\n"
            )

        # -----------------------------------------------------------
        # ALL MOMENTUM SAMPLES
        # -----------------------------------------------------------

        file.write("\n")

        file.write(
            "=" * 120 + "\n"
        )

        file.write(
            "ALL MOMENTUM SAMPLES\n"
        )

        file.write(
            "=" * 120 + "\n\n"
        )

        file.write(
            f"{'Sample':<12}"
        )

        for category in categories:

            name = category["name"]

            file.write(
                f"{name + ' Vrad':>20}"
                f"{name + ' Bond':>15}"
            )

        file.write(
            f"{'Selected':>12}\n"
        )

        file.write(
            "-" * 120 + "\n"
        )

        # -----------------------------------------------------------
        # Write every sampled momentum ID.
        # -----------------------------------------------------------

        for sample in range(
            momentum_samples
        ):

            file.write(
                f"{sample + 1:<12}"
            )

            for category in categories:

                name = category["name"]

                vrad = (
                    scores[name]["vrad"][sample]
                )

                bond = (
                    scores[name]["bond"][sample]
                )

                # Convert -inf into a readable value.
                if np.isneginf(vrad):

                    vrad_text = "N/A"

                else:

                    vrad_text = (
                        f"{vrad:.10e}"
                    )

                file.write(
                    f"{vrad_text:>20}"
                    f"{str(bond):>15}"
                )

            if sample in selected_lookup:

                file.write(
                    f"{'YES':>12}\n"
                )

            else:

                file.write(
                    f"{'NO':>12}\n"
                )

        # -----------------------------------------------------------
        # SUMMARY
        # -----------------------------------------------------------

        file.write("\n")

        file.write(
            "=" * 120 + "\n"
        )

        file.write(
            "SUMMARY\n"
        )

        file.write(
            "=" * 120 + "\n\n"
        )

        for category in categories:

            name = category["name"]

            requested = int(
                category["select"]
            )

            actual = sum(
                1
                for item in selected
                if item["category"] == name
            )

            file.write(
                f"{name}: "
                f"requested = {requested}, "
                f"selected = {actual}\n"
            )

        file.write(
            f"\nTotal selected trajectories = "
            f"{len(selected)}\n"
        )


# =====================================================================
# WRITE CSV FILE
# =====================================================================

def write_selection_csv(
        selected,
        scores,
        categories,
        momentum_samples,
        output_file):

    """
    Write the complete screening data into a CSV file.

    This file can be opened directly in Excel.

    One row = one momentum sample.
    """

    selected_lookup = {}

    for rep, item in enumerate(
        selected,
        start=1
    ):

        selected_lookup[
            item["sample"]
        ] = {
            "rep": rep,
            "category": item["category"],
            "bond": item["bond"]
        }

    with open(
        output_file,
        "w",
        newline=""
    ) as file:

        writer = csv.writer(
            file
        )

        # -----------------------------------------------------------
        # Header
        # -----------------------------------------------------------

        header = [
            "Sample"
        ]

        for category in categories:

            name = category["name"]

            header.append(
                f"{name}_Vrad"
            )

            header.append(
                f"{name}_bond"
            )

        header.extend([
            "Selected",
            "BOMD_rep",
            "Selected_category",
            "Selected_bond"
        ])

        writer.writerow(
            header
        )

        # -----------------------------------------------------------
        # Data
        # -----------------------------------------------------------

        for sample in range(
            momentum_samples
        ):

            row = [
                sample + 1
            ]

            for category in categories:

                name = category["name"]

                vrad = (
                    scores[name]["vrad"][sample]
                )

                bond = (
                    scores[name]["bond"][sample]
                )

                if np.isneginf(vrad):

                    vrad = ""

                row.append(
                    vrad
                )

                row.append(
                    bond
                )

            if sample in selected_lookup:

                info = selected_lookup[
                    sample
                ]

                row.extend([
                    "YES",
                    f"rep-{info['rep']}",
                    info["category"],
                    info["bond"]
                ])

            else:

                row.extend([
                    "NO",
                    "",
                    "",
                    ""
                ])

            writer.writerow(
                row
            )


# =====================================================================
# WRITE EXACT MOMENTA OF SELECTED TRAJECTORIES
# =====================================================================

def write_selected_momenta(
        selected,
        Px,
        Py,
        Pz,
        symbols,
        output_file):

    """
    Save the exact initial momenta that are actually passed to BOMD.

    This allows every selected trajectory to be traced back to its
    original momentum sample.
    """

    with open(
        output_file,
        "w"
    ) as file:

        file.write(
            "=" * 100 + "\n"
        )

        file.write(
            "INITIAL MOMENTA OF SELECTED BOMD TRAJECTORIES\n"
        )

        file.write(
            "=" * 100 + "\n\n"
        )

        for rep, item in enumerate(
            selected,
            start=1
        ):

            sample = item["sample"]

            file.write(
                "-" * 100 + "\n"
            )

            file.write(
                f"BOMD trajectory : rep-{rep}\n"
            )

            file.write(
                f"Original sample : {sample + 1}\n"
            )

            file.write(
                f"Category        : "
                f"{item['category']}\n"
            )

            file.write(
                f"Selected bond   : "
                f"{item['bond']}\n"
            )

            file.write(
                f"Positive Vrad   : "
                f"{item['vrad']:.10e}\n\n"
            )

            file.write(
                f"{'Atom':<8}"
                f"{'Element':<10}"
                f"{'Px':>22}"
                f"{'Py':>22}"
                f"{'Pz':>22}\n"
            )

            file.write(
                "-" * 100 + "\n"
            )

            for atom in range(
                len(symbols)
            ):

                file.write(
                    f"{atom + 1:<8}"
                    f"{symbols[atom]:<10}"
                    f"{Px[atom, sample]:>22.10e}"
                    f"{Py[atom, sample]:>22.10e}"
                    f"{Pz[atom, sample]:>22.10e}\n"
                )

            file.write("\n")


# =====================================================================
# MAIN SETUP PROCEDURE
# =====================================================================

if __name__ == "__main__":

    # ---------------------------------------------------------------
    # Read inputs.json
    # ---------------------------------------------------------------

    with open(
        "../inputs.json"
    ) as f:

        inputs = json.load(f)


    # ===============================================================
    # OBTAIN INITIAL GEOMETRY
    # ===============================================================

    if (
        inputs["Molecule_data"]["Geom_flg"]
        == "PubChem"
    ):

        symbols, coordinates = (
            get_geometry_pubchem(
                inputs["Molecule_data"]["Molecule"]
            )
        )

    elif (
        inputs["Molecule_data"]["Geom_flg"]
        == "Initial"
    ):

        symbols, coordinates = (
            get_geometry_file()
        )

        # -----------------------------------------------------------
        # Generate bondarr.txt from the initial molecular geometry.
        # -----------------------------------------------------------

        generate_bondarr(
            symbols,
            coordinates
        )


    # ===============================================================
    # INITIAL ELECTRONIC-STRUCTURE CALCULATION
    # ===============================================================

    if (
        inputs["Elec_structure"]["method"]
        == "QChem"
    ):

        import qchem as qc

        opt_coords, modes = (
            qc.initial_conditions(
                symbols,
                coordinates,
                inputs["HPC"]["cores"],
                inputs["Elec_structure"]["Basis"]
            )
        )

    elif (
        inputs["Elec_structure"]["method"]
        == "PySCF"
    ):

        import elec_pyscf as pyscf

        pyscf.initial_conditions(
            symbols,
            coordinates
        )


    # ===============================================================
    # ORGANISE NORMAL MODES
    # ===============================================================

    modes = modes["modes"]

    num_modes = len(modes)

    natoms = int(
        (num_modes + 6) / 3
    )

    modes = organise_modes(
        modes,
        natoms
    )

    masses = init.setup_masses(
        symbols
    )


    # ===============================================================
    # READ TRAJECTORY-SELECTION SETTINGS
    # ===============================================================

    selection = inputs.get(
        "Trajectory_selection",
        {}
    )

    selection_enabled = int(
        selection.get(
            "enabled",
            0
        )
    )


    # ===============================================================
    # MODE 1:
    #
    # TRAJECTORY SELECTION DISABLED
    #
    # This is the ORIGINAL BOMD behaviour.
    #
    # Nothing is screened.
    # HPC.repeats directly determines the number of momentum
    # samples and the number of trajectories.
    # ===============================================================

    if selection_enabled == 0:

        print(
            "Trajectory selection is DISABLED."
        )

        print(
            "Using the original BOMD workflow."
        )

        momentum_samples = int(
            inputs["HPC"]["repeats"]
        )

        # -----------------------------------------------------------
        # EXISTING MOMENTUM SAMPLING.
        #
        # This is exactly the same create_geom() function used
        # previously.
        # -----------------------------------------------------------

        Px, Py, Pz = create_geom(
            natoms,
            num_modes,
            inputs["Molecule_data"]["Temp"],
            modes,
            masses,
            momentum_samples
        )

        # -----------------------------------------------------------
        # Original momentum checks.
        # -----------------------------------------------------------

        momenta_checks(
            Px,
            Py,
            Pz,
            symbols,
            inputs["Molecule_data"]["Temp"]
        )

        # -----------------------------------------------------------
        # Write every sampled momentum as a Geometry file.
        # -----------------------------------------------------------

        for j in range(
            inputs["HPC"]["repeats"]
        ):

            with open(
                "../repetitions/rep-"
                + str(j + 1)
                + "/Geometry",
                "w"
            ) as file:

                file.write(
                    opt_coords
                )

                file.write(
                    "momentum\n"
                )

                for atom in range(
                    natoms
                ):

                    px_value = (
                        Px[atom, j]
                    )

                    py_value = (
                        Py[atom, j]
                    )

                    pz_value = (
                        Pz[atom, j]
                    )

                    file.write(
                        f"{px_value}  "
                        f"{py_value}  "
                        f"{pz_value}\n"
                    )

        print(
            "Original BOMD setup completed."
        )


    # ===============================================================
    # MODE 2:
    #
    # TRAJECTORY SELECTION ENABLED
    #
    # The code:
    #
    #   1. Generates a large number of momentum samples using the
    #      EXISTING create_geom() function.
    #
    #   2. Calculates positive Vrad for every requested bond.
    #
    #   3. Selects the highest-Vrad samples.
    #
    #   4. Writes only selected samples into rep-* directories.
    #
    #   5. Saves complete selection information.
    #
    # IMPORTANT:
    # No second momentum sampling is performed.
    # ===============================================================

    else:

        print(
            "Trajectory selection is ENABLED."
        )

        categories = (
            selection["categories"]
        )

        # -----------------------------------------------------------
        # Number of candidate momentum samples.
        #
        # Example:
        #
        #     1,000,000
        #
        # This is NOT the number of BOMD trajectories.
        # -----------------------------------------------------------

        momentum_samples = int(
            selection[
                "momentum_samples"
            ]
        )

        # -----------------------------------------------------------
        # Check HPC.repeats.
        #
        # A mismatch only produces a warning.
        # It does NOT stop the calculation.
        # -----------------------------------------------------------

        actual_repeats = (
            check_selection_repeats(
                inputs,
                categories
            )
        )

        print(
            f"Generating {momentum_samples} "
            "candidate momentum samples."
        )

        print(
            f"Requested selected trajectories: "
            f"{actual_repeats}"
        )

        # -----------------------------------------------------------
        # USE THE EXISTING MOMENTUM SAMPLING.
        #
        # This is the critical point:
        #
        # We are NOT creating a new sampling algorithm.
        #
        # The same random normal-mode sampling used by your original
        # BOMD code generates all candidate momentum distributions.
        # -----------------------------------------------------------

        Px, Py, Pz = create_geom(
            natoms,
            num_modes,
            inputs["Molecule_data"]["Temp"],
            modes,
            masses,
            momentum_samples
        )

        print(
            "Momentum sampling completed."
        )


        # ===========================================================
        # CONVERT OPTIMISED GEOMETRY TO NUMERICAL COORDINATES
        # ===========================================================

        coordinates = []

        for line in (
            opt_coords.splitlines()
        ):

            parts = line.split()

            if len(parts) >= 4:

                try:

                    coordinates.append([
                        float(parts[1]),
                        float(parts[2]),
                        float(parts[3])
                    ])

                except ValueError:

                    pass

        coordinates = np.array(
            coordinates,
            dtype=float
        )

        if len(coordinates) != natoms:

            print(
                "WARNING: Number of coordinates in "
                "optimised geometry does not match "
                "the number of atoms."
            )


        # ===========================================================
        # CALCULATE RADIAL VELOCITIES
        # ===========================================================

        print(
            "Calculating bond-directed radial velocities..."
        )

        scores = (
            calculate_selection_scores(
                Px,
                Py,
                Pz,
                coordinates,
                masses,
                categories
            )
        )

        print(
            "Vrad calculation completed."
        )


        # ===========================================================
        # SELECT THE MOMENTUM SAMPLES
        # ===========================================================

        selected = (
            select_momentum_samples(
                scores,
                categories
            )
        )

        print(
            f"Selected {len(selected)} "
            "momentum samples."
        )


        # ===========================================================
        # CREATE RESULTS DIRECTORY
        # ===========================================================

        os.makedirs(
            "../results",
            exist_ok=True
        )


        # ===========================================================
        # SAVE COMPLETE HUMAN-READABLE RESULTS
        # ===========================================================

        write_selection_txt(
            selected,
            scores,
            categories,
            momentum_samples,
            inputs["Molecule_data"]["Temp"],
            "../results/"
            "trajectory_selection.txt"
        )


        # ===========================================================
        # SAVE COMPLETE CSV FOR EXCEL
        # ===========================================================

        write_selection_csv(
            selected,
            scores,
            categories,
            momentum_samples,
            "../results/"
            "trajectory_selection.csv"
        )


        # ===========================================================
        # SAVE EXACT MOMENTA USED BY SELECTED TRAJECTORIES
        # ===========================================================

        write_selected_momenta(
            selected,
            Px,
            Py,
            Pz,
            symbols,
            "../results/"
            "selected_momenta.txt"
        )


        # ===========================================================
        # WRITE ONLY SELECTED MOMENTUM SAMPLES TO Geometry FILES
        # ===========================================================

        for rep, item in enumerate(
            selected,
            start=1
        ):

            # -------------------------------------------------------
            # Original sample number in the million-sample pool.
            # -------------------------------------------------------

            sample = item["sample"]

            rep_dir = (
                "../repetitions/rep-"
                + str(rep)
            )

            os.makedirs(
                rep_dir,
                exist_ok=True
            )

            geometry_file = os.path.join(
                rep_dir,
                "Geometry"
            )

            with open(
                geometry_file,
                "w"
            ) as file:

                # ---------------------------------------------------
                # Every trajectory starts from the same optimised
                # molecular geometry.
                # ---------------------------------------------------

                file.write(
                    opt_coords
                )

                file.write(
                    "momentum\n"
                )

                # ---------------------------------------------------
                # But every selected trajectory receives a different
                # initial momentum distribution.
                # ---------------------------------------------------

                for atom in range(
                    natoms
                ):

                    file.write(
                        f"{Px[atom, sample]}  "
                        f"{Py[atom, sample]}  "
                        f"{Pz[atom, sample]}\n"
                    )


        # ===========================================================
        # FINAL MESSAGE
        # ===========================================================

        print(
            "Selected Geometry files successfully created."
        )

        print(
            "Trajectory selection setup completed."
        )

        print(
            "Results saved to:"
        )

        print(
            "  ../results/trajectory_selection.txt"
        )

        print(
            "  ../results/trajectory_selection.csv"
        )

        print(
            "  ../results/selected_momenta.txt"
        )

        if os.path.exists(
            "../results/"
            "trajectory_selection_error.txt"
        ):

            print(
                "WARNING: Check:"
            )

            print(
                "  ../results/"
                "trajectory_selection_error.txt"
            )
