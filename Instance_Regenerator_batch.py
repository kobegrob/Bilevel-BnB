import random
from pathlib import Path


def read_instance_file(instance_path):
    """
    Expected format:
      line 1: n
      line 2: b1
      line 3: b2
      line 4: a1_values (n integers)
      line 5: lowers    (n integers)
      line 6: fixed_costs (n integers)

    Any trailing non-numeric lines are ignored.
    """
    path = Path(instance_path)

    with path.open('r') as f:
        numeric_lines = []
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if all(part.lstrip('-').isdigit() for part in parts):
                numeric_lines.append([int(part) for part in parts])

    if len(numeric_lines) < 6:
        raise ValueError(
            "The instance file must contain at least 6 numeric lines: "
            "n, b1, b2, a1_values, lowers, fixed_costs."
        )

    n = numeric_lines[0][0]
    b1 = numeric_lines[1][0]
    b2 = numeric_lines[2][0]
    a1_values = numeric_lines[3]
    lowers = numeric_lines[4]
    fixed_costs = numeric_lines[5]

    for name, values in {
        'a1_values': a1_values,
        'lowers': lowers,
        'fixed_costs': fixed_costs,
    }.items():
        if len(values) != n:
            raise ValueError(
                f"{name} must contain exactly {n} values, but found {len(values)}."
            )

    return n, b1, b2, a1_values, lowers, fixed_costs



def generate_files_from_instance(base_filename, instance_path, knp_dir, mps_dir):
    n, b1, b2, a1_values, lowers, fixed_costs = read_instance_file(instance_path)

    # Values not specified in the input file are still generated as in the original file.
    capacities = []
    lowers1 = []

    for low in lowers:
        if random.random() > 0.75:   
            capacities.append((low // 2) + low)
            lowers1.append((low // 2) + 1)
        else:                       
            capacities.append(low)
            lowers1.append(low)

    fixed_costs1 = [((fix*10) // 2) + 1 for fix in fixed_costs]
    variable_costs = [((fix*10) // (2*low)) + 1 for fix in fixed_costs for low in lowers]

    items = list(zip(a1_values, fixed_costs1, capacities, lowers1, variable_costs))
    base_name = Path(base_filename).name
    # 1. Generate .knp File
    knp_path = knp_dir / f"{base_name}.knp"
    with open(knp_path, 'w') as f:
        f.write(f"{n}\t{b1}\t{b2}\n\n")
        for item in items:
            f.write("\t".join(map(str, item)) + "\n")

    # 2. Generate .mps File
    mps_path = mps_dir / f"{base_name}.mps"
    with open(mps_path, 'w') as f:
        f.write(f"NAME          {Path(base_filename).name}\n")
        f.write("ROWS\n L  UP_BUD\n")
        for i in range(1, n + 1):
            f.write(f" L  LINK{i}\n")
        f.write(" L  LO_CAP\n")
        for i in range(1, n + 1):
            f.write(f" L  LB_{i}\n L  UB_{i}\n")
        f.write(" N  OBJ\nCOLUMNS\n")
        f.write("    INT1      'MARKER'                 'INTORG'\n")

        # X variables (Binary Selection Indicators)
        for i, (a1, fixed, cap, low, var) in enumerate(items, 1):
            f.write(f"    X{i:<9}UP_BUD            {a1}\n")
            f.write(f"    X{i:<9}LINK{i:<12}1\n")

        # Z variables (Fixed Charge Indicators)
        for i, (a1, fixed, cap, low, var) in enumerate(items, 1):
            f.write(f"    Z{i:<9}LINK{i:<12}1\n")
            f.write(f"    Z{i:<9}LB_{i:<14} {low}\n")
            f.write(f"    Z{i:<9}UB_{i:<12}-{cap}\n")
            f.write(f"    Z{i:<9}OBJ               {fixed}\n")

        # Y variables (Continuous Quantities)
        for i, (a1, fixed, cap, low, var) in enumerate(items, 1):
            f.write(f"    Y{i:<9}LO_CAP            1\n")
            f.write(f"    Y{i:<9}LB_{i:<14}-1\n")
            f.write(f"    Y{i:<9}UB_{i:<14}1\n")
            f.write(f"    Y{i:<9}OBJ               {var}\n")

        f.write("    INT1END   'MARKER'                 'INTEND'\n")

        # RHS Section
        f.write("RHS\n")
        f.write(f"    B          UP_BUD           {b1}\n")
        for i in range(1, n + 1):
            f.write(f"    B          LINK{i:<14}1\n")
        f.write(f"    B          LO_CAP           {b2}\n")

        # BOUNDS Section
        f.write("BOUNDS\n")
        for i in range(1, n + 1):
            f.write(f" UP BOUND      X{i:<17}1\n")
        for i in range(1, n + 1):
            f.write(f" UP BOUND      Z{i:<17}1\n")
        for i, (a1, fixed, cap, low, var) in enumerate(items, 1):
            f.write(f" UP BOUND      Y{i:<17}{cap}\n")
        f.write("ENDATA\n")

    # 3. Generate .aux File
    aux_path = mps_dir / f"{base_name}.aux"
    with open(aux_path, 'w') as f:
        f.write(f"N {n * 2}\nM {1 + 3 * n}\n")
        for i in range(n, n + (n * 2)):
            f.write(f"LC {i}\n")
        for i in range(1, (1 + 3 * n) + 1):
            f.write(f"LR {i}\n")
        for i in range(n):
            f.write(f"LO -{fixed_costs1[i]}\n")
        for i in range(n):
            f.write(f"LO -{variable_costs[i]}\n")
        f.write("OS 1\n")



def generate_files_for_folder(folder_path):
    folder = Path(folder_path)

    if not folder.is_dir():
        raise NotADirectoryError(f"Folder not found: {folder}")

    # Create output folders
    knp_dir = folder / "new_knp"
    mps_dir = folder / "new_mps"

    knp_dir.mkdir(exist_ok=True)
    mps_dir.mkdir(exist_ok=True)

    ki_files = sorted(folder.glob('*.ki'))

    if not ki_files:
        print(f"No .ki files found in {folder}")
        return

    for ki_file in ki_files:
        output_base = ki_file.with_suffix('')
        generate_files_from_instance(output_base, ki_file, knp_dir, mps_dir)
        print(f"Processed: {ki_file.name}")


if __name__ == '__main__':
    input_folder = Path('/home/kobe/Documents/BnB/Bilevel-BnB/Density_Checks')
    generate_files_for_folder(input_folder)
