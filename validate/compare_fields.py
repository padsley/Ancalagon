#!/usr/bin/env python3
"""Compare a MitrayXField C++ port against the reference values produced
directly by the matching real GEANT3 mitray_*.f source. Works for any of
the probe CSVs (bx_T/by_T/bz_T for the magnetic elements, ex_kVcm/ey_kVcm/
ez_kVcm for the electrostatic deflectors) -- field-component columns are
whatever isn't x_cm/y_cm/z_cm.
"""
import csv
import sys


def load(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    fortran = load(sys.argv[1] if len(sys.argv) > 1 else "q1_fortran_reference.csv")
    cpp = load(sys.argv[2] if len(sys.argv) > 2 else "q1_cpp_reference.csv")

    if len(fortran) != len(cpp):
        print(f"FAIL: row count mismatch ({len(fortran)} vs {len(cpp)})")
        sys.exit(1)

    non_field_keys = ("x_cm", "y_cm", "z_cm", "name")
    field_keys = [k for k in fortran[0].keys() if k not in non_field_keys]

    max_abs_err = 0.0
    max_rel_err = 0.0
    worst = None
    n_bad = 0

    for i, (fr, cp) in enumerate(zip(fortran, cpp)):
        if fr.get("name") != cp.get("name"):
            print(f"FAIL: element-name mismatch at row {i}: {fr.get('name')} vs {cp.get('name')}")
            sys.exit(1)
        for key in ("x_cm", "y_cm", "z_cm"):
            if abs(float(fr[key]) - float(cp[key])) > 1e-9:
                print(f"FAIL: grid point mismatch at row {i}: {fr} vs {cp}")
                sys.exit(1)

        for key in field_keys:
            f_val = float(fr[key])
            c_val = float(cp[key])
            abs_err = abs(f_val - c_val)
            scale = max(abs(f_val), abs(c_val), 1e-12)
            rel_err = abs_err / scale

            if abs_err > 1e-9 and rel_err > 1e-6:
                n_bad += 1
                if abs_err > max_abs_err:
                    max_abs_err = abs_err
                    max_rel_err = rel_err
                    worst = (i, key, f_val, c_val)

    print(f"Compared {len(fortran)} grid points x {len(field_keys)} components ({', '.join(field_keys)})")
    print(f"Worst absolute error: {max_abs_err:.3e} (relative: {max_rel_err:.3e})")
    if worst:
        print(f"  at row {worst[0]}, {worst[1]}: fortran={worst[2]:.8e}  cpp={worst[3]:.8e}")

    if n_bad == 0:
        print("PASS: C++ port matches the real GEANT3 source")
    else:
        print(f"FAIL: {n_bad} components disagree beyond tolerance")
        sys.exit(1)


if __name__ == "__main__":
    main()
