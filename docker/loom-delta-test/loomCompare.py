#!/usr/bin/env python3

import sys
import argparse
import loompy
import numpy as np


def main():
    description = """This script compares two loom files and checks that they contain identical data up to a constant"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--truth-loom", dest="truth_loom_path", required=True, help="Path to truth loom file")
    parser.add_argument("--check-loom", dest="check_loom_path", required=True, help="Path to loom file to check")
    parser.add_argument("--delta-cutoff", dest="delta_cutoff", required=True, help="Max delta value allowed", type=int)
    args = parser.parse_args()
    truth_loom = loompy.connect(args.truth_loom_path)
    check_loom = loompy.connect(args.check_loom_path)

    truth_loom_array = truth_loom[:, :]
    check_loom_array = check_loom[:, :]

    delta = np.sum(np.absolute(np.subtract(truth_loom_array, check_loom_array)))

    if delta < args.delta_cutoff:
        print("Matrices are identical")
        sys.exit(0)
    else:
        print("Matrices are not identical")
        sys.exit(1)


if __name__ == "__main__":
    main()

