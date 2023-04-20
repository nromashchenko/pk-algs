#! /usr/bin/env python3


__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
import subprocess
from typing import List


TIME_BIN = "/usr/bin/time"
XPAS_ALGS_BIN = "cmake-build-release/bin/xpas_algs"
OUTPUT_DIR = "revision-fig8"

ALGS = [
    "bb",
    "dc",
    "dccw",
]
ALG_ARGMAP = {
    "bb": [1,0,0],
    "dc": [0,1,0],
    "dccw": [0,0,1],
}


class Dataset:
    def __init__(self, name: str, ar_file: str, ghost_file: str) -> None:
        self.name = name
        self.ar_file = ar_file
        self.ghost_file = ghost_file


def check_exists(filename: str) -> None:
    if not os.path.exists(filename):
        raise RuntimeError("File does not exist: " + filename)


def check_necessary_files(datasets: List[Dataset]) -> None:
    necessary_files = [TIME_BIN, XPAS_ALGS_BIN,]
    necessary_files.extend([d.ar_file for d in datasets])
    necessary_files.extend([d.ghost_file for d in datasets])

    for f in necessary_files:
        check_exists(f)


def run_dataset(dataset: Dataset) -> None:

    for alg in ALGS:
        output_csv = f"{OUTPUT_DIR}/{dataset.name}-{alg}.csv"
        time_csv = f"{OUTPUT_DIR}/{dataset.name}-{alg}.time"
        command = [
            TIME_BIN,
            "-v",
            XPAS_ALGS_BIN,
            dataset.ar_file,
            dataset.ghost_file,]
        command.extend(ALG_ARGMAP[alg])
        command.append(output_csv)
        command = [str(v) for v in command]

        command_str = " ".join(str(v) for v in command)
        print(f"\nRunning: {command_str}...\n")

        with open(time_csv, "w") as f:
            return_code = subprocess.call(command, stderr=f, stdout=None)

            if return_code != 0:
                raise RuntimeError("xpas-algs returned error: " + str(return_code))


def run(datasets: List[Dataset]) -> None:
    check_necessary_files(datasets)

    for dataset in datasets:
        run_dataset(dataset)



if __name__ == "__main__":
    datasets = [
        Dataset(
            "neotrop",
            "/ngs/rappas/neotrop/temp/AR/extended_align.phylip.raxml.ancestralProbs",
            "data/neotrop.ghost_ids.txt"),
        Dataset(
            "D155",
            "/ngs/rappas/D155/temp/AR/extended_align.phylip.raxml.ancestralProbs",
            "data/d155.ghost_ids.txt"
        )
    ]

    run(datasets)