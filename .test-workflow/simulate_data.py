#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "msprime",
#     "numpy",
#     "pandas",
#     "pyslim",
#     "tskit",
# ]
# ///
import os
import subprocess
import tempfile
from itertools import combinations

import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit


# Compute pairwise distances on a unit torus between 2D points.
def torus_distance(points):
    n = points.shape[0]
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dx = abs(points[i, 0] - points[j, 0])
            dy = abs(points[i, 1] - points[j, 1])
            dx = min(dx, 1 - dx)
            dy = min(dy, 1 - dy)
            dist_matrix[i, j] = np.sqrt(dx**2 + dy**2)
    return dist_matrix


def simulate_data(
    seed: int,
    D: float,
    SD: float,
    SM: float,
    outvcf: str,
    outmap: str,
    outdist: str,
    contig_id: str,
    individual_names: list,
):
    # Simulate the population in a continuous space
    with tempfile.NamedTemporaryFile(suffix=".trees", delete=False) as temp_file:
        outpath = temp_file.name
        subprocess.run(
            [
                "slim",
                "-p",
                "-s",
                str(seed),
                "-d",
                f"NE={D}",
                "-d",
                f"SD={SD}",
                "-d",
                f"SM={SM}",
                "-d",
                f'OUTPATH="{outpath}"',
                "../simulations/constant_density.slim",
            ],
            check=True,
        )
    ts = tskit.load(outpath)
    os.remove(outpath)
    # Recapitation
    ts = pyslim.recapitate(ts, ancestral_Ne=D, recombination_rate=1e-8)
    n_samples = len(individual_names)
    sampled_inds = np.random.choice(
        np.arange(ts.num_individuals), size=n_samples, replace=False
    )
    nodes = np.array([ts.individual(i).nodes for i in sampled_inds]).flatten()
    ts = ts.simplify(samples=nodes)
    # Overlay neutral mutations
    mts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)
    # Write mutations to vcf.gz
    read_fd, write_fd = os.pipe()
    write_pipe = os.fdopen(write_fd, "w")
    with open(outvcf, "wb") as vcf_file:
        proc = subprocess.Popen(
            ["bcftools", "view", "-O", "z"], stdin=read_fd, stdout=vcf_file
        )
        mts.write_vcf(
            write_pipe, individual_names=individual_names, contig_id=contig_id
        )
        write_pipe.close()
        os.close(read_fd)
        proc.wait()

        if proc.returncode != 0:
            raise RuntimeError("bcftools failed with status:", proc.returncode)
    # Compute a toy genetic map in plink format
    with open(outmap, "w") as f:
        f.write(f"{contig_id}\trs\t0\t1\t\n{contig_id}\trs\t100\t100000000\t\n")
    # Compute distances between individuals in the torus
    points = mts.individual_locations
    dist_matrix = torus_distance(points)
    n = points.shape[0]
    rows = []
    for i, j in combinations(range(n), 2):
        rows.append((individual_names[i], individual_names[j], dist_matrix[i, j]))
    df = pd.DataFrame(rows, columns=["ID1", "ID2", "distance"])
    df.to_csv(outdist, index=False)


if __name__ == "__main__":
    D = 1000
    SD = 0.1
    SM = 0.01
    seed = 1000
    n_dip_indv = 20
    indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]
    simulate_data(
        seed=seed,
        D=D,
        SD=SD,
        SM=SM,
        outvcf=f"data_{seed}.vcf.gz",
        outmap=f"data_{seed}.map",
        outdist="pairwise_distances.csv",
        contig_id="chr1",
        individual_names=indv_names,
    )
