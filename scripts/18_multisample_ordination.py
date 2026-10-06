#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Compute Jaccard/Bray-Curtis distances, PCoA and UPGMA."
)
parser.add_argument("--breadth", required=True)
parser.add_argument("--tpm", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--threshold", type=float, required=True)
args = parser.parse_args()

BREADTH_FILE = Path(args.breadth)
TPM_FILE = Path(args.tpm)
OUTDIR = Path(args.outdir)
PRESENCE_THRESHOLD = args.threshold

OUTDIR.mkdir(parents=True, exist_ok=True)


# =========================================================
# Read matrix
# =========================================================

def read_matrix(path):
    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")

        header = next(reader)
        header = [x.rstrip("\r") for x in header]

        samples = header[1:]

        data = {}

        for row in reader:
            row = [x.rstrip("\r") for x in row]

            if not row:
                continue

            votu = row[0]

            data[votu] = {
                sample: float(value)
                for sample, value in zip(samples, row[1:])
            }

    return samples, data


samples, breadth = read_matrix(BREADTH_FILE)
tpm_samples, tpm = read_matrix(TPM_FILE)

if set(samples) != set(tpm_samples):
    raise RuntimeError(
        "Sample names differ between breadth and TPM matrices"
    )

votus = sorted(breadth)

print(f"Samples : {len(samples)}")
print(f"vOTUs   : {len(votus)}")


# =========================================================
# Jaccard distance
# =========================================================

presence_sets = {}

for sample in samples:
    presence_sets[sample] = {
        votu
        for votu in votus
        if breadth[votu].get(sample, 0.0) >= PRESENCE_THRESHOLD
    }


def jaccard_distance(a, b):
    union = a | b

    if not union:
        return 0.0

    return 1.0 - len(a & b) / len(union)


jaccard = []

for s1 in samples:
    row = []

    for s2 in samples:
        row.append(
            jaccard_distance(
                presence_sets[s1],
                presence_sets[s2]
            )
        )

    jaccard.append(row)


# =========================================================
# sqrt(TPM) Bray-Curtis distance
# =========================================================

sqrt_tpm = {}

for sample in samples:
    sqrt_tpm[sample] = []

    for votu in votus:
        value = tpm.get(votu, {}).get(sample, 0.0)
        sqrt_tpm[sample].append(
            math.sqrt(max(value, 0.0))
        )


def bray_curtis(x, y):
    numerator = sum(
        abs(a - b)
        for a, b in zip(x, y)
    )

    denominator = sum(
        a + b
        for a, b in zip(x, y)
    )

    if denominator == 0:
        return 0.0

    return numerator / denominator


bray = []

for s1 in samples:
    row = []

    for s2 in samples:
        row.append(
            bray_curtis(
                sqrt_tpm[s1],
                sqrt_tpm[s2]
            )
        )

    bray.append(row)


# =========================================================
# Write distance matrix
# =========================================================

def write_distance_matrix(path, labels, matrix):

    with path.open("w", newline="") as f:
        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(["sample"] + labels)

        for label, row in zip(labels, matrix):
            writer.writerow(
                [label] +
                [f"{x:.8f}" for x in row]
            )


write_distance_matrix(
    OUTDIR / "jaccard_distance.tsv",
    samples,
    jaccard
)

write_distance_matrix(
    OUTDIR / "braycurtis_sqrt_tpm_distance.tsv",
    samples,
    bray
)


# =========================================================
# Symmetric Jacobi eigendecomposition
# Pure Python — sufficient for 7 x 7 matrices
# =========================================================

def jacobi_eigen(matrix, tol=1e-12, max_iter=10000):

    n = len(matrix)

    A = [
        list(row)
        for row in matrix
    ]

    V = [
        [
            1.0 if i == j else 0.0
            for j in range(n)
        ]
        for i in range(n)
    ]

    for _ in range(max_iter):

        max_value = 0.0
        p = 0
        q = 1

        for i in range(n):
            for j in range(i + 1, n):

                value = abs(A[i][j])

                if value > max_value:
                    max_value = value
                    p = i
                    q = j

        if max_value < tol:
            break

        app = A[p][p]
        aqq = A[q][q]
        apq = A[p][q]

        phi = 0.5 * math.atan2(
            2.0 * apq,
            aqq - app
        )

        c = math.cos(phi)
        s = math.sin(phi)

        for i in range(n):

            if i == p or i == q:
                continue

            aip = A[i][p]
            aiq = A[i][q]

            A[i][p] = c * aip - s * aiq
            A[p][i] = A[i][p]

            A[i][q] = s * aip + c * aiq
            A[q][i] = A[i][q]

        A[p][p] = (
            c*c*app
            - 2*s*c*apq
            + s*s*aqq
        )

        A[q][q] = (
            s*s*app
            + 2*s*c*apq
            + c*c*aqq
        )

        A[p][q] = 0.0
        A[q][p] = 0.0

        for i in range(n):

            vip = V[i][p]
            viq = V[i][q]

            V[i][p] = c * vip - s * viq
            V[i][q] = s * vip + c * viq

    eigenvalues = [
        A[i][i]
        for i in range(n)
    ]

    return eigenvalues, V


# =========================================================
# Classical PCoA
# =========================================================

def pcoa(distance_matrix):

    n = len(distance_matrix)

    D2 = [
        [
            distance_matrix[i][j] ** 2
            for j in range(n)
        ]
        for i in range(n)
    ]

    row_mean = [
        sum(row) / n
        for row in D2
    ]

    grand_mean = (
        sum(sum(row) for row in D2)
        / (n * n)
    )

    B = []

    for i in range(n):

        row = []

        for j in range(n):

            value = -0.5 * (
                D2[i][j]
                - row_mean[i]
                - row_mean[j]
                + grand_mean
            )

            row.append(value)

        B.append(row)

    eigenvalues, eigenvectors = jacobi_eigen(B)

    order = sorted(
        range(n),
        key=lambda k: eigenvalues[k],
        reverse=True
    )

    eigenvalues = [
        eigenvalues[k]
        for k in order
    ]

    eigenvectors = [
        [
            eigenvectors[i][k]
            for k in order
        ]
        for i in range(n)
    ]

    positive_sum = sum(
        x
        for x in eigenvalues
        if x > 0
    )

    coordinates = []

    for i in range(n):

        row = []

        for k in range(n):

            lam = eigenvalues[k]

            if lam > 0:
                value = (
                    eigenvectors[i][k]
                    * math.sqrt(lam)
                )
            else:
                value = 0.0

            row.append(value)

        coordinates.append(row)

    explained = [
        (
            100.0 * x / positive_sum
            if x > 0 and positive_sum > 0
            else 0.0
        )
        for x in eigenvalues
    ]

    return eigenvalues, explained, coordinates


jac_eig, jac_exp, jac_coord = pcoa(jaccard)
bray_eig, bray_exp, bray_coord = pcoa(bray)


# =========================================================
# Write PCoA
# =========================================================

def write_pcoa(path, labels, eig, explained, coord):

    with path.open("w", newline="") as f:

        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "sample",
            "PCoA1",
            "PCoA2"
        ])

        for sample, values in zip(labels, coord):

            writer.writerow([
                sample,
                f"{values[0]:.8f}",
                f"{values[1]:.8f}",
            ])

    return {
        "PCoA1": explained[0],
        "PCoA2": explained[1],
        "negative_eigenvalues": sum(
            1 for x in eig if x < -1e-10
        )
    }


jac_info = write_pcoa(
    OUTDIR / "jaccard_pcoa.tsv",
    samples,
    jac_eig,
    jac_exp,
    jac_coord
)

bray_info = write_pcoa(
    OUTDIR / "braycurtis_sqrt_tpm_pcoa.tsv",
    samples,
    bray_eig,
    bray_exp,
    bray_coord
)


# =========================================================
# Average-linkage / UPGMA clustering
# =========================================================

def upgma(labels, matrix):

    clusters = {}

    for i, label in enumerate(labels):
        clusters[i] = {
            "name": label,
            "size": 1,
            "height": 0.0,
            "newick": label,
        }

    distances = {}

    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            distances[frozenset((i, j))] = matrix[i][j]

    next_id = len(labels)
    merges = []

    while len(clusters) > 1:

        pair, distance = min(
            distances.items(),
            key=lambda x: x[1]
        )

        a, b = tuple(pair)

        A = clusters[a]
        B = clusters[b]

        new_height = distance / 2.0

        branch_a = max(
            new_height - A["height"],
            0.0
        )

        branch_b = max(
            new_height - B["height"],
            0.0
        )

        newick = (
            f"({A['newick']}:{branch_a:.6f},"
            f"{B['newick']}:{branch_b:.6f})"
        )

        new_size = A["size"] + B["size"]

        merges.append([
            A["name"],
            B["name"],
            distance,
            new_size
        ])

        other_ids = [
            k
            for k in clusters
            if k not in (a, b)
        ]

        new_distances = {}

        for k in other_ids:

            da = distances[
                frozenset((a, k))
            ]

            db = distances[
                frozenset((b, k))
            ]

            dnew = (
                da * A["size"]
                + db * B["size"]
            ) / new_size

            new_distances[k] = dnew

        # remove distances involving a or b
        distances = {
            key: value
            for key, value in distances.items()
            if a not in key and b not in key
        }

        del clusters[a]
        del clusters[b]

        clusters[next_id] = {
            "name": f"cluster_{next_id}",
            "size": new_size,
            "height": new_height,
            "newick": newick,
        }

        for k, value in new_distances.items():
            distances[
                frozenset((next_id, k))
            ] = value

        next_id += 1

    final_cluster = next(iter(clusters.values()))

    return final_cluster["newick"] + ";", merges


jac_tree, jac_merges = upgma(
    samples,
    jaccard
)

bray_tree, bray_merges = upgma(
    samples,
    bray
)


def write_tree(prefix, tree, merges):

    with (
        OUTDIR / f"{prefix}_upgma.newick"
    ).open("w") as f:
        f.write(tree + "\n")

    with (
        OUTDIR / f"{prefix}_upgma_merges.tsv"
    ).open("w", newline="") as f:

        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "cluster_1",
            "cluster_2",
            "distance",
            "new_cluster_size"
        ])

        for a, b, d, size in merges:
            writer.writerow([
                a,
                b,
                f"{d:.8f}",
                size
            ])


write_tree(
    "jaccard",
    jac_tree,
    jac_merges
)

write_tree(
    "braycurtis_sqrt_tpm",
    bray_tree,
    bray_merges
)


# =========================================================
# Summary
# =========================================================

summary = OUTDIR / "ordination_summary.tsv"

with summary.open("w", newline="") as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "analysis",
        "axis",
        "explained_positive_variance_pct"
    ])

    writer.writerow([
        "Jaccard",
        "PCoA1",
        f"{jac_info['PCoA1']:.4f}"
    ])

    writer.writerow([
        "Jaccard",
        "PCoA2",
        f"{jac_info['PCoA2']:.4f}"
    ])

    writer.writerow([
        "BrayCurtis_sqrtTPM",
        "PCoA1",
        f"{bray_info['PCoA1']:.4f}"
    ])

    writer.writerow([
        "BrayCurtis_sqrtTPM",
        "PCoA2",
        f"{bray_info['PCoA2']:.4f}"
    ])


print()
print("====================================")
print("Ordination + clustering complete")
print("====================================")

print(
    "Jaccard PCoA1/2 explained: "
    f"{jac_info['PCoA1']:.2f}% / "
    f"{jac_info['PCoA2']:.2f}%"
)

print(
    "Bray-Curtis PCoA1/2 explained: "
    f"{bray_info['PCoA1']:.2f}% / "
    f"{bray_info['PCoA2']:.2f}%"
)

print(
    "Jaccard negative eigenvalues:",
    jac_info["negative_eigenvalues"]
)

print(
    "Bray-Curtis negative eigenvalues:",
    bray_info["negative_eigenvalues"]
)

print()
print("Jaccard tree:")
print(jac_tree)

print()
print("Bray-Curtis tree:")
print(bray_tree)

print()
print("Output directory:")
print(OUTDIR)
