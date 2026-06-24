#!/usr/bin/env python3
"""
Find centromeric haplotype marker k-mers.

The preferred input map contains at least:
    k-mer    chr    pos    strand

If no map is supplied, the script can generate a core-CEN k-mer map from a
reference FASTA and core-centromere coordinates.

For each k-mer, this script counts hits in each core centromere, selects
k-mers that are strongly biased toward one core centromere, then collapses
overlapping selected target hits into non-redundant physical blocks. The raw
selected k-mer count and the collapsed block count are both useful: the former
captures sequence-marker richness, while the latter avoids counting a SNP's
overlapping 23-mer halo as many independent observations.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


CORE_PREFIX = "CEN"
PERICEN_PREFIX = "PERICEN"
CORE_REGION_NAMES = {"cen", "centromere", "core", "core_centromere"}
PERICEN_REGION_NAMES = {"pericen_left", "pericen_right"}
REVCOMP_TABLE = str.maketrans("ACGT", "TGCA")
VALID_DNA_RE = re.compile("^[ACGT]+$")


@dataclass
class CoreInterval:
    label: str
    chrom: str
    start: int
    end: int
    region_class: str = "CEN"
    region_name: str = "centromere"


@dataclass
class KmerStats:
    total_hits: int = 0
    core_hits: Counter = field(default_factory=Counter)
    target_positions: dict[str, list[int]] = field(default_factory=lambda: defaultdict(list))


@dataclass
class SelectedKmer:
    kmer: str
    assigned_cen: str
    target_hits: int
    total_core_hits: int
    other_core_hits: int
    max_other_core_hits: int
    target_core_fraction: float
    target_vs_max_other_plus1: float
    total_map_hits: int


class DisjointSet:
    def __init__(self, items: Iterable[str]):
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Identify k-mers biased to one core centromere and collapse "
            "overlapping target hits into candidate cenhap strength blocks."
        )
    )
    p.add_argument(
        "--map-tsv",
        default=None,
        help="TSV from map_kmers_ac.2.1.py. If omitted, --fasta is used to generate a core-CEN map.",
    )
    p.add_argument("--coords", required=True, help="Centromere coordinate TSV")
    p.add_argument(
        "--include-pericen",
        action="store_true",
        help=(
            "In chr/region/start/end coordinate files, include pericen_left and "
            "pericen_right rows as one logical PERICENn region per chromosome."
        ),
    )
    p.add_argument(
        "--fasta",
        default=None,
        help=(
            "Reference FASTA used for validation/reporting. Required when --map-tsv "
            "is omitted, in which case a core-CEN k-mer map is generated from this FASTA."
        ),
    )
    p.add_argument("--prefix", required=True, help="Output prefix")
    p.add_argument(
        "--kmer-size",
        type=int,
        default=0,
        help="K-mer size. Default: infer from the first map row, or 23 for FASTA-generated maps.",
    )
    p.add_argument(
        "--kmer-step",
        type=int,
        default=10,
        help="Step between k-mers when generating a map from FASTA. Default: 10.",
    )
    p.add_argument(
        "--min-cen-count",
        type=int,
        default=10,
        help="Minimum total core-CEN count when generating a map from FASTA. Default: 10.",
    )
    p.add_argument(
        "--max-outside-ratio",
        type=float,
        default=0.25,
        help=(
            "When generating a map from FASTA, reject k-mers whose outside-CEN count "
            "is greater than core_CEN_count times this ratio. Default: 0.25."
        ),
    )
    p.add_argument(
        "--no-canonical",
        action="store_true",
        help="When generating a map from FASTA, do not collapse reverse complements.",
    )
    p.add_argument(
        "--min-target-hits",
        type=int,
        default=3,
        help="Minimum hits for a k-mer inside its best core centromere.",
    )
    p.add_argument(
        "--max-other-core-hits",
        type=int,
        default=1,
        help="Maximum hits allowed across all other core centromeres.",
    )
    p.add_argument(
        "--min-target-core-fraction",
        type=float,
        default=0.90,
        help="Minimum fraction of core-centromere hits in the best centromere.",
    )
    p.add_argument(
        "--min-target-enrichment",
        type=float,
        default=3.0,
        help="Minimum target_hits / (max_other_core_hits + 1).",
    )
    p.add_argument(
        "--max-map-hits",
        type=int,
        default=0,
        help="Optional maximum total physical hits for a k-mer. 0 disables this filter.",
    )
    p.add_argument(
        "--merge-gap",
        type=int,
        default=0,
        help="Merge selected target intervals separated by no more than this many bases.",
    )
    p.add_argument(
        "--window-size",
        type=int,
        default=100_000,
        help="Running-window size for local cenhap strength.",
    )
    p.add_argument(
        "--window-step",
        type=int,
        default=25_000,
        help="Running-window step for local cenhap strength.",
    )
    p.add_argument(
        "--bin-size",
        type=int,
        default=200_000,
        help="Fixed-bin size for the simple regional cenhap k-mer count plot.",
    )
    p.add_argument(
        "--min-shared-cens",
        type=int,
        default=2,
        help="Minimum number of core CENs for a k-mer to enter shared-kmer diagnostics.",
    )
    p.add_argument(
        "--min-shared-hits-per-cen",
        type=int,
        default=1,
        help="Minimum hits in a CEN for that CEN to count as present in shared-kmer diagnostics.",
    )
    p.add_argument(
        "--write-window-plot",
        action="store_true",
        help="Also write the older overlapping-window local strength SVG.",
    )
    p.add_argument(
        "--write-all-kmers",
        action="store_true",
        help="Write all k-mers to the summary table. Default: write selected plus failed-threshold rows only.",
    )
    p.add_argument(
        "--skip-plot",
        action="store_true",
        help="Do not write SVG plots.",
    )
    return p.parse_args()


def clean_int(value: str) -> int:
    return int(str(value).replace("_", "").replace(",", ""))


def normalize_chrom(value: str) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        suffix = text[3:]
    else:
        suffix = text
    return f"Chr{suffix}"


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    suffix = normalize_chrom(chrom)[3:]
    if suffix.isdigit():
        return (0, int(suffix))
    return (1, suffix)


def label_sort_key(label: str) -> tuple[int, int | str, int, str]:
    if label.startswith(PERICEN_PREFIX):
        suffix = label[len(PERICEN_PREFIX) :]
        region_order = 1
    elif label.startswith(CORE_PREFIX):
        suffix = label[len(CORE_PREFIX) :]
        region_order = 0
    else:
        suffix = label
        region_order = 2
    chrom_key = (0, int(suffix)) if suffix.isdigit() else (1, suffix)
    return (chrom_key[0], chrom_key[1], region_order, label)


def unique_labels(intervals: list[CoreInterval]) -> list[str]:
    return sorted({interval.label for interval in intervals}, key=label_sort_key)


def intervals_by_label(intervals: list[CoreInterval]) -> dict[str, list[CoreInterval]]:
    grouped: dict[str, list[CoreInterval]] = defaultdict(list)
    for interval in intervals:
        grouped[interval.label].append(interval)
    for rows in grouped.values():
        rows.sort(key=lambda item: (chrom_sort_key(item.chrom), item.start, item.end))
    return dict(grouped)


def label_chrom(label: str, grouped: dict[str, list[CoreInterval]]) -> str:
    intervals = grouped.get(label, [])
    return intervals[0].chrom if intervals else ""


def label_extent(label: str, grouped: dict[str, list[CoreInterval]]) -> tuple[int | str, int | str]:
    intervals = grouped.get(label, [])
    if not intervals:
        return "", ""
    return min(interval.start for interval in intervals), max(interval.end for interval in intervals)


def label_length(label: str, grouped: dict[str, list[CoreInterval]]) -> int:
    return sum(interval.end - interval.start + 1 for interval in grouped.get(label, []))


def label_region_class(label: str) -> str:
    if label.startswith(PERICEN_PREFIX):
        return "PERICEN"
    return "CEN"


def revcomp(seq: str) -> str:
    return seq.translate(REVCOMP_TABLE)[::-1]


def canonical_kmer(seq: str, canonical: bool) -> tuple[str, str]:
    if not canonical:
        return seq, "+"
    rc = revcomp(seq)
    if rc < seq:
        return rc, "-"
    return seq, "+"


def open_text(path: str | Path):
    text = str(path)
    if text.endswith(".gz"):
        return gzip.open(text, "rt")
    return open(text, "rt")


def fasta_records(path: str | Path):
    name = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def read_fasta_sequences(path: str | Path) -> dict[str, str]:
    sequences = {}
    for name, seq in fasta_records(path):
        sequences[normalize_chrom(name)] = seq
    if not sequences:
        raise SystemExit(f"No FASTA records found in {path}")
    return sequences


def read_core_intervals(path: str, include_pericen: bool = False) -> list[CoreInterval]:
    intervals: list[CoreInterval] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"No header found in coordinate file: {path}")
        fields = {name.lower(): name for name in reader.fieldnames}

        if {"chr", "left", "right"}.issubset(fields):
            for row in reader:
                chrom = normalize_chrom(row[fields["chr"]])
                suffix = chrom[3:]
                intervals.append(
                    CoreInterval(
                        label=f"{CORE_PREFIX}{suffix}",
                        chrom=chrom,
                        start=clean_int(row[fields["left"]]),
                        end=clean_int(row[fields["right"]]),
                        region_class="CEN",
                        region_name="centromere",
                    )
                )
        elif {"chr", "region", "start", "end"}.issubset(fields):
            for row in reader:
                region = row[fields["region"]].strip().lower()
                chrom = normalize_chrom(row[fields["chr"]])
                suffix = chrom[3:]
                if region in CORE_REGION_NAMES:
                    label = f"{CORE_PREFIX}{suffix}"
                    region_class = "CEN"
                elif include_pericen and region in PERICEN_REGION_NAMES:
                    label = f"{PERICEN_PREFIX}{suffix}"
                    region_class = "PERICEN"
                else:
                    continue
                intervals.append(
                    CoreInterval(
                        label=label,
                        chrom=chrom,
                        start=clean_int(row[fields["start"]]),
                        end=clean_int(row[fields["end"]]),
                        region_class=region_class,
                        region_name=region,
                    )
                )
        else:
            raise SystemExit(
                "Coordinate file must have either Chr/Left/Right or chr/region/start/end columns."
            )

    if not intervals:
        if include_pericen:
            raise SystemExit("No accepted CEN or periCEN intervals were found.")
        raise SystemExit("No core centromere intervals were found.")
    intervals.sort(key=lambda item: (chrom_sort_key(item.chrom), item.start, item.end, item.label))
    return intervals


def read_fasta_lengths(path: str | None) -> dict[str, int]:
    if not path:
        return {}
    lengths: dict[str, int] = {}
    for name, seq in fasta_records(path):
        lengths[normalize_chrom(name)] = len(seq)
    return lengths


def validate_intervals(intervals: Iterable[CoreInterval], lengths: dict[str, int]) -> list[str]:
    warnings = []
    if not lengths:
        return warnings
    for interval in intervals:
        length = lengths.get(interval.chrom)
        if length is None:
            warnings.append(f"missing_fasta_chrom\t{interval.chrom}")
        elif interval.end > length:
            warnings.append(f"interval_beyond_fasta\t{interval.label}\t{interval.end}\t{length}")
    return warnings


def build_interval_index(intervals: list[CoreInterval]) -> dict[str, list[CoreInterval]]:
    by_chrom: dict[str, list[CoreInterval]] = defaultdict(list)
    for interval in intervals:
        by_chrom[interval.chrom].append(interval)
    return by_chrom


def find_interval(pos: int, intervals: list[CoreInterval]) -> CoreInterval | None:
    for interval in intervals:
        if interval.start <= pos <= interval.end:
            return interval
    return None


def next_step_aligned_start(start0: int, step: int) -> int:
    remainder = start0 % step
    return start0 if remainder == 0 else start0 + (step - remainder)


def iter_interval_windows(
    seq: str,
    interval: CoreInterval,
    kmer_size: int,
    step: int,
    canonical: bool,
):
    start0 = next_step_aligned_start(interval.start - 1, step)
    final_start0 = interval.end - kmer_size
    for pos0 in range(start0, final_start0 + 1, step):
        window = seq[pos0 : pos0 + kmer_size]
        if len(window) != kmer_size or not VALID_DNA_RE.match(window):
            continue
        kmer, strand = canonical_kmer(window, canonical)
        yield kmer, pos0 + 1, strand


def iter_genome_windows(seq: str, kmer_size: int, step: int, canonical: bool):
    for pos0 in range(0, len(seq) - kmer_size + 1, step):
        window = seq[pos0 : pos0 + kmer_size]
        if not VALID_DNA_RE.match(window):
            continue
        kmer, strand = canonical_kmer(window, canonical)
        yield kmer, pos0 + 1, strand


def generate_core_kmer_map(
    fasta: str,
    intervals: list[CoreInterval],
    prefix: Path,
    args: argparse.Namespace,
) -> tuple[Path, Path, dict[str, int]]:
    if args.kmer_size <= 0:
        args.kmer_size = 23
    if args.kmer_step <= 0:
        raise SystemExit("--kmer-step must be positive")
    if args.min_cen_count <= 0:
        raise SystemExit("--min-cen-count must be positive")
    if args.max_outside_ratio < 0:
        raise SystemExit("--max-outside-ratio must be non-negative")

    canonical = not args.no_canonical
    sequences = read_fasta_sequences(fasta)
    inside_counts: Counter = Counter()
    positions_by_kmer: dict[str, list[dict[str, object]]] = defaultdict(list)
    intervals_missing_sequence = 0
    core_windows_seen = 0
    core_windows_valid = 0

    for interval in intervals:
        seq = sequences.get(interval.chrom)
        if seq is None:
            intervals_missing_sequence += 1
            continue
        for kmer, pos, strand in iter_interval_windows(
            seq, interval, args.kmer_size, args.kmer_step, canonical
        ):
            core_windows_seen += 1
            core_windows_valid += 1
            inside_counts[kmer] += 1
            positions_by_kmer[kmer].append(
                {
                    "k-mer": kmer,
                    "chr": interval.chrom,
                    "pos": pos,
                    "strand": strand,
                    "source_cen": interval.label,
                }
            )

    candidate_set = set(inside_counts)
    total_counts: Counter = Counter()
    genome_windows_seen = 0
    genome_windows_matching_candidates = 0
    for chrom, seq in sequences.items():
        for kmer, _pos, _strand in iter_genome_windows(
            seq, args.kmer_size, args.kmer_step, canonical
        ):
            genome_windows_seen += 1
            if kmer in candidate_set:
                total_counts[kmer] += 1
                genome_windows_matching_candidates += 1

    selected = []
    rejected_low = 0
    rejected_outside = 0
    for kmer, inside_count in inside_counts.items():
        total_count = total_counts.get(kmer, 0)
        outside_count = max(0, total_count - inside_count)
        low = inside_count < args.min_cen_count
        outside_high = outside_count > inside_count * args.max_outside_ratio
        if low:
            rejected_low += 1
        if outside_high:
            rejected_outside += 1
        if not low and not outside_high:
            selected.append(kmer)
    selected.sort(key=lambda kmer: (-inside_counts[kmer], kmer))
    selected_set = set(selected)

    map_path = Path(str(prefix) + ".generated_core_kmer_map.tsv")
    stats_path = Path(str(prefix) + ".generated_core_kmer_map.stats.tsv")
    with open(map_path, "w", newline="") as out:
        fieldnames = [
            "k-mer",
            "chr",
            "pos",
            "strand",
            "source_cen",
            "core_cen_count",
            "genome_total_count",
            "outside_cen_count",
            "outside_to_core_ratio",
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for kmer in selected:
            inside_count = inside_counts[kmer]
            total_count = total_counts.get(kmer, 0)
            outside_count = max(0, total_count - inside_count)
            ratio = outside_count / inside_count if inside_count else 0.0
            for row in positions_by_kmer[kmer]:
                writer.writerow(
                    {
                        **row,
                        "core_cen_count": inside_count,
                        "genome_total_count": total_count,
                        "outside_cen_count": outside_count,
                        "outside_to_core_ratio": f"{ratio:.6f}",
                    }
                )

    generated_rows = sum(len(positions_by_kmer[kmer]) for kmer in selected)
    with open(stats_path, "w") as out:
        out.write(f"fasta\t{fasta}\n")
        out.write(f"coords\t{args.coords}\n")
        out.write(f"kmer_size\t{args.kmer_size}\n")
        out.write(f"kmer_step\t{args.kmer_step}\n")
        out.write(f"canonical\t{int(canonical)}\n")
        out.write(f"min_cen_count\t{args.min_cen_count}\n")
        out.write(f"max_outside_ratio\t{args.max_outside_ratio}\n")
        out.write(f"core_candidate_kmers\t{len(inside_counts)}\n")
        out.write(f"selected_generated_kmers\t{len(selected)}\n")
        out.write(f"generated_map_rows\t{generated_rows}\n")
        out.write(f"rejected_low_cen_count\t{rejected_low}\n")
        out.write(f"rejected_outside_gt_core_times_ratio\t{rejected_outside}\n")
        out.write(f"core_windows_seen\t{core_windows_seen}\n")
        out.write(f"core_windows_valid\t{core_windows_valid}\n")
        out.write(f"genome_windows_seen\t{genome_windows_seen}\n")
        out.write(f"genome_windows_matching_candidates\t{genome_windows_matching_candidates}\n")
        out.write(f"intervals_missing_sequence\t{intervals_missing_sequence}\n")
        out.write(f"generated_map_tsv\t{map_path}\n")

    return (
        map_path,
        stats_path,
        {
            "core_candidate_kmers": len(inside_counts),
            "selected_generated_kmers": len(selected),
            "generated_map_rows": generated_rows,
            "rejected_low_cen_count": rejected_low,
            "rejected_outside_gt_core_times_ratio": rejected_outside,
        },
    )


def infer_kmer_size(map_tsv: str) -> int:
    with open(map_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            return len(row["k-mer"])
    raise SystemExit(f"No map rows found in {map_tsv}")


def read_map_counts(
    map_tsv: str, interval_index: dict[str, list[CoreInterval]]
) -> tuple[dict[str, KmerStats], int, int]:
    stats: dict[str, KmerStats] = defaultdict(KmerStats)
    total_counts_from_map: dict[str, int] = {}
    rows = 0
    core_rows = 0
    with open(map_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"k-mer", "chr", "pos"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Map file is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            rows += 1
            kmer = row["k-mer"]
            chrom = normalize_chrom(row["chr"])
            pos = clean_int(row["pos"])
            record = stats[kmer]
            record.total_hits += 1
            if row.get("genome_total_count"):
                total_counts_from_map[kmer] = max(
                    total_counts_from_map.get(kmer, 0),
                    clean_int(row["genome_total_count"]),
                )
            interval = find_interval(pos, interval_index.get(chrom, []))
            if interval is not None:
                core_rows += 1
                record.core_hits[interval.label] += 1
                record.target_positions[interval.label].append(pos)
    for kmer, total_count in total_counts_from_map.items():
        stats[kmer].total_hits = max(stats[kmer].total_hits, total_count)
    return stats, rows, core_rows


def select_kmers(
    stats: dict[str, KmerStats], labels: list[str], args: argparse.Namespace
) -> tuple[dict[str, SelectedKmer], Counter]:
    selected: dict[str, SelectedKmer] = {}
    fail_counts: Counter = Counter()

    for kmer, record in stats.items():
        total_core_hits = sum(record.core_hits.values())
        ranked = sorted(
            ((label, record.core_hits.get(label, 0)) for label in labels),
            key=lambda item: (-item[1], item[0]),
        )
        assigned_cen, target_hits = ranked[0]
        other_counts = [count for label, count in ranked[1:]]
        other_core_hits = total_core_hits - target_hits
        max_other = max(other_counts) if other_counts else 0
        target_fraction = target_hits / total_core_hits if total_core_hits else 0.0
        enrichment = target_hits / (max_other + 1)
        reasons = fail_reasons(
            record.total_hits,
            target_hits,
            other_core_hits,
            target_fraction,
            enrichment,
            args,
        )
        if reasons:
            fail_counts.update(reasons)
            continue
        selected[kmer] = SelectedKmer(
            kmer=kmer,
            assigned_cen=assigned_cen,
            target_hits=target_hits,
            total_core_hits=total_core_hits,
            other_core_hits=other_core_hits,
            max_other_core_hits=max_other,
            target_core_fraction=target_fraction,
            target_vs_max_other_plus1=enrichment,
            total_map_hits=record.total_hits,
        )
    return selected, fail_counts


def fail_reasons(
    total_map_hits: int,
    target_hits: int,
    other_core_hits: int,
    target_fraction: float,
    enrichment: float,
    args: argparse.Namespace,
) -> list[str]:
    reasons = []
    if target_hits < args.min_target_hits:
        reasons.append("low_target_hits")
    if other_core_hits > args.max_other_core_hits:
        reasons.append("other_core_hits")
    if target_fraction < args.min_target_core_fraction:
        reasons.append("low_target_core_fraction")
    if enrichment < args.min_target_enrichment:
        reasons.append("low_target_enrichment")
    if args.max_map_hits and total_map_hits > args.max_map_hits:
        reasons.append("high_total_map_hits")
    return reasons


def write_summary(
    path: Path,
    stats: dict[str, KmerStats],
    selected: dict[str, SelectedKmer],
    labels: list[str],
    args: argparse.Namespace,
) -> None:
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "kmer",
                "selected",
                "fail_reasons",
                "assigned_cen",
                "total_map_hits",
                "total_core_hits",
                "target_core_hits",
                "other_core_hits",
                "max_other_core_hits",
                "target_core_fraction",
                "target_vs_max_other_plus1",
                *[f"{label}_hits" for label in labels],
            ]
        )
        for kmer in sorted(stats):
            record = stats[kmer]
            total_core_hits = sum(record.core_hits.values())
            ranked = sorted(
                ((label, record.core_hits.get(label, 0)) for label in labels),
                key=lambda item: (-item[1], item[0]),
            )
            assigned_cen, target_hits = ranked[0]
            other_core_hits = total_core_hits - target_hits
            max_other = max((count for _, count in ranked[1:]), default=0)
            target_fraction = target_hits / total_core_hits if total_core_hits else 0.0
            enrichment = target_hits / (max_other + 1)
            reasons = fail_reasons(
                record.total_hits,
                target_hits,
                other_core_hits,
                target_fraction,
                enrichment,
                args,
            )
            if not args.write_all_kmers and kmer not in selected and "low_target_hits" in reasons:
                continue
            writer.writerow(
                [
                    kmer,
                    int(kmer in selected),
                    ",".join(reasons),
                    assigned_cen,
                    record.total_hits,
                    total_core_hits,
                    target_hits,
                    other_core_hits,
                    max_other,
                    f"{target_fraction:.6f}",
                    f"{enrichment:.6f}",
                    *[record.core_hits.get(label, 0) for label in labels],
                ]
            )


def idf_weight(total_cens: int, present_cens: int) -> float:
    return math.log((total_cens + 1) / (present_cens + 1)) + 1.0


def split_core_pericen_labels(labels: list[str]) -> tuple[list[str], list[str]]:
    core_labels = [label for label in labels if label.startswith(CORE_PREFIX)]
    pericen_labels = [label for label in labels if label.startswith(PERICEN_PREFIX)]
    return core_labels, pericen_labels


def classify_core_pericen_presence(
    present_core: list[str],
    present_pericen: list[str],
    total_core: int,
    total_pericen: int,
) -> str:
    has_all_core = bool(total_core) and len(present_core) == total_core
    has_all_pericen = bool(total_pericen) and len(present_pericen) == total_pericen
    has_core = bool(present_core)
    has_pericen = bool(present_pericen)
    if has_all_core and has_all_pericen:
        return "all_core_cens_and_all_pericens"
    if has_all_core and not has_pericen:
        return "all_core_cens_only"
    if has_all_pericen and not has_core:
        return "all_pericens_only"
    if has_all_core:
        return "all_core_cens_plus_pericen_subset"
    if has_all_pericen:
        return "all_pericens_plus_core_cen_subset"
    if has_core and has_pericen:
        return "mixed_core_cen_pericen_subset"
    if has_core:
        return "core_cen_subset"
    if has_pericen:
        return "pericen_subset"
    return "none"


def build_shared_kmer_rows(
    stats: dict[str, KmerStats],
    labels: list[str],
    selected: dict[str, SelectedKmer],
    min_shared_cens: int,
    min_hits_per_cen: int,
) -> list[dict[str, object]]:
    if min_shared_cens <= 0:
        raise SystemExit("--min-shared-cens must be a positive integer.")
    if min_hits_per_cen <= 0:
        raise SystemExit("--min-shared-hits-per-cen must be a positive integer.")

    total_cens = len(labels)
    all_shared_class = (
        "all_analysis_regions" if any(label.startswith(PERICEN_PREFIX) for label in labels) else "all_core_cens"
    )
    subset_shared_class = (
        "subset_analysis_regions"
        if any(label.startswith(PERICEN_PREFIX) for label in labels)
        else "subset_core_cens"
    )
    core_labels, pericen_labels = split_core_pericen_labels(labels)
    rows = []
    for kmer, record in stats.items():
        counts = {label: int(record.core_hits.get(label, 0)) for label in labels}
        present = [label for label in labels if counts[label] >= min_hits_per_cen]
        if len(present) < min_shared_cens:
            continue
        present_core = [label for label in core_labels if counts[label] >= min_hits_per_cen]
        present_pericen = [label for label in pericen_labels if counts[label] >= min_hits_per_cen]
        present_counts = [counts[label] for label in present]
        selected_call = selected.get(kmer)
        rows.append(
            {
                "kmer": kmer,
                "shared_class": all_shared_class if len(present) == total_cens else subset_shared_class,
                "present_cens_count": len(present),
                "present_cens_fraction": len(present) / total_cens if total_cens else 0.0,
                "present_cens": ",".join(present),
                "present_region_class": classify_core_pericen_presence(
                    present_core,
                    present_pericen,
                    len(core_labels),
                    len(pericen_labels),
                ),
                "present_core_cens_count": len(present_core),
                "present_core_cens_fraction": (
                    len(present_core) / len(core_labels) if core_labels else 0.0
                ),
                "present_all_core_cens": int(bool(core_labels) and len(present_core) == len(core_labels)),
                "present_core_cens": ",".join(present_core),
                "present_pericens_count": len(present_pericen),
                "present_pericens_fraction": (
                    len(present_pericen) / len(pericen_labels) if pericen_labels else 0.0
                ),
                "present_all_pericens": int(
                    bool(pericen_labels) and len(present_pericen) == len(pericen_labels)
                ),
                "present_pericens": ",".join(present_pericen),
                "total_core_hits": sum(counts.values()),
                "total_core_cen_hits": sum(counts[label] for label in core_labels),
                "total_pericen_hits": sum(counts[label] for label in pericen_labels),
                "min_present_cen_hits": min(present_counts) if present_counts else 0,
                "max_present_cen_hits": max(present_counts) if present_counts else 0,
                "mean_present_cen_hits": (
                    sum(present_counts) / len(present_counts) if present_counts else 0.0
                ),
                "idf_weight": idf_weight(total_cens, len(present)),
                "selected_cenhap_defining": int(selected_call is not None),
                "selected_assigned_cen": selected_call.assigned_cen if selected_call else "",
                **{f"{label}_hits": counts[label] for label in labels},
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["present_cens_count"]),
            str(row["present_cens"]),
            -int(row["total_core_hits"]),
            str(row["kmer"]),
        )
    )
    return rows


def build_shared_kmer_set_rows(shared_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in shared_rows:
        grouped[str(row["present_cens"])].append(row)

    rows = []
    for cen_set, members in grouped.items():
        total_hits = sum(int(row["total_core_hits"]) for row in members)
        selected_members = sum(int(row["selected_cenhap_defining"]) for row in members)
        idf_values = [float(row["idf_weight"]) for row in members]
        rows.append(
            {
                "present_cens": cen_set,
                "present_cens_count": int(members[0]["present_cens_count"]) if members else 0,
                "present_region_class": str(members[0]["present_region_class"]) if members else "",
                "present_core_cens_count": (
                    int(members[0]["present_core_cens_count"]) if members else 0
                ),
                "present_core_cens_fraction": (
                    float(members[0]["present_core_cens_fraction"]) if members else 0.0
                ),
                "present_all_core_cens": int(members[0]["present_all_core_cens"]) if members else 0,
                "present_core_cens": str(members[0]["present_core_cens"]) if members else "",
                "present_pericens_count": (
                    int(members[0]["present_pericens_count"]) if members else 0
                ),
                "present_pericens_fraction": (
                    float(members[0]["present_pericens_fraction"]) if members else 0.0
                ),
                "present_all_pericens": int(members[0]["present_all_pericens"]) if members else 0,
                "present_pericens": str(members[0]["present_pericens"]) if members else "",
                "shared_class": str(members[0]["shared_class"]) if members else "",
                "kmers": len(members),
                "selected_cenhap_defining_kmers": selected_members,
                "fraction_selected_cenhap_defining": (
                    selected_members / len(members) if members else 0.0
                ),
                "total_core_hits": total_hits,
                "total_core_cen_hits": sum(int(row["total_core_cen_hits"]) for row in members),
                "total_pericen_hits": sum(int(row["total_pericen_hits"]) for row in members),
                "mean_total_core_hits_per_kmer": total_hits / len(members) if members else 0.0,
                "mean_idf_weight": sum(idf_values) / len(idf_values) if idf_values else 0.0,
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["present_cens_count"]),
            -int(row["kmers"]),
            str(row["present_cens"]),
        )
    )
    return rows


def build_idf_relatedness_rows(
    stats: dict[str, KmerStats],
    labels: list[str],
) -> tuple[list[dict[str, object]], list[list[object]]]:
    weights = {}
    for kmer, record in stats.items():
        present = sum(1 for label in labels if int(record.core_hits.get(label, 0)) > 0)
        if present:
            weights[kmer] = idf_weight(len(labels), present)

    pair_rows = []
    matrix = [["assigned_cen", *labels]]
    for left in labels:
        matrix_row = [left]
        for right in labels:
            numerator = 0.0
            denominator = 0.0
            shared_kmers = 0
            left_kmers = 0
            right_kmers = 0
            for kmer, record in stats.items():
                left_count = int(record.core_hits.get(left, 0))
                right_count = int(record.core_hits.get(right, 0))
                if left_count:
                    left_kmers += 1
                if right_count:
                    right_kmers += 1
                if left_count and right_count:
                    shared_kmers += 1
                if left_count or right_count:
                    weight = weights.get(kmer, 1.0)
                    numerator += weight * min(left_count, right_count)
                    denominator += weight * max(left_count, right_count)
            score = numerator / denominator if denominator else 0.0
            matrix_row.append(score)
            if left < right:
                pair_rows.append(
                    {
                        "cen_a": left,
                        "cen_b": right,
                        "idf_weighted_jaccard": score,
                        "shared_kmers": shared_kmers,
                        "cen_a_kmers": left_kmers,
                        "cen_b_kmers": right_kmers,
                        "binary_jaccard": (
                            shared_kmers / (left_kmers + right_kmers - shared_kmers)
                            if left_kmers + right_kmers - shared_kmers
                            else 0.0
                        ),
                    }
                )
        matrix.append(matrix_row)
    return pair_rows, matrix


def write_shared_kmer_rows(
    path: Path,
    rows: list[dict[str, object]],
    labels: list[str],
) -> None:
    fieldnames = [
        "kmer",
        "shared_class",
        "present_cens_count",
        "present_cens_fraction",
        "present_cens",
        "present_region_class",
        "present_core_cens_count",
        "present_core_cens_fraction",
        "present_all_core_cens",
        "present_core_cens",
        "present_pericens_count",
        "present_pericens_fraction",
        "present_all_pericens",
        "present_pericens",
        "total_core_hits",
        "total_core_cen_hits",
        "total_pericen_hits",
        "min_present_cen_hits",
        "max_present_cen_hits",
        "mean_present_cen_hits",
        "idf_weight",
        "selected_cenhap_defining",
        "selected_assigned_cen",
        *[f"{label}_hits" for label in labels],
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            formatted["present_cens_fraction"] = f"{float(row['present_cens_fraction']):.6f}"
            formatted["present_core_cens_fraction"] = (
                f"{float(row['present_core_cens_fraction']):.6f}"
            )
            formatted["present_pericens_fraction"] = (
                f"{float(row['present_pericens_fraction']):.6f}"
            )
            formatted["mean_present_cen_hits"] = f"{float(row['mean_present_cen_hits']):.6f}"
            formatted["idf_weight"] = f"{float(row['idf_weight']):.6f}"
            writer.writerow(formatted)


def write_shared_kmer_set_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "present_cens",
        "present_cens_count",
        "present_region_class",
        "present_core_cens_count",
        "present_core_cens_fraction",
        "present_all_core_cens",
        "present_core_cens",
        "present_pericens_count",
        "present_pericens_fraction",
        "present_all_pericens",
        "present_pericens",
        "shared_class",
        "kmers",
        "selected_cenhap_defining_kmers",
        "fraction_selected_cenhap_defining",
        "total_core_hits",
        "total_core_cen_hits",
        "total_pericen_hits",
        "mean_total_core_hits_per_kmer",
        "mean_idf_weight",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            for key in [
                "present_core_cens_fraction",
                "present_pericens_fraction",
                "fraction_selected_cenhap_defining",
                "mean_total_core_hits_per_kmer",
                "mean_idf_weight",
            ]:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def write_relatedness_pairs(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "cen_a",
        "cen_b",
        "idf_weighted_jaccard",
        "shared_kmers",
        "cen_a_kmers",
        "cen_b_kmers",
        "binary_jaccard",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            formatted["idf_weighted_jaccard"] = f"{float(row['idf_weighted_jaccard']):.6f}"
            formatted["binary_jaccard"] = f"{float(row['binary_jaccard']):.6f}"
            writer.writerow(formatted)


def write_relatedness_matrix(path: Path, matrix: list[list[object]]) -> None:
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        for row_idx, row in enumerate(matrix):
            if row_idx == 0:
                writer.writerow(row)
            else:
                writer.writerow([row[0], *[f"{float(value):.6f}" for value in row[1:]]])


def build_blocks(
    stats: dict[str, KmerStats],
    selected: dict[str, SelectedKmer],
    interval_groups: dict[str, list[CoreInterval]],
    kmer_size: int,
    merge_gap: int,
) -> dict[str, list[dict[str, object]]]:
    by_cen: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for kmer, call in selected.items():
        for pos in stats[kmer].target_positions.get(call.assigned_cen, []):
            by_cen[call.assigned_cen].append((pos, pos + kmer_size - 1, kmer))

    blocks_by_cen: dict[str, list[dict[str, object]]] = {}
    for label, intervals in by_cen.items():
        intervals.sort(key=lambda item: (item[0], item[1], item[2]))
        blocks: list[dict[str, object]] = []
        current_start = None
        current_end = None
        current_kmers: set[str] = set()
        current_hits = 0

        def flush() -> None:
            if current_start is None or current_end is None:
                return
            blocks.append(
                {
                    "assigned_cen": label,
                    "chr": label_chrom(label, interval_groups),
                    "block_start": current_start,
                    "block_end": current_end,
                    "block_length_bp": current_end - current_start + 1,
                    "distinct_kmers": len(current_kmers),
                    "target_map_hits": current_hits,
                    "representative_kmer": sorted(current_kmers)[0] if current_kmers else "",
                    "_kmers": set(current_kmers),
                }
            )

        for start, end, kmer in intervals:
            if current_start is None:
                current_start = start
                current_end = end
                current_kmers = {kmer}
                current_hits = 1
                continue
            if start <= current_end + merge_gap + 1:
                current_end = max(current_end, end)
                current_kmers.add(kmer)
                current_hits += 1
            else:
                flush()
                current_start = start
                current_end = end
                current_kmers = {kmer}
                current_hits = 1
        flush()
        blocks_by_cen[label] = blocks
    return blocks_by_cen


def write_blocks(path: Path, blocks_by_cen: dict[str, list[dict[str, object]]]) -> None:
    fieldnames = [
        "assigned_cen",
        "unit_ids",
        "chr",
        "block_start",
        "block_end",
        "block_length_bp",
        "distinct_kmers",
        "target_map_hits",
        "representative_kmer",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for label in sorted(blocks_by_cen):
            for block in blocks_by_cen[label]:
                row = {key: block[key] for key in fieldnames if key in block}
                row["unit_ids"] = ",".join(sorted(block.get("_unit_ids", set())))
                writer.writerow(row)


def build_units(
    labels: list[str],
    selected: dict[str, SelectedKmer],
    blocks_by_cen: dict[str, list[dict[str, object]]],
) -> dict[str, list[dict[str, object]]]:
    kmers_by_cen: dict[str, list[str]] = defaultdict(list)
    for kmer, call in selected.items():
        kmers_by_cen[call.assigned_cen].append(kmer)

    units_by_cen: dict[str, list[dict[str, object]]] = {}
    for label in labels:
        kmers = kmers_by_cen.get(label, [])
        dsu = DisjointSet(kmers)
        for block in blocks_by_cen.get(label, []):
            block_kmers = sorted(block.get("_kmers", set()))
            if not block_kmers:
                continue
            first = block_kmers[0]
            for kmer in block_kmers[1:]:
                dsu.union(first, kmer)

        grouped: dict[str, list[str]] = defaultdict(list)
        for kmer in kmers:
            grouped[dsu.find(kmer)].append(kmer)

        unit_rows = []
        for unit_kmers in grouped.values():
            unit_set = set(unit_kmers)
            unit_blocks = [
                block
                for block in blocks_by_cen.get(label, [])
                if unit_set.intersection(block.get("_kmers", set()))
            ]
            physical_blocks = len(unit_blocks)
            target_hits = sum(selected[kmer].target_hits for kmer in unit_kmers)
            span_start = min((int(block["block_start"]) for block in unit_blocks), default="")
            span_end = max((int(block["block_end"]) for block in unit_blocks), default="")
            span_bp = span_end - span_start + 1 if span_start != "" and span_end != "" else ""
            unit_rows.append(
                {
                    "assigned_cen": label,
                    "unit_id": "",
                    "chr": unit_blocks[0]["chr"] if unit_blocks else "",
                    "unit_span_start": span_start,
                    "unit_span_end": span_end,
                    "unit_span_bp": span_bp,
                    "distinct_kmers": len(unit_kmers),
                    "target_map_hits": target_hits,
                    "physical_blocks": physical_blocks,
                    "representative_kmer": sorted(unit_kmers)[0],
                    "_kmers": unit_set,
                }
            )
        unit_rows.sort(
            key=lambda row: (
                -int(row["distinct_kmers"]),
                -int(row["target_map_hits"]),
                str(row["representative_kmer"]),
            )
        )
        for idx, row in enumerate(unit_rows, start=1):
            row["unit_id"] = f"{label}_unit_{idx:06d}"
        unit_by_kmer = {
            kmer: str(row["unit_id"])
            for row in unit_rows
            for kmer in row.get("_kmers", set())
        }
        for block in blocks_by_cen.get(label, []):
            block["_unit_ids"] = {
                unit_by_kmer[kmer] for kmer in block.get("_kmers", set()) if kmer in unit_by_kmer
            }
        units_by_cen[label] = unit_rows
    return units_by_cen


def write_units(path: Path, units_by_cen: dict[str, list[dict[str, object]]]) -> None:
    fieldnames = [
        "assigned_cen",
        "unit_id",
        "chr",
        "unit_span_start",
        "unit_span_end",
        "unit_span_bp",
        "distinct_kmers",
        "target_map_hits",
        "physical_blocks",
        "representative_kmer",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for label in sorted(units_by_cen):
            for row in units_by_cen[label]:
                writer.writerow({key: row[key] for key in fieldnames})


def write_units_bed(path: Path, units_by_cen: dict[str, list[dict[str, object]]]) -> None:
    max_kmers = max(
        (int(row["distinct_kmers"]) for rows in units_by_cen.values() for row in rows),
        default=1,
    )
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for label in sorted(units_by_cen):
            for row in units_by_cen[label]:
                if row["unit_span_start"] == "" or row["unit_span_end"] == "":
                    continue
                start0 = int(row["unit_span_start"]) - 1
                end = int(row["unit_span_end"])
                score = round((int(row["distinct_kmers"]) / max_kmers) * 1000)
                writer.writerow(
                    [
                        row["chr"],
                        start0,
                        end,
                        row["unit_id"],
                        score,
                        ".",
                        int(row["unit_span_start"]) - 1,
                        end,
                        "47,111,115",
                    ]
                )


def quantile(values: list[int], q: float) -> float:
    if not values:
        return 0.0
    ranked = sorted(values)
    if len(ranked) == 1:
        return float(ranked[0])
    pos = (len(ranked) - 1) * q
    lower = int(pos)
    upper = min(lower + 1, len(ranked) - 1)
    fraction = pos - lower
    return ranked[lower] * (1 - fraction) + ranked[upper] * fraction


def unit_size_bucket(size: int) -> str:
    if size <= 1:
        return "1"
    if size == 2:
        return "2"
    if size <= 5:
        return "3-5"
    if size <= 10:
        return "6-10"
    if size <= 20:
        return "11-20"
    if size <= 50:
        return "21-50"
    if size <= 100:
        return "51-100"
    if size <= 500:
        return "101-500"
    return "501+"


def unit_size_bucket_order() -> list[str]:
    return ["1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", "101-500", "501+"]


def build_unit_size_summary(
    labels: list[str],
    units_by_cen: dict[str, list[dict[str, object]]],
) -> list[dict[str, object]]:
    rows = []
    for label in labels:
        units = units_by_cen.get(label, [])
        sizes = [int(unit["distinct_kmers"]) for unit in units]
        hits = [int(unit["target_map_hits"]) for unit in units]
        blocks = [int(unit["physical_blocks"]) for unit in units]
        total_size = sum(sizes)
        largest = max(sizes) if sizes else 0
        rows.append(
            {
                "assigned_cen": label,
                "units": len(units),
                "singleton_units": sum(1 for size in sizes if size == 1),
                "fraction_singleton_units": (
                    sum(1 for size in sizes if size == 1) / len(units) if units else 0.0
                ),
                "total_distinct_kmers_in_units": total_size,
                "min_distinct_kmers_per_unit": min(sizes) if sizes else 0,
                "q25_distinct_kmers_per_unit": quantile(sizes, 0.25),
                "median_distinct_kmers_per_unit": quantile(sizes, 0.50),
                "q75_distinct_kmers_per_unit": quantile(sizes, 0.75),
                "max_distinct_kmers_per_unit": largest,
                "mean_distinct_kmers_per_unit": total_size / len(units) if units else 0.0,
                "largest_unit_fraction_of_kmers": largest / total_size if total_size else 0.0,
                "mean_target_map_hits_per_unit": sum(hits) / len(units) if units else 0.0,
                "mean_physical_blocks_per_unit": sum(blocks) / len(units) if units else 0.0,
            }
        )
    return rows


def build_unit_size_distribution(
    labels: list[str],
    units_by_cen: dict[str, list[dict[str, object]]],
) -> list[dict[str, object]]:
    bucket_order = unit_size_bucket_order()
    rows = []
    for label in labels:
        units = units_by_cen.get(label, [])
        counts = Counter(unit_size_bucket(int(unit["distinct_kmers"])) for unit in units)
        total_units = len(units)
        for bucket in bucket_order:
            count = counts.get(bucket, 0)
            rows.append(
                {
                    "assigned_cen": label,
                    "size_bucket_distinct_kmers": bucket,
                    "units": count,
                    "fraction_units": count / total_units if total_units else 0.0,
                }
            )
    return rows


def write_unit_size_summary(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "assigned_cen",
        "units",
        "singleton_units",
        "fraction_singleton_units",
        "total_distinct_kmers_in_units",
        "min_distinct_kmers_per_unit",
        "q25_distinct_kmers_per_unit",
        "median_distinct_kmers_per_unit",
        "q75_distinct_kmers_per_unit",
        "max_distinct_kmers_per_unit",
        "mean_distinct_kmers_per_unit",
        "largest_unit_fraction_of_kmers",
        "mean_target_map_hits_per_unit",
        "mean_physical_blocks_per_unit",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            for key in [
                "fraction_singleton_units",
                "q25_distinct_kmers_per_unit",
                "median_distinct_kmers_per_unit",
                "q75_distinct_kmers_per_unit",
                "mean_distinct_kmers_per_unit",
                "largest_unit_fraction_of_kmers",
                "mean_target_map_hits_per_unit",
                "mean_physical_blocks_per_unit",
            ]:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def write_unit_size_distribution(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "assigned_cen",
        "size_bucket_distinct_kmers",
        "units",
        "fraction_units",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            formatted["fraction_units"] = f"{float(formatted['fraction_units']):.6f}"
            writer.writerow(formatted)


def write_unit_size_plot(path: Path, rows: list[dict[str, object]]) -> str:
    labels = sorted({str(row["assigned_cen"]) for row in rows})
    buckets = unit_size_bucket_order()
    by_key = {
        (str(row["assigned_cen"]), str(row["size_bucket_distinct_kmers"])): int(row["units"])
        for row in rows
    }
    width = 1320
    height = 530
    margin_left = 78
    margin_right = 36
    margin_top = 74
    margin_bottom = 72
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    max_count = max((int(row["units"]) for row in rows), default=1)
    colors = ["#2f6f73", "#9b5d2e", "#5f6f95", "#8a6d3b", "#4f7f4f"]

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f1f1f}",
        ".title{font-size:22px;font-weight:700}",
        ".tick{font-size:11px;fill:#555}",
        ".legend{font-size:12px}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#d8d8d8;stroke-width:1}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="660" y="34" text-anchor="middle">'
        "CU Size Distribution</text>",
        '<text class="tick" x="660" y="56" text-anchor="middle">'
        "Size is distinct selected k-mers per CU</text>",
    ]

    baseline = margin_top + plot_height
    for idx in range(5):
        tick = round(max_count * idx / 4)
        y = baseline - (tick / max_count) * plot_height if max_count else baseline
        svg.append(
            f'<line class="grid" x1="{margin_left}" y1="{y:.1f}" '
            f'x2="{margin_left + plot_width}" y2="{y:.1f}"/>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left - 8}" y="{y + 4:.1f}" '
            f'text-anchor="end">{tick:,}</text>'
        )

    svg.append(
        f'<line class="axis" x1="{margin_left}" y1="{baseline}" '
        f'x2="{margin_left + plot_width}" y2="{baseline}"/>'
    )
    svg.append(
        f'<line class="axis" x1="{margin_left}" y1="{margin_top}" '
        f'x2="{margin_left}" y2="{baseline}"/>'
    )

    bucket_gap = 22
    group_width = (plot_width - bucket_gap * (len(buckets) - 1)) / len(buckets)
    bar_gap = 3
    bar_width = (group_width - bar_gap * (len(labels) - 1)) / max(1, len(labels))
    for bucket_idx, bucket in enumerate(buckets):
        group_x = margin_left + bucket_idx * (group_width + bucket_gap)
        for label_idx, label in enumerate(labels):
            count = by_key.get((label, bucket), 0)
            bar_h = (count / max_count) * plot_height if max_count else 0
            x = group_x + label_idx * (bar_width + bar_gap)
            y = baseline - bar_h
            color = colors[label_idx % len(colors)]
            svg.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" '
                f'height="{bar_h:.1f}" fill="{color}" opacity="0.86">'
                f'<title>{html.escape(label)} {html.escape(bucket)}: {count:,} CU</title>'
                "</rect>"
            )
        svg.append(
            f'<text class="tick" x="{group_x + group_width / 2:.1f}" y="{baseline + 20}" '
            f'text-anchor="middle">{html.escape(bucket)}</text>'
        )

    legend_x = margin_left
    legend_y = height - 20
    for label_idx, label in enumerate(labels):
        x = legend_x + label_idx * 96
        color = colors[label_idx % len(colors)]
        svg.append(f'<rect x="{x}" y="{legend_y - 12}" width="12" height="12" fill="{color}"/>')
        svg.append(
            f'<text class="legend" x="{x + 18}" y="{legend_y - 2}">{html.escape(label)}</text>'
        )

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")
    return ""


def write_cen_strength(
    path: Path,
    labels: list[str],
    selected: dict[str, SelectedKmer],
    blocks_by_cen: dict[str, list[dict[str, object]]],
    units_by_cen: dict[str, list[dict[str, object]]],
) -> None:
    selected_by_cen = Counter(call.assigned_cen for call in selected.values())
    target_hits_by_cen = Counter()
    for call in selected.values():
        target_hits_by_cen[call.assigned_cen] += call.target_hits

    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "assigned_cen",
                "selected_distinct_kmers",
                "target_map_hits",
                "cenhap_strength_units",
                "cenhap_strength_blocks",
                "block_bp",
                "mean_distinct_kmers_per_block",
            ]
        )
        for label in labels:
            blocks = blocks_by_cen.get(label, [])
            block_bp = sum(int(block["block_length_bp"]) for block in blocks)
            distinct_sum = sum(int(block["distinct_kmers"]) for block in blocks)
            mean_distinct = distinct_sum / len(blocks) if blocks else 0.0
            writer.writerow(
                [
                    label,
                    selected_by_cen.get(label, 0),
                    target_hits_by_cen.get(label, 0),
                    len(units_by_cen.get(label, [])),
                    len(blocks),
                    block_bp,
                    f"{mean_distinct:.3f}",
                ]
            )


def collect_strength_rows(
    labels: list[str],
    selected: dict[str, SelectedKmer],
    blocks_by_cen: dict[str, list[dict[str, object]]],
    units_by_cen: dict[str, list[dict[str, object]]],
) -> list[dict[str, object]]:
    selected_by_cen = Counter(call.assigned_cen for call in selected.values())
    target_hits_by_cen = Counter()
    for call in selected.values():
        target_hits_by_cen[call.assigned_cen] += call.target_hits

    rows = []
    for label in labels:
        blocks = blocks_by_cen.get(label, [])
        block_bp = sum(int(block["block_length_bp"]) for block in blocks)
        distinct_sum = sum(int(block["distinct_kmers"]) for block in blocks)
        mean_distinct = distinct_sum / len(blocks) if blocks else 0.0
        rows.append(
            {
                "assigned_cen": label,
                "selected_distinct_kmers": selected_by_cen.get(label, 0),
                "target_map_hits": target_hits_by_cen.get(label, 0),
                "cenhap_strength_units": len(units_by_cen.get(label, [])),
                "cenhap_strength_blocks": len(blocks),
                "block_bp": block_bp,
                "mean_distinct_kmers_per_block": mean_distinct,
            }
        )
    return rows


def build_region_strength_rows(
    strength_rows: list[dict[str, object]],
    interval_groups: dict[str, list[CoreInterval]],
) -> list[dict[str, object]]:
    rows = []
    for row in strength_rows:
        label = str(row["assigned_cen"])
        start, end = label_extent(label, interval_groups)
        length_bp = label_length(label, interval_groups)
        mb = length_bp / 1_000_000 if length_bp else 0.0
        selected_count = int(row["selected_distinct_kmers"])
        target_hits = int(row["target_map_hits"])
        units = int(row["cenhap_strength_units"])
        blocks = int(row["cenhap_strength_blocks"])
        block_bp = int(row["block_bp"])
        rows.append(
            {
                "assigned_cen": label,
                "region_class": label_region_class(label),
                "chr": label_chrom(label, interval_groups),
                "region_start": start,
                "region_end": end,
                "region_length_bp": length_bp,
                "selected_distinct_kmers": selected_count,
                "selected_distinct_kmers_per_mb": selected_count / mb if mb else 0.0,
                "target_map_hits": target_hits,
                "target_map_hits_per_mb": target_hits / mb if mb else 0.0,
                "cenhap_strength_units": units,
                "cenhap_strength_units_per_mb": units / mb if mb else 0.0,
                "cenhap_strength_blocks": blocks,
                "cenhap_strength_blocks_per_mb": blocks / mb if mb else 0.0,
                "block_bp": block_bp,
                "block_bp_fraction": block_bp / length_bp if length_bp else 0.0,
                "mean_distinct_kmers_per_block": float(row["mean_distinct_kmers_per_block"]),
            }
        )
    return rows


def write_region_strength(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "assigned_cen",
        "region_class",
        "chr",
        "region_start",
        "region_end",
        "region_length_bp",
        "selected_distinct_kmers",
        "selected_distinct_kmers_per_mb",
        "target_map_hits",
        "target_map_hits_per_mb",
        "cenhap_strength_units",
        "cenhap_strength_units_per_mb",
        "cenhap_strength_blocks",
        "cenhap_strength_blocks_per_mb",
        "block_bp",
        "block_bp_fraction",
        "mean_distinct_kmers_per_block",
    ]
    float_fields = {
        "selected_distinct_kmers_per_mb",
        "target_map_hits_per_mb",
        "cenhap_strength_units_per_mb",
        "cenhap_strength_blocks_per_mb",
        "block_bp_fraction",
        "mean_distinct_kmers_per_block",
    }
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            for key in float_fields:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def label_suffix(label: str) -> str:
    if label.startswith(PERICEN_PREFIX):
        return label[len(PERICEN_PREFIX) :]
    if label.startswith(CORE_PREFIX):
        return label[len(CORE_PREFIX) :]
    return label


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator:
        return numerator / denominator
    return math.inf if numerator else 0.0


def build_paired_cen_pericen_rows(
    region_strength_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    by_label = {str(row["assigned_cen"]): row for row in region_strength_rows}
    suffixes = sorted(
        {
            label_suffix(label)
            for label in by_label
            if label.startswith(CORE_PREFIX) or label.startswith(PERICEN_PREFIX)
        },
        key=lambda suffix: (0, int(suffix)) if suffix.isdigit() else (1, suffix),
    )
    rows = []
    for suffix in suffixes:
        cen_label = f"{CORE_PREFIX}{suffix}"
        pericen_label = f"{PERICEN_PREFIX}{suffix}"
        if cen_label not in by_label or pericen_label not in by_label:
            continue
        cen = by_label[cen_label]
        pericen = by_label[pericen_label]
        cen_units = float(cen["cenhap_strength_units_per_mb"])
        pericen_units = float(pericen["cenhap_strength_units_per_mb"])
        cen_selected = float(cen["selected_distinct_kmers_per_mb"])
        pericen_selected = float(pericen["selected_distinct_kmers_per_mb"])
        rows.append(
            {
                "chr": str(cen["chr"]),
                "cen_region": cen_label,
                "pericen_region": pericen_label,
                "cen_region_length_bp": int(cen["region_length_bp"]),
                "pericen_region_length_bp": int(pericen["region_length_bp"]),
                "cen_cenhap_strength_units_per_mb": cen_units,
                "pericen_cenhap_strength_units_per_mb": pericen_units,
                "pericen_to_cen_units_per_mb_ratio": safe_ratio(pericen_units, cen_units),
                "cen_selected_distinct_kmers_per_mb": cen_selected,
                "pericen_selected_distinct_kmers_per_mb": pericen_selected,
                "pericen_to_cen_selected_kmers_per_mb_ratio": safe_ratio(
                    pericen_selected, cen_selected
                ),
                "stronger_region_by_units_per_mb": (
                    pericen_label if pericen_units > cen_units else cen_label if cen_units > pericen_units else "tie"
                ),
            }
        )
    return rows


def write_paired_cen_pericen_strength(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "chr",
        "cen_region",
        "pericen_region",
        "cen_region_length_bp",
        "pericen_region_length_bp",
        "cen_cenhap_strength_units_per_mb",
        "pericen_cenhap_strength_units_per_mb",
        "pericen_to_cen_units_per_mb_ratio",
        "cen_selected_distinct_kmers_per_mb",
        "pericen_selected_distinct_kmers_per_mb",
        "pericen_to_cen_selected_kmers_per_mb_ratio",
        "stronger_region_by_units_per_mb",
    ]
    float_fields = {
        "cen_cenhap_strength_units_per_mb",
        "pericen_cenhap_strength_units_per_mb",
        "pericen_to_cen_units_per_mb_ratio",
        "cen_selected_distinct_kmers_per_mb",
        "pericen_selected_distinct_kmers_per_mb",
        "pericen_to_cen_selected_kmers_per_mb_ratio",
    }
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            for key in float_fields:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def relatedness_lookup(
    relatedness_pair_rows: list[dict[str, object]],
) -> dict[tuple[str, str], dict[str, object]]:
    lookup = {}
    for row in relatedness_pair_rows:
        left = str(row["cen_a"])
        right = str(row["cen_b"])
        lookup[tuple(sorted((left, right)))] = row
    return lookup


def relatedness_score(
    lookup: dict[tuple[str, str], dict[str, object]], left: str, right: str
) -> float:
    if left == right:
        return 1.0
    row = lookup.get(tuple(sorted((left, right))))
    return float(row["idf_weighted_jaccard"]) if row else 0.0


def mean_or_zero(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def build_cen_pericen_independence_rows(
    labels: list[str],
    relatedness_pair_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    lookup = relatedness_lookup(relatedness_pair_rows)
    cen_labels = [label for label in labels if label.startswith(CORE_PREFIX)]
    pericen_labels = [label for label in labels if label.startswith(PERICEN_PREFIX)]
    rows = []
    for cen_label in cen_labels:
        suffix = label_suffix(cen_label)
        pericen_label = f"{PERICEN_PREFIX}{suffix}"
        if pericen_label not in pericen_labels:
            continue
        same = relatedness_score(lookup, cen_label, pericen_label)
        cen_to_other_pericens = [
            relatedness_score(lookup, cen_label, label)
            for label in pericen_labels
            if label != pericen_label
        ]
        pericen_to_other_cens = [
            relatedness_score(lookup, pericen_label, label)
            for label in cen_labels
            if label != cen_label
        ]
        pericen_to_other_pericens = [
            relatedness_score(lookup, pericen_label, label)
            for label in pericen_labels
            if label != pericen_label
        ]
        cen_to_other_cens = [
            relatedness_score(lookup, cen_label, label) for label in cen_labels if label != cen_label
        ]
        mean_cen_other_pericen = mean_or_zero(cen_to_other_pericens)
        mean_pericen_other_cen = mean_or_zero(pericen_to_other_cens)
        rows.append(
            {
                "chr": f"Chr{suffix}",
                "cen_region": cen_label,
                "pericen_region": pericen_label,
                "same_chrom_cen_pericen_relatedness": same,
                "mean_cen_to_other_pericens_relatedness": mean_cen_other_pericen,
                "same_chrom_vs_cen_to_other_pericens_ratio": safe_ratio(
                    same, mean_cen_other_pericen
                ),
                "mean_pericen_to_other_cens_relatedness": mean_pericen_other_cen,
                "same_chrom_vs_pericen_to_other_cens_ratio": safe_ratio(
                    same, mean_pericen_other_cen
                ),
                "mean_pericen_to_other_pericens_relatedness": mean_or_zero(
                    pericen_to_other_pericens
                ),
                "mean_cen_to_other_cens_relatedness": mean_or_zero(cen_to_other_cens),
            }
        )
    return rows


def write_cen_pericen_independence(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "chr",
        "cen_region",
        "pericen_region",
        "same_chrom_cen_pericen_relatedness",
        "mean_cen_to_other_pericens_relatedness",
        "same_chrom_vs_cen_to_other_pericens_ratio",
        "mean_pericen_to_other_cens_relatedness",
        "same_chrom_vs_pericen_to_other_cens_ratio",
        "mean_pericen_to_other_pericens_relatedness",
        "mean_cen_to_other_cens_relatedness",
    ]
    float_fields = set(fieldnames[3:])
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = row.copy()
            for key in float_fields:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def write_strength_plot(path: Path, strength_rows: list[dict[str, object]]) -> str:
    labels = [str(row["assigned_cen"]) for row in strength_rows]
    includes_pericen = any(label.startswith(PERICEN_PREFIX) for label in labels)
    label_count = len(labels)
    rotate_x_labels = label_count >= 8
    label_angle = -45 if rotate_x_labels else 0
    label_anchor = "end" if rotate_x_labels else "middle"
    label_y_offset = 24 if rotate_x_labels else 19
    axis_title_y_offset = 82 if rotate_x_labels else 43
    metrics = [
        ("cenhap_strength_units", "Non-redundant cenhap units"),
        ("selected_distinct_kmers", "Selected distinct k-mers"),
        ("target_map_hits", "Assigned core-CEN map hits"),
    ]
    colors = ["#2f6f73", "#9b5d2e", "#5f6f95"]
    width = 1320
    height = 520 if rotate_x_labels else 470
    margin_left = 62
    margin_top = 74
    panel_gap = 46
    panel_width = (width - margin_left - 32 - panel_gap * (len(metrics) - 1)) / len(metrics)
    panel_height = 285
    baseline = margin_top + panel_height

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f1f1f}",
        ".title{font-size:22px;font-weight:700}",
        ".panel-title{font-size:15px;font-weight:700}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#d8d8d8;stroke-width:1}",
        ".tick{font-size:11px;fill:#555}",
        ".value{font-size:11px;font-weight:700}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="660" y="34" text-anchor="middle">'
        f'{"Cenhap Strength by Analysis Region" if includes_pericen else "Cenhap Strength by Core Centromere"}</text>',
    ]

    for metric_idx, ((metric, title), color) in enumerate(zip(metrics, colors)):
        panel_x = margin_left + metric_idx * (panel_width + panel_gap)
        values = [int(row[metric]) for row in strength_rows]
        top = max(values) if values else 1
        tick_step = max(1, round(top / 4))
        y_ticks = [0, tick_step, tick_step * 2, tick_step * 3, top]
        y_ticks = sorted(set(y_ticks))

        svg.append(
            f'<text class="panel-title" x="{panel_x + panel_width / 2:.1f}" y="61" '
            f'text-anchor="middle">{html.escape(title)}</text>'
        )
        svg.append(
            f'<line class="axis" x1="{panel_x:.1f}" y1="{baseline}" '
            f'x2="{panel_x + panel_width:.1f}" y2="{baseline}"/>'
        )
        svg.append(
            f'<line class="axis" x1="{panel_x:.1f}" y1="{margin_top}" '
            f'x2="{panel_x:.1f}" y2="{baseline}"/>'
        )
        for tick in y_ticks:
            y = baseline - (tick / top) * panel_height if top else baseline
            svg.append(
                f'<line class="grid" x1="{panel_x:.1f}" y1="{y:.1f}" '
                f'x2="{panel_x + panel_width:.1f}" y2="{y:.1f}"/>'
            )
            svg.append(
                f'<text class="tick" x="{panel_x - 8:.1f}" y="{y + 4:.1f}" '
                f'text-anchor="end">{tick:,}</text>'
            )

        bar_gap = 12
        bar_width = (panel_width - bar_gap * (len(labels) + 1)) / len(labels)
        for idx, (label, value) in enumerate(zip(labels, values)):
            bar_x = panel_x + bar_gap + idx * (bar_width + bar_gap)
            bar_height = (value / top) * panel_height if top else 0
            bar_y = baseline - bar_height
            label_x = bar_x + bar_width / 2
            svg.append(
                f'<rect x="{bar_x:.1f}" y="{bar_y:.1f}" width="{bar_width:.1f}" '
                f'height="{bar_height:.1f}" fill="{color}" stroke="#222" stroke-width="0.7"/>'
            )
            svg.append(
                f'<text class="value" x="{label_x:.1f}" y="{bar_y - 7:.1f}" '
                f'text-anchor="middle">{value:,}</text>'
            )
            transform = (
                f' transform="rotate({label_angle} {label_x:.1f} {baseline + label_y_offset:.1f})"'
                if rotate_x_labels
                else ""
            )
            svg.append(
                f'<text class="tick" x="{label_x:.1f}" y="{baseline + label_y_offset}" '
                f'text-anchor="{label_anchor}"{transform}>{html.escape(label)}</text>'
            )
        svg.append(
            f'<text class="tick" x="{panel_x + panel_width / 2:.1f}" y="{baseline + axis_title_y_offset}" '
            f'text-anchor="middle">{"Analysis region" if includes_pericen else "Core centromere"}</text>'
        )

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")
    return ""


def format_ratio(value: float) -> str:
    if math.isinf(value):
        return "inf"
    return f"{value:.2f}x"


def write_paired_strength_plot(path: Path, rows: list[dict[str, object]]) -> str:
    if not rows:
        return ""
    width = 1120
    row_height = 54
    margin_left = 126
    margin_right = 178
    margin_top = 78
    margin_bottom = 52
    plot_width = width - margin_left - margin_right
    height = margin_top + row_height * len(rows) + margin_bottom
    max_value = max(
        max(
            float(row["cen_cenhap_strength_units_per_mb"]),
            float(row["pericen_cenhap_strength_units_per_mb"]),
        )
        for row in rows
    )
    top = max_value if max_value else 1.0
    cen_color = "#2f6f73"
    pericen_color = "#b97833"

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f1f1f}",
        ".title{font-size:22px;font-weight:700}",
        ".tick{font-size:11px;fill:#555}",
        ".label{font-size:13px;font-weight:700}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#dddddd;stroke-width:1}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="560" y="34" text-anchor="middle">'
        "CEN vs periCEN Cenhap Strength</text>",
        '<text class="tick" x="560" y="56" text-anchor="middle">'
        "Non-redundant units per Mb; ratio is periCEN/CEN</text>",
    ]

    for tick_idx in range(5):
        value = top * tick_idx / 4
        x = margin_left + (value / top) * plot_width
        svg.append(
            f'<line class="grid" x1="{x:.1f}" y1="{margin_top - 16}" '
            f'x2="{x:.1f}" y2="{height - margin_bottom + 12}"/>'
        )
        svg.append(
            f'<text class="tick" x="{x:.1f}" y="{height - margin_bottom + 32}" '
            f'text-anchor="middle">{value:.1f}</text>'
        )
    svg.append(
        f'<line class="axis" x1="{margin_left}" y1="{height - margin_bottom + 12}" '
        f'x2="{margin_left + plot_width}" y2="{height - margin_bottom + 12}"/>'
    )

    for idx, row in enumerate(rows):
        y = margin_top + idx * row_height
        cen_value = float(row["cen_cenhap_strength_units_per_mb"])
        pericen_value = float(row["pericen_cenhap_strength_units_per_mb"])
        ratio = float(row["pericen_to_cen_units_per_mb_ratio"])
        chrom = str(row["chr"])
        svg.append(
            f'<text class="label" x="{margin_left - 16}" y="{y + 22}" '
            f'text-anchor="end">{html.escape(chrom)}</text>'
        )
        for sub_idx, (name, value, color) in enumerate(
            [("CEN", cen_value, cen_color), ("PERICEN", pericen_value, pericen_color)]
        ):
            bar_y = y + 2 + sub_idx * 21
            bar_w = (value / top) * plot_width if top else 0.0
            svg.append(
                f'<text class="tick" x="{margin_left - 10}" y="{bar_y + 12}" '
                f'text-anchor="end">{name}</text>'
            )
            svg.append(
                f'<rect x="{margin_left}" y="{bar_y}" width="{bar_w:.1f}" height="15" '
                f'fill="{color}" opacity="0.86">'
                f'<title>{html.escape(str(row["cen_region"] if name == "CEN" else row["pericen_region"]))}: '
                f'{value:.3f} units per Mb</title></rect>'
            )
            if value:
                svg.append(
                    f'<text class="tick" x="{margin_left + bar_w + 5:.1f}" y="{bar_y + 12}" '
                    f'text-anchor="start">{value:.1f}</text>'
                )
        svg.append(
            f'<text class="label" x="{margin_left + plot_width + 22}" y="{y + 29}" '
            f'text-anchor="start">{format_ratio(ratio)}</text>'
        )

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")
    return ""


def write_cen_pericen_scatter_plot(path: Path, rows: list[dict[str, object]]) -> str:
    if not rows:
        return ""
    width = 760
    height = 700
    margin_left = 86
    margin_right = 38
    margin_top = 74
    margin_bottom = 76
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    max_value = max(
        max(
            float(row["cen_cenhap_strength_units_per_mb"]),
            float(row["pericen_cenhap_strength_units_per_mb"]),
        )
        for row in rows
    )
    top = max_value * 1.08 if max_value else 1.0

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f1f1f}",
        ".title{font-size:22px;font-weight:700}",
        ".tick{font-size:11px;fill:#555}",
        ".label{font-size:12px;font-weight:700}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#dddddd;stroke-width:1}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="380" y="34" text-anchor="middle">'
        "CEN vs periCEN Strength Scatter</text>",
        '<text class="tick" x="380" y="56" text-anchor="middle">'
        "Points above the diagonal have stronger periCEN signal</text>",
    ]
    baseline = margin_top + plot_height
    for tick_idx in range(5):
        value = top * tick_idx / 4
        x = margin_left + (value / top) * plot_width
        y = baseline - (value / top) * plot_height
        svg.append(
            f'<line class="grid" x1="{x:.1f}" y1="{margin_top}" '
            f'x2="{x:.1f}" y2="{baseline}"/>'
        )
        svg.append(
            f'<line class="grid" x1="{margin_left}" y1="{y:.1f}" '
            f'x2="{margin_left + plot_width}" y2="{y:.1f}"/>'
        )
        svg.append(
            f'<text class="tick" x="{x:.1f}" y="{baseline + 20}" '
            f'text-anchor="middle">{value:.1f}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left - 8}" y="{y + 4:.1f}" '
            f'text-anchor="end">{value:.1f}</text>'
        )
    svg.append(
        f'<line class="axis" x1="{margin_left}" y1="{baseline}" '
        f'x2="{margin_left + plot_width}" y2="{baseline}"/>'
    )
    svg.append(
        f'<line class="axis" x1="{margin_left}" y1="{margin_top}" '
        f'x2="{margin_left}" y2="{baseline}"/>'
    )
    svg.append(
        f'<line x1="{margin_left}" y1="{baseline}" '
        f'x2="{margin_left + plot_width}" y2="{margin_top}" '
        'stroke="#777" stroke-width="1.5" stroke-dasharray="5 5"/>'
    )
    for row in rows:
        cen_value = float(row["cen_cenhap_strength_units_per_mb"])
        pericen_value = float(row["pericen_cenhap_strength_units_per_mb"])
        x = margin_left + (cen_value / top) * plot_width
        y = baseline - (pericen_value / top) * plot_height
        suffix = label_suffix(str(row["cen_region"]))
        svg.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="#5f6f95" opacity="0.88">'
            f'<title>{html.escape(str(row["chr"]))}: CEN {cen_value:.3f}, '
            f'PERICEN {pericen_value:.3f} units per Mb</title></circle>'
        )
        svg.append(
            f'<text class="label" x="{x + 8:.1f}" y="{y - 8:.1f}" '
            f'text-anchor="start">{html.escape(suffix)}</text>'
        )
    svg.append(
        f'<text class="tick" x="{margin_left + plot_width / 2:.1f}" y="{height - 24}" '
        'text-anchor="middle">CEN units per Mb</text>'
    )
    svg.append(
        f'<text class="tick" x="22" y="{margin_top + plot_height / 2:.1f}" '
        'text-anchor="middle" transform="rotate(-90 22 '
        f'{margin_top + plot_height / 2:.1f})">periCEN units per Mb</text>'
    )
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")
    return ""


def overlaps(left_start: int, left_end: int, right_start: int, right_end: int) -> bool:
    return left_start <= right_end and right_start <= left_end


def window_starts(start: int, end: int, window_size: int, window_step: int) -> list[int]:
    if window_size <= 0 or window_step <= 0:
        raise SystemExit("--window-size and --window-step must be positive integers.")
    if window_size >= end - start + 1:
        return [start]
    starts = list(range(start, end - window_size + 2, window_step))
    final_start = end - window_size + 1
    if starts[-1] != final_start:
        starts.append(final_start)
    return starts


def build_window_rows(
    intervals: list[CoreInterval],
    stats: dict[str, KmerStats],
    selected: dict[str, SelectedKmer],
    blocks_by_cen: dict[str, list[dict[str, object]]],
    kmer_size: int,
    window_size: int,
    window_step: int,
) -> list[dict[str, object]]:
    target_hits_by_cen: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for kmer, call in selected.items():
        for pos in stats[kmer].target_positions.get(call.assigned_cen, []):
            target_hits_by_cen[call.assigned_cen].append((pos, pos + kmer_size - 1, kmer))

    rows = []
    for interval in intervals:
        label = interval.label
        target_hits = target_hits_by_cen.get(label, [])
        blocks = blocks_by_cen.get(label, [])
        for start in window_starts(interval.start, interval.end, window_size, window_step):
            end = min(start + window_size - 1, interval.end)
            selected_kmers = set()
            target_map_hits = 0
            for hit_start, hit_end, kmer in target_hits:
                if overlaps(hit_start, hit_end, start, end):
                    selected_kmers.add(kmer)
                    target_map_hits += 1

            block_count = 0
            block_bp = 0
            unit_ids = set()
            for block in blocks:
                block_start = int(block["block_start"])
                block_end = int(block["block_end"])
                if overlaps(block_start, block_end, start, end):
                    block_count += 1
                    block_bp += min(block_end, end) - max(block_start, start) + 1
                    unit_ids.update(block.get("_unit_ids", set()))

            window_length = end - start + 1
            scale = 100_000 / window_length
            rows.append(
                {
                    "assigned_cen": label,
                    "region_class": label_region_class(label),
                    "region_part": interval.region_name,
                    "chr": interval.chrom,
                    "window_start": start,
                    "window_end": end,
                    "window_midpoint": (start + end) // 2,
                    "window_length_bp": window_length,
                    "selected_distinct_kmers": len(selected_kmers),
                    "target_map_hits": target_map_hits,
                    "cenhap_strength_units": len(unit_ids),
                    "cenhap_strength_blocks": block_count,
                    "block_bp": block_bp,
                    "selected_kmers_per_100kb": len(selected_kmers) * scale,
                    "target_map_hits_per_100kb": target_map_hits * scale,
                    "units_per_100kb": len(unit_ids) * scale,
                    "blocks_per_100kb": block_count * scale,
                }
            )
    return rows


def build_bin_rows(
    intervals: list[CoreInterval],
    stats: dict[str, KmerStats],
    selected: dict[str, SelectedKmer],
    kmer_size: int,
    bin_size: int,
) -> list[dict[str, object]]:
    if bin_size <= 0:
        raise SystemExit("--bin-size must be a positive integer.")

    interval_groups = intervals_by_label(intervals)
    bins_by_label: dict[str, list[tuple[int, CoreInterval, int, int]]] = defaultdict(list)
    for label in unique_labels(intervals):
        bin_index = 1
        for interval in interval_groups.get(label, []):
            for bin_start, bin_end in make_bins(interval, bin_size):
                bins_by_label[label].append((bin_index, interval, bin_start, bin_end))
                bin_index += 1
    kmers_by_bin: dict[tuple[str, int], set[str]] = defaultdict(set)
    hits_by_bin: Counter = Counter()

    for kmer, call in selected.items():
        bins = bins_by_label[call.assigned_cen]
        for pos in stats[kmer].target_positions.get(call.assigned_cen, []):
            for bin_index, interval, bin_start, bin_end in bins:
                hit_start = max(pos, interval.start)
                hit_end = min(pos + kmer_size - 1, interval.end)
                if hit_start > hit_end:
                    continue
                if overlaps(hit_start, hit_end, bin_start, bin_end):
                    kmers_by_bin[(call.assigned_cen, bin_index)].add(kmer)
                    hits_by_bin[(call.assigned_cen, bin_index)] += 1

    rows = []
    for label in unique_labels(intervals):
        for bin_index, interval, start, end in bins_by_label[label]:
            length = end - start + 1
            selected_count = len(kmers_by_bin.get((label, bin_index), set()))
            target_hits = hits_by_bin.get((label, bin_index), 0)
            rows.append(
                {
                    "assigned_cen": label,
                    "region_class": label_region_class(label),
                    "region_part": interval.region_name,
                    "chr": interval.chrom,
                    "bin_index": bin_index,
                    "bin_start": start,
                    "bin_end": end,
                    "bin_midpoint": (start + end) // 2,
                    "bin_length_bp": length,
                    "selected_distinct_kmers": selected_count,
                    "target_map_hits": target_hits,
                    "selected_kmers_per_100kb": selected_count * (100_000 / length),
                    "target_map_hits_per_100kb": target_hits * (100_000 / length),
                }
            )
    return rows


def make_bins(interval: CoreInterval, bin_size: int) -> list[tuple[int, int]]:
    bins = []
    start = interval.start
    while start <= interval.end:
        end = min(start + bin_size - 1, interval.end)
        bins.append((start, end))
        start = end + 1

    if len(bins) > 1:
        last_start, last_end = bins[-1]
        if last_end - last_start + 1 < bin_size * 0.10:
            prev_start, _ = bins[-2]
            bins[-2] = (prev_start, last_end)
            bins.pop()
    return bins


def gini(values: list[int]) -> float:
    if not values or sum(values) == 0:
        return 0.0
    ranked = sorted(values)
    n = len(ranked)
    weighted = sum((idx + 1) * value for idx, value in enumerate(ranked))
    return (2 * weighted) / (n * sum(ranked)) - (n + 1) / n


def top_fraction(values: list[int], fraction: float) -> float:
    total = sum(values)
    if not values or total == 0:
        return 0.0
    n_top = max(1, int(len(values) * fraction + 0.999999))
    return sum(sorted(values, reverse=True)[:n_top]) / total


def build_dispersion_rows(
    intervals: list[CoreInterval],
    window_rows: list[dict[str, object]],
    window_size: int,
    window_step: int,
) -> list[dict[str, object]]:
    rows_by_cen: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in window_rows:
            rows_by_cen[str(row["assigned_cen"])].append(row)

    dispersion = []
    interval_groups = intervals_by_label(intervals)
    for label in unique_labels(intervals):
        rows = rows_by_cen.get(label, [])
        unit_values = [int(row["cenhap_strength_units"]) for row in rows]
        block_values = [int(row["cenhap_strength_blocks"]) for row in rows]
        windows = len(rows)
        unit_total = sum(unit_values)
        block_total = sum(block_values)
        start, end = label_extent(label, interval_groups)
        dispersion.append(
            {
                "assigned_cen": label,
                "region_class": label_region_class(label),
                "chr": label_chrom(label, interval_groups),
                "core_start": start,
                "core_end": end,
                "core_length_bp": label_length(label, interval_groups),
                "window_size": window_size,
                "window_step": window_step,
                "windows": windows,
                "windows_with_units": sum(1 for value in unit_values if value > 0),
                "fraction_windows_with_units": (
                    sum(1 for value in unit_values if value > 0) / windows if windows else 0.0
                ),
                "total_window_units": unit_total,
                "max_window_units": max(unit_values) if unit_values else 0,
                "mean_window_units": unit_total / windows if windows else 0.0,
                "top_10pct_windows_fraction_of_units": top_fraction(unit_values, 0.10),
                "gini_window_units": gini(unit_values),
                "windows_with_blocks": sum(1 for value in block_values if value > 0),
                "fraction_windows_with_blocks": (
                    sum(1 for value in block_values if value > 0) / windows if windows else 0.0
                ),
                "total_window_blocks": block_total,
                "max_window_blocks": max(block_values) if block_values else 0,
                "mean_window_blocks": block_total / windows if windows else 0.0,
                "top_10pct_windows_fraction_of_blocks": top_fraction(block_values, 0.10),
                "gini_window_blocks": gini(block_values),
            }
        )
    return dispersion


def write_window_rows(path: Path, window_rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "assigned_cen",
        "region_class",
        "region_part",
        "chr",
        "window_start",
        "window_end",
        "window_midpoint",
        "window_length_bp",
        "selected_distinct_kmers",
        "target_map_hits",
        "cenhap_strength_units",
        "cenhap_strength_blocks",
        "block_bp",
        "selected_kmers_per_100kb",
        "target_map_hits_per_100kb",
        "units_per_100kb",
        "blocks_per_100kb",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in window_rows:
            formatted = row.copy()
            for key in [
                "selected_kmers_per_100kb",
                "target_map_hits_per_100kb",
                "units_per_100kb",
                "blocks_per_100kb",
            ]:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def write_bin_rows(path: Path, bin_rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "assigned_cen",
        "region_class",
        "region_part",
        "chr",
        "bin_index",
        "bin_start",
        "bin_end",
        "bin_midpoint",
        "bin_length_bp",
        "selected_distinct_kmers",
        "target_map_hits",
        "selected_kmers_per_100kb",
        "target_map_hits_per_100kb",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in bin_rows:
            formatted = row.copy()
            formatted["selected_kmers_per_100kb"] = (
                f"{float(formatted['selected_kmers_per_100kb']):.6f}"
            )
            formatted["target_map_hits_per_100kb"] = (
                f"{float(formatted['target_map_hits_per_100kb']):.6f}"
            )
            writer.writerow(formatted)


def write_dispersion_rows(path: Path, dispersion_rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "assigned_cen",
        "region_class",
        "chr",
        "core_start",
        "core_end",
        "core_length_bp",
        "window_size",
        "window_step",
        "windows",
        "windows_with_units",
        "fraction_windows_with_units",
        "total_window_units",
        "max_window_units",
        "mean_window_units",
        "top_10pct_windows_fraction_of_units",
        "gini_window_units",
        "windows_with_blocks",
        "fraction_windows_with_blocks",
        "total_window_blocks",
        "max_window_blocks",
        "mean_window_blocks",
        "top_10pct_windows_fraction_of_blocks",
        "gini_window_blocks",
    ]
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in dispersion_rows:
            formatted = row.copy()
            for key in [
                "fraction_windows_with_units",
                "mean_window_units",
                "top_10pct_windows_fraction_of_units",
                "gini_window_units",
                "fraction_windows_with_blocks",
                "mean_window_blocks",
                "top_10pct_windows_fraction_of_blocks",
                "gini_window_blocks",
            ]:
                formatted[key] = f"{float(formatted[key]):.6f}"
            writer.writerow(formatted)


def write_window_plot(
    path: Path,
    intervals: list[CoreInterval],
    window_rows: list[dict[str, object]],
) -> str:
    rows_by_cen: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in window_rows:
        rows_by_cen[str(row["assigned_cen"])].append(row)

    width = 1320
    row_height = 110
    labels = unique_labels(intervals)
    interval_groups = intervals_by_label(intervals)
    margin_left = 88
    margin_right = 32
    margin_top = 70
    plot_width = width - margin_left - margin_right
    height = margin_top + row_height * len(labels) + 52
    max_units = max((int(row["cenhap_strength_units"]) for row in window_rows), default=1)
    max_blocks = max((int(row["cenhap_strength_blocks"]) for row in window_rows), default=1)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f1f1f}",
        ".title{font-size:22px;font-weight:700}",
        ".label{font-size:14px;font-weight:700}",
        ".tick{font-size:11px;fill:#555}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#d8d8d8;stroke-width:1}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="660" y="34" text-anchor="middle">'
        "Local Cenhap Strength Across Analysis Regions</text>",
        '<text class="tick" x="660" y="56" text-anchor="middle">'
        "Blue bars: non-redundant units per window; amber line: physical blocks per window</text>",
    ]

    for idx, label in enumerate(labels):
        row_y = margin_top + idx * row_height
        baseline = row_y + 70
        top_y = row_y + 8
        rows = rows_by_cen.get(label, [])
        region_start, region_end = label_extent(label, interval_groups)
        region_span = int(region_end) - int(region_start) + 1
        region_len = label_length(label, interval_groups)
        svg.append(
            f'<text class="label" x="{margin_left - 16}" y="{baseline - 28}" '
            f'text-anchor="end">{html.escape(label)}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left - 16}" y="{baseline - 10}" '
            f'text-anchor="end">{region_len / 1_000_000:.2f} Mb</text>'
        )
        svg.append(
            f'<line class="axis" x1="{margin_left}" y1="{baseline}" '
            f'x2="{margin_left + plot_width}" y2="{baseline}"/>'
        )
        svg.append(
            f'<line class="grid" x1="{margin_left}" y1="{top_y}" '
            f'x2="{margin_left + plot_width}" y2="{top_y}"/>'
        )
        points = []
        for row in rows:
            start = int(row["window_start"])
            end = int(row["window_end"])
            midpoint = int(row["window_midpoint"])
            units = int(row["cenhap_strength_units"])
            blocks = int(row["cenhap_strength_blocks"])
            x = margin_left + ((start - int(region_start)) / region_span) * plot_width
            x_end = margin_left + ((end - int(region_start) + 1) / region_span) * plot_width
            bar_w = max(1.0, x_end - x)
            bar_h = (units / max_units) * 58 if max_units else 0
            svg.append(
                f'<rect x="{x:.1f}" y="{baseline - bar_h:.1f}" width="{bar_w:.1f}" '
                f'height="{bar_h:.1f}" fill="#2f6f73" opacity="0.72"/>'
            )
            point_x = margin_left + ((midpoint - int(region_start)) / region_span) * plot_width
            point_y = baseline - ((blocks / max_blocks) * 58 if max_blocks else 0)
            points.append(f"{point_x:.1f},{point_y:.1f}")
        if points:
            svg.append(
                f'<polyline points="{" ".join(points)}" fill="none" stroke="#b97833" '
                'stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/>'
            )
        svg.append(
            f'<text class="tick" x="{margin_left}" y="{baseline + 18}" text-anchor="start">'
            f'{int(region_start):,}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left + plot_width}" y="{baseline + 18}" '
            f'text-anchor="end">{int(region_end):,}</text>'
        )

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")
    return ""


def write_bin_plot(
    path: Path,
    intervals: list[CoreInterval],
    bin_rows: list[dict[str, object]],
    bin_size: int,
) -> str:
    rows_by_cen: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in bin_rows:
        rows_by_cen[str(row["assigned_cen"])].append(row)

    width = 1320
    row_height = 108
    labels = unique_labels(intervals)
    interval_groups = intervals_by_label(intervals)
    margin_left = 88
    margin_right = 36
    margin_top = 72
    plot_width = width - margin_left - margin_right
    height = margin_top + row_height * len(labels) + 54
    max_count = max((int(row["selected_distinct_kmers"]) for row in bin_rows), default=1)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#1f1f1f}",
        ".title{font-size:22px;font-weight:700}",
        ".label{font-size:14px;font-weight:700}",
        ".tick{font-size:11px;fill:#555}",
        ".axis{stroke:#333;stroke-width:1}",
        ".grid{stroke:#d8d8d8;stroke-width:1}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="660" y="34" text-anchor="middle">'
        "Cenhap-Defining K-mers by Fixed Region Bin</text>",
        f'<text class="tick" x="660" y="56" text-anchor="middle">'
        f"Bar height is distinct selected k-mers per {bin_size:,} bp bin</text>",
    ]

    for idx, label in enumerate(labels):
        row_y = margin_top + idx * row_height
        baseline = row_y + 68
        top_y = row_y + 8
        rows = rows_by_cen.get(label, [])
        region_start, region_end = label_extent(label, interval_groups)
        region_span = int(region_end) - int(region_start) + 1
        region_len = label_length(label, interval_groups)
        svg.append(
            f'<text class="label" x="{margin_left - 16}" y="{baseline - 30}" '
            f'text-anchor="end">{html.escape(label)}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left - 16}" y="{baseline - 12}" '
            f'text-anchor="end">{region_len / 1_000_000:.2f} Mb</text>'
        )
        svg.append(
            f'<line class="axis" x1="{margin_left}" y1="{baseline}" '
            f'x2="{margin_left + plot_width}" y2="{baseline}"/>'
        )
        svg.append(
            f'<line class="grid" x1="{margin_left}" y1="{top_y}" '
            f'x2="{margin_left + plot_width}" y2="{top_y}"/>'
        )
        for row in rows:
            start = int(row["bin_start"])
            end = int(row["bin_end"])
            count = int(row["selected_distinct_kmers"])
            x = margin_left + ((start - int(region_start)) / region_span) * plot_width
            x_end = margin_left + ((end - int(region_start) + 1) / region_span) * plot_width
            bar_w = max(1.0, x_end - x - 1.0)
            bar_h = (count / max_count) * 58 if max_count else 0
            svg.append(
                f'<rect x="{x:.1f}" y="{baseline - bar_h:.1f}" width="{bar_w:.1f}" '
                f'height="{bar_h:.1f}" fill="#2f6f73" opacity="0.82">'
                f'<title>{html.escape(str(row["assigned_cen"]))} bin {row["bin_index"]}: '
                f'{count:,} distinct selected k-mers, {int(row["target_map_hits"]):,} map hits</title>'
                "</rect>"
            )
            if count:
                svg.append(
                    f'<text class="tick" x="{x + bar_w / 2:.1f}" y="{baseline - bar_h - 4:.1f}" '
                    f'text-anchor="middle">{count:,}</text>'
                )
        svg.append(
            f'<text class="tick" x="{margin_left}" y="{baseline + 18}" text-anchor="start">'
            f'{int(region_start):,}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left + plot_width}" y="{baseline + 18}" '
            f'text-anchor="end">{int(region_end):,}</text>'
        )

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")
    return ""


def write_selected_kmers(path: Path, selected: dict[str, SelectedKmer]) -> None:
    with open(path, "w") as out:
        for kmer in sorted(selected):
            out.write(kmer + "\n")


def write_assigned_core_map(
    map_tsv: str,
    path: Path,
    selected: dict[str, SelectedKmer],
    interval_groups: dict[str, list[CoreInterval]],
) -> int:
    rows_written = 0
    with open(map_tsv, newline="") as f, open(path, "w", newline="") as out:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            kmer = row["k-mer"]
            call = selected.get(kmer)
            if call is None:
                continue
            pos = clean_int(row["pos"])
            chrom = normalize_chrom(row["chr"])
            if any(
                interval.chrom == chrom and interval.start <= pos <= interval.end
                for interval in interval_groups.get(call.assigned_cen, [])
            ):
                writer.writerow(row)
                rows_written += 1
    return rows_written


def write_stats(
    path: Path,
    args: argparse.Namespace,
    intervals: list[CoreInterval],
    fasta_lengths: dict[str, int],
    validation_warnings: list[str],
    map_rows: int,
    core_rows: int,
    kmer_count: int,
    selected: dict[str, SelectedKmer],
    blocks_by_cen: dict[str, list[dict[str, object]]],
    units_by_cen: dict[str, list[dict[str, object]]],
    fail_counts: Counter,
    assigned_map_rows: int,
    plot_status: str,
) -> None:
    with open(path, "w") as out:
        out.write(f"input_mode\t{getattr(args, 'input_mode', 'map_tsv')}\n")
        out.write(f"map_tsv\t{args.map_tsv}\n")
        if getattr(args, "generated_map_stats_tsv", ""):
            out.write(f"generated_map_stats_tsv\t{args.generated_map_stats_tsv}\n")
        out.write(f"coords\t{args.coords}\n")
        out.write(f"include_pericen\t{int(args.include_pericen)}\n")
        out.write(f"fasta\t{args.fasta or ''}\n")
        out.write(f"kmer_size\t{args.kmer_size}\n")
        out.write(f"kmer_step\t{args.kmer_step}\n")
        out.write(f"min_cen_count\t{args.min_cen_count}\n")
        out.write(f"max_outside_ratio\t{args.max_outside_ratio}\n")
        out.write(f"canonical\t{int(not args.no_canonical)}\n")
        out.write(f"min_target_hits\t{args.min_target_hits}\n")
        out.write(f"max_other_core_hits\t{args.max_other_core_hits}\n")
        out.write(f"min_target_core_fraction\t{args.min_target_core_fraction}\n")
        out.write(f"min_target_enrichment\t{args.min_target_enrichment}\n")
        out.write(f"max_map_hits\t{args.max_map_hits}\n")
        out.write(f"merge_gap\t{args.merge_gap}\n")
        out.write(f"window_size\t{args.window_size}\n")
        out.write(f"window_step\t{args.window_step}\n")
        out.write(f"bin_size\t{args.bin_size}\n")
        out.write(f"min_shared_cens\t{args.min_shared_cens}\n")
        out.write(f"min_shared_hits_per_cen\t{args.min_shared_hits_per_cen}\n")
        out.write(f"map_rows\t{map_rows}\n")
        out.write(f"core_centromere_map_rows\t{core_rows}\n")
        out.write(f"analysis_regions\t{len(unique_labels(intervals))}\n")
        out.write(f"analysis_intervals\t{len(intervals)}\n")
        out.write(f"distinct_kmers\t{kmer_count}\n")
        out.write(f"selected_kmers\t{len(selected)}\n")
        out.write(
            f"cenhap_strength_units_total\t{sum(len(units) for units in units_by_cen.values())}\n"
        )
        out.write(
            f"cenhap_strength_blocks_total\t{sum(len(blocks) for blocks in blocks_by_cen.values())}\n"
        )
        out.write(f"assigned_core_map_rows\t{assigned_map_rows}\n")
        if plot_status:
            out.write(f"{plot_status}\n")
        for interval in intervals:
            out.write(
                "core_interval\t"
                f"{interval.label}\t{interval.region_class}\t{interval.region_name}\t"
                f"{interval.chrom}\t{interval.start}\t{interval.end}\n"
            )
            if fasta_lengths.get(interval.chrom):
                out.write(f"fasta_length\t{interval.chrom}\t{fasta_lengths[interval.chrom]}\n")
        for warning in validation_warnings:
            out.write(f"warning\t{warning}\n")
        for reason, count in sorted(fail_counts.items()):
            out.write(f"failed_{reason}\t{count}\n")


def main() -> None:
    args = parse_args()
    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    intervals = read_core_intervals(args.coords, args.include_pericen)
    labels = unique_labels(intervals)
    interval_groups = intervals_by_label(intervals)
    interval_index = build_interval_index(intervals)

    if args.map_tsv:
        args.input_mode = "map_tsv"
        args.generated_map_stats_tsv = ""
        if args.kmer_size <= 0:
            args.kmer_size = infer_kmer_size(args.map_tsv)
    else:
        if not args.fasta:
            raise SystemExit("Provide either --map-tsv or --fasta with --coords.")
        args.input_mode = "generated_from_fasta"
        generated_map, generated_stats, _generated_summary = generate_core_kmer_map(
            args.fasta,
            intervals,
            prefix,
            args,
        )
        args.map_tsv = str(generated_map)
        args.generated_map_stats_tsv = str(generated_stats)

    fasta_lengths = read_fasta_lengths(args.fasta)
    validation_warnings = validate_intervals(intervals, fasta_lengths)

    stats, map_rows, core_rows = read_map_counts(args.map_tsv, interval_index)
    selected, fail_counts = select_kmers(stats, labels, args)
    shared_kmer_rows = build_shared_kmer_rows(
        stats,
        labels,
        selected,
        args.min_shared_cens,
        args.min_shared_hits_per_cen,
    )
    shared_kmer_set_rows = build_shared_kmer_set_rows(shared_kmer_rows)
    relatedness_pair_rows, relatedness_matrix = build_idf_relatedness_rows(stats, labels)
    blocks_by_cen = build_blocks(
        stats,
        selected,
        interval_groups,
        args.kmer_size,
        args.merge_gap,
    )
    units_by_cen = build_units(labels, selected, blocks_by_cen)
    unit_size_summary_rows = build_unit_size_summary(labels, units_by_cen)
    unit_size_distribution_rows = build_unit_size_distribution(labels, units_by_cen)
    bin_rows = build_bin_rows(
        intervals,
        stats,
        selected,
        args.kmer_size,
        args.bin_size,
    )
    window_rows = build_window_rows(
        intervals,
        stats,
        selected,
        blocks_by_cen,
        args.kmer_size,
        args.window_size,
        args.window_step,
    )
    dispersion_rows = build_dispersion_rows(
        intervals,
        window_rows,
        args.window_size,
        args.window_step,
    )

    selected_path = Path(str(prefix) + ".selected_kmers.txt")
    summary_path = Path(str(prefix) + ".kmer_summary.tsv")
    blocks_path = Path(str(prefix) + ".cenhap_blocks.tsv")
    units_path = Path(str(prefix) + ".cenhap_units.tsv")
    units_bed_path = Path(str(prefix) + ".cenhap_units.bed")
    unit_size_summary_path = Path(str(prefix) + ".cenhap_unit_size_summary.tsv")
    unit_size_distribution_path = Path(str(prefix) + ".cenhap_unit_size_distribution.tsv")
    unit_size_plot_path = Path(str(prefix) + ".cenhap_unit_size_distribution.svg")
    shared_kmers_path = Path(str(prefix) + ".cen_shared_kmers.tsv")
    shared_sets_path = Path(str(prefix) + ".cen_shared_kmer_sets.tsv")
    relatedness_pairs_path = Path(str(prefix) + ".cen_relatedness.idf_weighted_pairs.tsv")
    relatedness_matrix_path = Path(str(prefix) + ".cen_relatedness.idf_weighted_matrix.tsv")
    strength_path = Path(str(prefix) + ".cenhap_strength.tsv")
    region_strength_path = Path(str(prefix) + ".region_strength.tsv")
    paired_strength_path = Path(str(prefix) + ".paired_cen_pericen_strength.tsv")
    independence_path = Path(str(prefix) + ".cen_pericen_independence.tsv")
    plot_path = Path(str(prefix) + ".cenhap_strength_histogram.svg")
    paired_plot_path = Path(str(prefix) + ".paired_cen_pericen_strength.svg")
    scatter_plot_path = Path(str(prefix) + ".cen_pericen_strength_scatter.svg")
    bin_path = Path(str(prefix) + ".cenhap_bins.tsv")
    bin_plot_path = Path(str(prefix) + ".cenhap_bin_counts.svg")
    window_path = Path(str(prefix) + ".cenhap_windows.tsv")
    dispersion_path = Path(str(prefix) + ".cenhap_window_dispersion.tsv")
    window_plot_path = Path(str(prefix) + ".cenhap_local_strength.svg")
    assigned_map_path = Path(str(prefix) + ".assigned_core_map.tsv")
    stats_path = Path(str(prefix) + ".stats.txt")
    strength_rows = collect_strength_rows(labels, selected, blocks_by_cen, units_by_cen)
    region_strength_rows = build_region_strength_rows(strength_rows, interval_groups)
    paired_strength_rows = build_paired_cen_pericen_rows(region_strength_rows)
    independence_rows = build_cen_pericen_independence_rows(labels, relatedness_pair_rows)

    write_selected_kmers(selected_path, selected)
    write_summary(summary_path, stats, selected, labels, args)
    write_blocks(blocks_path, blocks_by_cen)
    write_units(units_path, units_by_cen)
    write_units_bed(units_bed_path, units_by_cen)
    write_unit_size_summary(unit_size_summary_path, unit_size_summary_rows)
    write_unit_size_distribution(unit_size_distribution_path, unit_size_distribution_rows)
    write_shared_kmer_rows(shared_kmers_path, shared_kmer_rows, labels)
    write_shared_kmer_set_rows(shared_sets_path, shared_kmer_set_rows)
    write_relatedness_pairs(relatedness_pairs_path, relatedness_pair_rows)
    write_relatedness_matrix(relatedness_matrix_path, relatedness_matrix)
    write_cen_strength(strength_path, labels, selected, blocks_by_cen, units_by_cen)
    write_region_strength(region_strength_path, region_strength_rows)
    write_paired_cen_pericen_strength(paired_strength_path, paired_strength_rows)
    write_cen_pericen_independence(independence_path, independence_rows)
    plot_status = "plot_skipped\t--skip-plot" if args.skip_plot else write_strength_plot(
        plot_path, strength_rows
    )
    unit_size_plot_status = "" if args.skip_plot else write_unit_size_plot(
        unit_size_plot_path,
        unit_size_distribution_rows,
    )
    plot_status = plot_status or unit_size_plot_status
    write_bin_rows(bin_path, bin_rows)
    bin_plot_status = "" if args.skip_plot else write_bin_plot(
        bin_plot_path,
        intervals,
        bin_rows,
        args.bin_size,
    )
    plot_status = plot_status or bin_plot_status
    paired_plot_status = "" if args.skip_plot or not paired_strength_rows else write_paired_strength_plot(
        paired_plot_path, paired_strength_rows
    )
    plot_status = plot_status or paired_plot_status
    scatter_plot_status = "" if args.skip_plot or not paired_strength_rows else write_cen_pericen_scatter_plot(
        scatter_plot_path, paired_strength_rows
    )
    plot_status = plot_status or scatter_plot_status
    write_window_rows(window_path, window_rows)
    write_dispersion_rows(dispersion_path, dispersion_rows)
    window_plot_status = "" if args.skip_plot or not args.write_window_plot else write_window_plot(
        window_plot_path, intervals, window_rows
    )
    plot_status = plot_status or window_plot_status
    assigned_map_rows = write_assigned_core_map(
        args.map_tsv,
        assigned_map_path,
        selected,
        interval_groups,
    )
    write_stats(
        stats_path,
        args,
        intervals,
        fasta_lengths,
        validation_warnings,
        map_rows,
        core_rows,
        len(stats),
        selected,
        blocks_by_cen,
        units_by_cen,
        fail_counts,
        assigned_map_rows,
        plot_status,
    )

    print(f"Wrote {selected_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {blocks_path}")
    print(f"Wrote {units_path}")
    print(f"Wrote {units_bed_path}")
    print(f"Wrote {unit_size_summary_path}")
    print(f"Wrote {unit_size_distribution_path}")
    print(f"Wrote {shared_kmers_path}")
    print(f"Wrote {shared_sets_path}")
    print(f"Wrote {relatedness_pairs_path}")
    print(f"Wrote {relatedness_matrix_path}")
    print(f"Wrote {strength_path}")
    print(f"Wrote {region_strength_path}")
    print(f"Wrote {paired_strength_path}")
    print(f"Wrote {independence_path}")
    if not plot_status:
        print(f"Wrote {plot_path}")
        print(f"Wrote {unit_size_plot_path}")
        print(f"Wrote {bin_plot_path}")
        if paired_strength_rows:
            print(f"Wrote {paired_plot_path}")
            print(f"Wrote {scatter_plot_path}")
        if args.write_window_plot:
            print(f"Wrote {window_plot_path}")
    print(f"Wrote {bin_path}")
    print(f"Wrote {window_path}")
    print(f"Wrote {dispersion_path}")
    print(f"Wrote {assigned_map_path}")
    print(f"Wrote {stats_path}")


if __name__ == "__main__":
    main()
