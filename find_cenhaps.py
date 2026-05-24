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
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


CORE_PREFIX = "CEN"
REVCOMP_TABLE = str.maketrans("ACGT", "TGCA")
VALID_DNA_RE = re.compile("^[ACGT]+$")


@dataclass
class CoreInterval:
    label: str
    chrom: str
    start: int
    end: int


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


def read_core_intervals(path: str) -> list[CoreInterval]:
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
                    )
                )
        elif {"chr", "region", "start", "end"}.issubset(fields):
            for row in reader:
                region = row[fields["region"]].strip().lower()
                if region not in {"cen", "centromere", "core", "core_centromere"}:
                    continue
                chrom = normalize_chrom(row[fields["chr"]])
                suffix = chrom[3:]
                intervals.append(
                    CoreInterval(
                        label=f"{CORE_PREFIX}{suffix}",
                        chrom=chrom,
                        start=clean_int(row[fields["start"]]),
                        end=clean_int(row[fields["end"]]),
                    )
                )
        else:
            raise SystemExit(
                "Coordinate file must have either Chr/Left/Right or chr/region/start/end columns."
            )

    if not intervals:
        raise SystemExit("No core centromere intervals were found.")
    intervals.sort(key=lambda item: (item.chrom, item.start, item.end))
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


def build_blocks(
    stats: dict[str, KmerStats],
    selected: dict[str, SelectedKmer],
    interval_by_label: dict[str, CoreInterval],
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
                    "chr": interval_by_label[label].chrom,
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
                writer.writerow({key: block[key] for key in fieldnames})


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
            physical_blocks = sum(
                1
                for block in blocks_by_cen.get(label, [])
                if unit_set.intersection(block.get("_kmers", set()))
            )
            target_hits = sum(selected[kmer].target_hits for kmer in unit_kmers)
            unit_rows.append(
                {
                    "assigned_cen": label,
                    "unit_id": "",
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


def write_strength_plot(path: Path, strength_rows: list[dict[str, object]]) -> str:
    labels = [str(row["assigned_cen"]) for row in strength_rows]
    metrics = [
        ("cenhap_strength_units", "Non-redundant cenhap units"),
        ("selected_distinct_kmers", "Selected distinct k-mers"),
        ("target_map_hits", "Assigned core-CEN map hits"),
    ]
    colors = ["#2f6f73", "#9b5d2e", "#5f6f95"]
    width = 1320
    height = 470
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
        "Cenhap Strength by Core Centromere</text>",
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
            svg.append(
                f'<rect x="{bar_x:.1f}" y="{bar_y:.1f}" width="{bar_width:.1f}" '
                f'height="{bar_height:.1f}" fill="{color}" stroke="#222" stroke-width="0.7"/>'
            )
            svg.append(
                f'<text class="value" x="{bar_x + bar_width / 2:.1f}" y="{bar_y - 7:.1f}" '
                f'text-anchor="middle">{value:,}</text>'
            )
            svg.append(
                f'<text class="tick" x="{bar_x + bar_width / 2:.1f}" y="{baseline + 19}" '
                f'text-anchor="middle">{html.escape(label)}</text>'
            )
        svg.append(
            f'<text class="tick" x="{panel_x + panel_width / 2:.1f}" y="{baseline + 43}" '
            'text-anchor="middle">Core centromere</text>'
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

    interval_by_label = {interval.label: interval for interval in intervals}
    bins_by_label = {interval.label: make_bins(interval, bin_size) for interval in intervals}
    kmers_by_bin: dict[tuple[str, int], set[str]] = defaultdict(set)
    hits_by_bin: Counter = Counter()

    for kmer, call in selected.items():
        interval = interval_by_label[call.assigned_cen]
        bins = bins_by_label[call.assigned_cen]
        for pos in stats[kmer].target_positions.get(call.assigned_cen, []):
            hit_start = max(pos, interval.start)
            hit_end = min(pos + kmer_size - 1, interval.end)
            for bin_index, (bin_start, bin_end) in enumerate(bins):
                if overlaps(hit_start, hit_end, bin_start, bin_end):
                    kmers_by_bin[(call.assigned_cen, bin_index)].add(kmer)
                    hits_by_bin[(call.assigned_cen, bin_index)] += 1

    rows = []
    for interval in intervals:
        label = interval.label
        for bin_index, (start, end) in enumerate(bins_by_label[label]):
            length = end - start + 1
            selected_count = len(kmers_by_bin.get((label, bin_index), set()))
            target_hits = hits_by_bin.get((label, bin_index), 0)
            rows.append(
                {
                    "assigned_cen": label,
                    "chr": interval.chrom,
                    "bin_index": bin_index + 1,
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
    for interval in intervals:
        rows = rows_by_cen.get(interval.label, [])
        unit_values = [int(row["cenhap_strength_units"]) for row in rows]
        block_values = [int(row["cenhap_strength_blocks"]) for row in rows]
        windows = len(rows)
        unit_total = sum(unit_values)
        block_total = sum(block_values)
        dispersion.append(
            {
                "assigned_cen": interval.label,
                "chr": interval.chrom,
                "core_start": interval.start,
                "core_end": interval.end,
                "core_length_bp": interval.end - interval.start + 1,
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
    margin_left = 88
    margin_right = 32
    margin_top = 70
    plot_width = width - margin_left - margin_right
    height = margin_top + row_height * len(intervals) + 52
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
        "Local Cenhap Strength Across Core Centromeres</text>",
        '<text class="tick" x="660" y="56" text-anchor="middle">'
        "Blue bars: non-redundant units per window; amber line: physical blocks per window</text>",
    ]

    for idx, interval in enumerate(intervals):
        row_y = margin_top + idx * row_height
        baseline = row_y + 70
        top_y = row_y + 8
        rows = rows_by_cen.get(interval.label, [])
        core_len = interval.end - interval.start + 1
        svg.append(
            f'<text class="label" x="{margin_left - 16}" y="{baseline - 28}" '
            f'text-anchor="end">{html.escape(interval.label)}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left - 16}" y="{baseline - 10}" '
            f'text-anchor="end">{core_len / 1_000_000:.2f} Mb</text>'
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
            x = margin_left + ((start - interval.start) / core_len) * plot_width
            x_end = margin_left + ((end - interval.start + 1) / core_len) * plot_width
            bar_w = max(1.0, x_end - x)
            bar_h = (units / max_units) * 58 if max_units else 0
            svg.append(
                f'<rect x="{x:.1f}" y="{baseline - bar_h:.1f}" width="{bar_w:.1f}" '
                f'height="{bar_h:.1f}" fill="#2f6f73" opacity="0.72"/>'
            )
            point_x = margin_left + ((midpoint - interval.start) / core_len) * plot_width
            point_y = baseline - ((blocks / max_blocks) * 58 if max_blocks else 0)
            points.append(f"{point_x:.1f},{point_y:.1f}")
        if points:
            svg.append(
                f'<polyline points="{" ".join(points)}" fill="none" stroke="#b97833" '
                'stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/>'
            )
        svg.append(
            f'<text class="tick" x="{margin_left}" y="{baseline + 18}" text-anchor="start">'
            f'{interval.start:,}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left + plot_width}" y="{baseline + 18}" '
            f'text-anchor="end">{interval.end:,}</text>'
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
    margin_left = 88
    margin_right = 36
    margin_top = 72
    plot_width = width - margin_left - margin_right
    height = margin_top + row_height * len(intervals) + 54
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
        "Cenhap-Defining K-mers by Fixed Core-CEN Bin</text>",
        f'<text class="tick" x="660" y="56" text-anchor="middle">'
        f"Bar height is distinct selected k-mers per {bin_size:,} bp bin</text>",
    ]

    for idx, interval in enumerate(intervals):
        row_y = margin_top + idx * row_height
        baseline = row_y + 68
        top_y = row_y + 8
        rows = rows_by_cen.get(interval.label, [])
        core_len = interval.end - interval.start + 1
        svg.append(
            f'<text class="label" x="{margin_left - 16}" y="{baseline - 30}" '
            f'text-anchor="end">{html.escape(interval.label)}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left - 16}" y="{baseline - 12}" '
            f'text-anchor="end">{core_len / 1_000_000:.2f} Mb</text>'
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
            x = margin_left + ((start - interval.start) / core_len) * plot_width
            x_end = margin_left + ((end - interval.start + 1) / core_len) * plot_width
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
            f'{interval.start:,}</text>'
        )
        svg.append(
            f'<text class="tick" x="{margin_left + plot_width}" y="{baseline + 18}" '
            f'text-anchor="end">{interval.end:,}</text>'
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
    interval_by_label: dict[str, CoreInterval],
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
            interval = interval_by_label[call.assigned_cen]
            if normalize_chrom(row["chr"]) != interval.chrom:
                continue
            pos = clean_int(row["pos"])
            if interval.start <= pos <= interval.end:
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
        out.write(f"map_rows\t{map_rows}\n")
        out.write(f"core_centromere_map_rows\t{core_rows}\n")
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
                f"core_interval\t{interval.label}\t{interval.chrom}\t{interval.start}\t{interval.end}\n"
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

    intervals = read_core_intervals(args.coords)
    labels = [interval.label for interval in intervals]
    interval_by_label = {interval.label: interval for interval in intervals}
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
    blocks_by_cen = build_blocks(
        stats,
        selected,
        interval_by_label,
        args.kmer_size,
        args.merge_gap,
    )
    units_by_cen = build_units(labels, selected, blocks_by_cen)
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
    strength_path = Path(str(prefix) + ".cenhap_strength.tsv")
    plot_path = Path(str(prefix) + ".cenhap_strength_histogram.svg")
    bin_path = Path(str(prefix) + ".cenhap_bins.tsv")
    bin_plot_path = Path(str(prefix) + ".cenhap_bin_counts.svg")
    window_path = Path(str(prefix) + ".cenhap_windows.tsv")
    dispersion_path = Path(str(prefix) + ".cenhap_window_dispersion.tsv")
    window_plot_path = Path(str(prefix) + ".cenhap_local_strength.svg")
    assigned_map_path = Path(str(prefix) + ".assigned_core_map.tsv")
    stats_path = Path(str(prefix) + ".stats.txt")
    strength_rows = collect_strength_rows(labels, selected, blocks_by_cen, units_by_cen)

    write_selected_kmers(selected_path, selected)
    write_summary(summary_path, stats, selected, labels, args)
    write_blocks(blocks_path, blocks_by_cen)
    write_units(units_path, units_by_cen)
    write_cen_strength(strength_path, labels, selected, blocks_by_cen, units_by_cen)
    plot_status = "plot_skipped\t--skip-plot" if args.skip_plot else write_strength_plot(
        plot_path, strength_rows
    )
    write_bin_rows(bin_path, bin_rows)
    bin_plot_status = "" if args.skip_plot else write_bin_plot(
        bin_plot_path,
        intervals,
        bin_rows,
        args.bin_size,
    )
    plot_status = plot_status or bin_plot_status
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
        interval_by_label,
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
    print(f"Wrote {strength_path}")
    if not plot_status:
        print(f"Wrote {plot_path}")
        print(f"Wrote {bin_plot_path}")
        if args.write_window_plot:
            print(f"Wrote {window_plot_path}")
    print(f"Wrote {bin_path}")
    print(f"Wrote {window_path}")
    print(f"Wrote {dispersion_path}")
    print(f"Wrote {assigned_map_path}")
    print(f"Wrote {stats_path}")


if __name__ == "__main__":
    main()
