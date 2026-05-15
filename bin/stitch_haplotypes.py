#!/usr/bin/env python3
"""
stitch_haplotypes.py — Post-hoc stitching of DEVIDER per-window haplotypes.

Why window-stitching is needed
------------------------------
DEVIDER v0.0.1 (the only tagged release) reconstructs haplotypes inside
overlapping genomic windows.  Each window emits its own set of candidate
haplotype sequences with locally consistent SNV calls.  v0.0.1 has **no
``--merge-windows`` flag**: it does not perform cross-window stitching.  For a
9.6 kb HCV genome that is window-wise reconstructed in chunks of (typically)
~2-3 kb, this means the user is left with a collection of per-window haplotypes
that need to be joined into longer-range haplotypes before downstream analysis
(phylogenetics, drug-resistance prediction, etc.).

This script implements that stitching as a post-hoc step on top of DEVIDER's
own outputs.  See ``docs/data_flow.md`` Step 5.20 and ``CLAUDE.md`` decision
D11 for the rationale; see ``docs/architecture_reasoning.md`` §9 for why we do
not try to force DEVIDER into one giant window.

What "spanning read" means here
-------------------------------
A spanning read is a single long ONT read whose primary alignment touches
**both** of two adjacent DEVIDER windows.  DEVIDER's ``--output-reads`` flag
emits a haplotype-tagged BAM in which each read carries the haplotype ID it
was assigned **within each window it spans**.  We look up that tag for both
windows of every adjacent pair and tally the ``(hap_in_W_i, hap_in_W_{i+1})``
pairings.  When ``count >= --min-reads`` we declare that ``hap_X → hap_Y`` is
a supported link between the two windows.

Greedy assembly + limitations
-----------------------------
We build merged haplotypes by greedily extending chains left-to-right through
the windows.  At each junction we follow the strongest supported link.  When
multiple links from the same upstream haplotype meet the threshold (e.g.
``hap A → hap X`` with 28 reads AND ``hap A → hap Y`` with 14 reads, both
above ``--min-reads``), we enumerate **all** of them as separate chains
rather than picking one arbitrarily.  This is conservative: it can over-emit
where the true biology has branching haplotype evolution, but it avoids
silently dropping minority-supported paths that may correspond to real
low-frequency variants.

The algorithm is NOT optimal — for true branching evolution with shared
upstream segments, downstream phylogenetic analysis is the right tool to
collapse near-identical merged haplotypes.  The abundance reported per
chain is a **lower bound**: it is the minimum of the per-junction support
fractions (``link_reads / total_spanning_reads``) across all junctions in
the chain.

What ``--min-reads`` controls
-----------------------------
``--min-reads`` (= ``params.stitch_min_reads``, default 5) is the minimum
number of spanning reads that must support a single ``hap_X → hap_Y`` link
for it to be accepted.  Tuning:

  * Low coverage / very short reads → lower (e.g. 3): you get more partial
    chains but at the cost of weakly-supported links that may be artefacts.
  * Deep coverage / very long reads → higher (e.g. 10-15): you get fewer
    but more confident chains.

DEVIDER failure handling
------------------------
If the DEVIDER process emitted ``devider.failed`` (or no haplotype FASTAs
exist) we emit an empty merged_haplotypes.fasta and a JSON report with
``fallback_used = true``.  The script always exits 0 in these cases — it is
not the stitcher's job to fail the sample.

DEVIDER output discovery
------------------------
DEVIDER's exact output filenames evolve between releases.  We do NOT
hard-code names.  Instead we discover the structure by globbing:
  *.fasta / *.fa under the devider dir → per-window haplotype FASTAs
  *.bam under the devider dir          → haplotype-tagged BAM
Window coordinates are parsed defensively from FASTA record headers.

Header parsing
--------------
We try several known DEVIDER-style header formats:
  >WINDOW_START-WINDOW_END_HAPLOTYPE_INDEX
  >region:START-END|hap:N
  >hap_N_window_START_END
  >chrom:START-END_hapN
If a header cannot be parsed we still emit the sequence as a chain-of-one
(a fallback) and record it in the report under ``unparsed_headers``.

Haplotype-tag parsing
---------------------
For each alignment in the haplotype-tagged BAM we look (in order of
preference) at:
  * ``HP``  standard haplotype tag (integer-valued)
  * ``YH``  DEVIDER-specific tag (string-valued, ``window:hap`` format)
  * read-name suffix (``..._hapN`` / ``..._h:N``) — last-resort
Whichever is present is used; if none match the read is skipped for the
purposes of building junction evidence (it cannot tell us which haplotype
it belongs to).

Outputs
-------
1. ``merged_haplotypes.fasta``: one record per emitted chain.  Header format
   ``>{sample_id}_{genotype}_haplotype_{N}`` where ``N`` is the chain index
   (zero-based).  Chain index 0 is the highest-coverage / longest chain; ties
   are broken by abundance then by chain length.

2. ``stitch_report.json``: machine-readable record of every window, every
   junction, every link, and the final emitted chains.  Schema in the
   ``write_report`` function.

Dependencies: ``pysam``, Python standard library only (no NumPy / pandas).

See also
--------
* ``docs/data_flow.md``                Steps 5.19, 5.20
* ``docs/implementation_prompts.md``   Prompt 21
* ``docs/architecture_reasoning.md``   §9
* ``CLAUDE.md``                        D11, params.stitch_min_reads
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

# pysam is imported lazily inside `collect_window_assignments` so the
# graceful-failure fallback path (no DEVIDER outputs / devider.failed marker)
# still runs even if pysam is unavailable for some reason at runtime.
# The Nextflow container (LoFreq biocontainer) bundles pysam, so production
# invocations always have it.


# --------------------------------------------------------------------------- #
# Logging                                                                     #
# --------------------------------------------------------------------------- #

def log(msg: str) -> None:
    """Structured log line to stderr (matches the convention used elsewhere)."""
    print(f"[stitch_haplotypes] {msg}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# Data classes                                                                #
# --------------------------------------------------------------------------- #

@dataclass
class WindowHaplotype:
    """One haplotype sequence inside one DEVIDER window."""
    window_index: int           # ordinal: 0, 1, 2, ... after sorting by start
    window_start: int           # 1-based genomic start (inclusive), or 0 if unknown
    window_end: int             # 1-based genomic end (inclusive), or 0 if unknown
    hap_index: int              # haplotype ordinal within the window (0, 1, 2, ...)
    raw_id: str                 # original FASTA header (without leading '>')
    sequence: str               # nucleotide sequence (uppercase, gap-free)
    reported_abundance: Optional[float] = None  # DEVIDER-reported abundance fraction (0-1)

    @property
    def key(self) -> str:
        """Stable identifier used in the report (e.g. 'W0_hap1')."""
        return f"W{self.window_index}_hap{self.hap_index}"


@dataclass
class Link:
    """A read-supported link between adjacent windows."""
    upstream_window: int        # window_index of the earlier window
    upstream_hap: int           # hap_index in that window
    downstream_window: int      # window_index of the later (adjacent) window
    downstream_hap: int         # hap_index in that window
    reads: int                  # number of spanning reads
    supported: bool             # reads >= min_reads?

    def to_dict(self) -> dict[str, Any]:
        return {
            "from": f"W{self.upstream_window}_hap{self.upstream_hap}",
            "to":   f"W{self.downstream_window}_hap{self.downstream_hap}",
            "reads": self.reads,
            "supported": self.supported,
        }


@dataclass
class Junction:
    """All links observed at one adjacent-window boundary."""
    upstream_window: int
    downstream_window: int
    spanning_reads: int                 # total spanning reads across all link pairings
    links: list[Link] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "window_pair": [f"W{self.upstream_window}", f"W{self.downstream_window}"],
            "spanning_reads": self.spanning_reads,
            "links": [lk.to_dict() for lk in self.links],
        }


@dataclass
class Chain:
    """One emitted merged haplotype: a path through the windows."""
    windows: list[int]                  # window_index list (in order)
    hap_indices: list[int]              # hap_index per window (same length)
    sequence: str                       # concatenated sequence
    spanning_reads_per_junction: list[int]
    abundance_lower_bound: float        # min of (link_reads / total_spanning) across junctions


# --------------------------------------------------------------------------- #
# DEVIDER output discovery                                                    #
# --------------------------------------------------------------------------- #

# Header-parser regex patterns — tried in order.  Each must capture
# (start, end, hap_idx) as groups 1/2/3.  All accept underscores or hyphens.
_HEADER_PATTERNS = [
    # >hap0_window_123_4567   (DEVIDER 0.0.x typical: hap then window)
    re.compile(r"hap(?:lotype)?[_]?(?P<hap>\d+).*?window[_:](?P<start>\d+)[_:\-](?P<end>\d+)", re.IGNORECASE),
    # >chrom:123-4567_hap0  or  >region:123-4567|hap:0
    re.compile(r"[:_](?P<start>\d+)[-_](?P<end>\d+)[_|].*?hap(?:lotype)?[_:]?(?P<hap>\d+)", re.IGNORECASE),
    # >window_123_4567_hap0  /  >123-4567_hap_0
    re.compile(r"(?:window[_:])?(?P<start>\d+)[_\-](?P<end>\d+)[_\-]?hap(?:lotype)?[_:]?(?P<hap>\d+)", re.IGNORECASE),
    # >123-4567_0   (start-end_hap)
    re.compile(r"^(?P<start>\d+)[_\-](?P<end>\d+)[_\-](?P<hap>\d+)\b"),
]

# Matches DEVIDER 0.0.1 whole-genome (single-window) output:
#   >Contig:...,Range:ALL-ALL,Haplotype:N,...
_ALLALL_RE = re.compile(r"Range:ALL-ALL.*?Haplotype:(?P<hap>\d+)", re.IGNORECASE)
# Abundance field in DEVIDER ALL-ALL headers: Abundance:12.55 (percentage)
_ABUNDANCE_RE = re.compile(r"Abundance:(?P<ab>[\d.]+)", re.IGNORECASE)


def parse_devider_header(header: str) -> Optional[tuple[int, int, int]]:
    """
    Parse a DEVIDER haplotype FASTA header into (start, end, hap_index).

    DEVIDER 0.0.1 ALL-ALL mode (whole-genome single window) emits headers like:
        >Contig:...,Range:ALL-ALL,Haplotype:N,Abundance:X,Depth:Y
    All haplotypes are assigned to virtual window (0, 0) so they are grouped
    as a single window and emitted as individual chains without stitching.

    Returns None if no pattern matches; the caller falls back to ordinal-only
    handling (each unrecognised header becomes its own pseudo-window chain).
    """
    name = header.split(None, 1)[0]

    # DEVIDER 0.0.1 ALL-ALL: single virtual window at (0, 0).
    m = _ALLALL_RE.search(name)
    if m:
        try:
            return (0, 0, int(m.group("hap")))
        except (KeyError, ValueError):
            pass

    for pat in _HEADER_PATTERNS:
        m = pat.search(name)
        if m:
            try:
                return (int(m.group("start")), int(m.group("end")), int(m.group("hap")))
            except (KeyError, ValueError):
                continue
    return None


def _parse_devider_abundance(header: str) -> Optional[float]:
    """Extract DEVIDER-reported abundance percentage from header, as a fraction 0-1."""
    m = _ABUNDANCE_RE.search(header)
    if m:
        try:
            return float(m.group("ab")) / 100.0
        except ValueError:
            pass
    return None


def _iter_fasta_records(path: Path):
    """
    Tiny dependency-free FASTA reader.  Yields (header_str, sequence_str).
    Sequence is uppercased and stripped of whitespace.
    """
    header: Optional[str] = None
    seq_parts: list[str] = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts).upper()


def discover_devider_outputs(devider_dir: Path) -> tuple[list[Path], list[Path], bool]:
    """
    Discover DEVIDER outputs in ``devider_dir`` (recursively).

    Returns (fasta_paths, bam_paths, devider_failed).
    """
    if not devider_dir.exists() or not devider_dir.is_dir():
        return [], [], True  # treat absent dir as failure

    # Failure marker takes precedence over content.
    if (devider_dir / "devider.failed").exists():
        return [], [], True

    # Only look at top-level files — do NOT recurse into intermediate/.
    # DEVIDER 0.0.1 emits two FASTA files at the top level:
    #   majority_vote_haplotypes.fasta — full-length genome sequences (preferred)
    #   snp_haplotypes.fasta           — variable-site-only sequences (261 bp)
    # We want the full-genome sequences for the merged haplotype output.
    preferred = devider_dir / "majority_vote_haplotypes.fasta"
    if preferred.is_file():
        fastas = [preferred]
    else:
        fastas = sorted(
            list(devider_dir.glob("*.fasta")) + list(devider_dir.glob("*.fa"))
        )
    bams = sorted(devider_dir.glob("*.bam"))

    return fastas, bams, False


def load_haplotypes(fasta_paths: list[Path]) -> tuple[list[WindowHaplotype], list[str]]:
    """
    Parse all haplotype FASTAs and group them by (start, end) window.

    Returns:
        haplotypes — sorted by (window_index, hap_index).
        unparsed_headers — list of headers we could not parse coordinates from
                           (still emitted in the output FASTA as chains-of-one).
    """
    # Group records by (start, end).  Records whose header we cannot parse are
    # collected in a side-bucket and treated as their own pseudo-window.
    by_window: dict[tuple[int, int], list[tuple[int, str, str]]] = defaultdict(list)
    unparsed: list[tuple[str, str]] = []  # (raw_id, sequence)

    for fasta in fasta_paths:
        for header, sequence in _iter_fasta_records(fasta):
            if not sequence:
                continue
            parsed = parse_devider_header(header)
            if parsed is None:
                unparsed.append((header, sequence))
                continue
            start, end, hap = parsed
            abund = _parse_devider_abundance(header)
            by_window[(start, end)].append((hap, header, sequence, abund))

    # Assign deterministic ordinals: sort windows by start position, then end.
    sorted_windows = sorted(by_window.keys(), key=lambda se: (se[0], se[1]))
    haplotypes: list[WindowHaplotype] = []
    for w_idx, (start, end) in enumerate(sorted_windows):
        # Sort haplotypes within the window by their declared hap_index for stability.
        for hap, header, seq, abund in sorted(by_window[(start, end)], key=lambda x: x[0]):
            haplotypes.append(WindowHaplotype(
                window_index=w_idx,
                window_start=start,
                window_end=end,
                hap_index=hap,
                raw_id=header,
                sequence=seq,
                reported_abundance=abund,
            ))

    # Append unparsed records as their own pseudo-windows at the end.
    next_w_idx = len(sorted_windows)
    for header, seq in unparsed:
        haplotypes.append(WindowHaplotype(
            window_index=next_w_idx,
            window_start=0,
            window_end=0,
            hap_index=0,
            raw_id=header,
            sequence=seq,
        ))
        next_w_idx += 1

    unparsed_headers = [h for h, _ in unparsed]
    return haplotypes, unparsed_headers


# --------------------------------------------------------------------------- #
# BAM walking                                                                 #
# --------------------------------------------------------------------------- #

# Read-name suffix fallback: matches "_hap3" / "_h:3" at the end of a read name.
_NAME_HAP_SUFFIX = re.compile(r"_(?:h(?:ap)?[:_]?)(?P<hap>\d+)\s*$", re.IGNORECASE)


def extract_read_hap(read) -> Optional[int]:  # read: pysam.AlignedSegment
    """
    Determine the haplotype index of a tagged BAM read.

    Tag preference order:
        HP — standard haplotype tag (integer)
        YH — DEVIDER-specific tag.  Observed formats vary; we accept either
             an integer or a string of the form "window:hap" / "hapN".
    Falls back to a read-name suffix regex.

    Returns None when no tag is found — such reads cannot contribute to
    junction evidence and are skipped.
    """
    # HP (samtools / longshot convention) — always an integer.
    if read.has_tag("HP"):
        try:
            return int(read.get_tag("HP"))
        except (TypeError, ValueError):
            pass

    # YH (DEVIDER) — try integer, then parse as 'window:hap' / 'hapN'.
    if read.has_tag("YH"):
        val = read.get_tag("YH")
        if isinstance(val, int):
            return val
        if isinstance(val, str):
            # 'window:hap' format → take last numeric component.
            for piece in re.findall(r"\d+", val):
                try:
                    return int(piece)
                except ValueError:
                    continue

    # Read-name suffix fallback.
    if read.query_name:
        m = _NAME_HAP_SUFFIX.search(read.query_name)
        if m:
            try:
                return int(m.group("hap"))
            except ValueError:
                pass

    return None


def collect_window_assignments(
    bam_path: Path,
    windows: list[tuple[int, int]],
) -> dict[str, dict[int, int]]:
    """
    For each read in the BAM, determine which window(s) its alignment spans
    and record the haplotype tag.

    A read "spans" window W if its alignment ``[ref_start, ref_end)`` overlaps
    [W_start, W_end].  We use a permissive overlap because DEVIDER windows
    themselves usually overlap and we want to catch any read that could
    contribute to junction evidence.

    Returns:
        ``{read_name: {window_index: hap_index, ...}}``

    Only reads with a non-None haplotype tag are included; only reads that
    overlap at least one window are included.
    """
    if not bam_path.exists():
        return {}

    # Lazy import so the no-DEVIDER fallback path doesn't require pysam.
    try:
        import pysam  # noqa: PLC0415
    except ImportError as exc:
        log(f"pysam not available ({exc}); cannot walk BAM — no junction evidence")
        return {}

    assignments: dict[str, dict[int, int]] = defaultdict(dict)
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except (OSError, ValueError) as exc:
        log(f"could not open BAM {bam_path}: {exc}")
        return {}

    try:
        for read in bam.fetch(until_eof=True):
            # Skip unmapped, secondary, supplementary and duplicates — only
            # primary alignments carry the canonical haplotype tag.
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate:
                continue

            hap = extract_read_hap(read)
            if hap is None:
                continue

            # pysam reference coordinates are 0-based half-open; the window
            # coordinates we parsed from FASTA headers are typically 1-based
            # inclusive.  We treat both as inclusive integers for overlap-
            # checking and accept any positive overlap.
            r_start = read.reference_start + 1  # convert to 1-based inclusive
            r_end = read.reference_end          # already 1-based inclusive end

            if r_end is None:
                continue

            for w_idx, (w_start, w_end) in enumerate(windows):
                if w_end <= 0:
                    continue  # synthetic / unparsed window
                if r_start <= w_end and r_end >= w_start:
                    # Record the haplotype tag the read carries while it
                    # overlaps this window.  If a read appears in multiple
                    # primary records (it shouldn't), the last one wins.
                    assignments[read.query_name][w_idx] = hap

    finally:
        bam.close()

    return assignments


# --------------------------------------------------------------------------- #
# Junction evidence and greedy chain assembly                                 #
# --------------------------------------------------------------------------- #

def build_junctions(
    assignments: dict[str, dict[int, int]],
    n_windows: int,
    min_reads: int,
) -> list[Junction]:
    """
    For each adjacent window pair ``(W_i, W_{i+1})`` count the number of
    spanning reads that support each ``(hap_in_W_i, hap_in_W_{i+1})`` pair.
    Return one Junction per adjacent pair.
    """
    junctions: list[Junction] = []
    for i in range(n_windows - 1):
        j = i + 1
        # count[(hap_i, hap_j)] = number of spanning reads with that pairing
        pair_counts: dict[tuple[int, int], int] = defaultdict(int)
        total_spanning = 0
        for read_name, win_haps in assignments.items():
            if i in win_haps and j in win_haps:
                pair_counts[(win_haps[i], win_haps[j])] += 1
                total_spanning += 1

        links = [
            Link(
                upstream_window=i,
                upstream_hap=hap_i,
                downstream_window=j,
                downstream_hap=hap_j,
                reads=reads,
                supported=(reads >= min_reads),
            )
            for (hap_i, hap_j), reads in sorted(pair_counts.items())
        ]
        # Stable, useful sort: most-supported first.
        links.sort(key=lambda lk: (-lk.reads, lk.upstream_hap, lk.downstream_hap))
        junctions.append(Junction(
            upstream_window=i,
            downstream_window=j,
            spanning_reads=total_spanning,
            links=links,
        ))
    return junctions


def assemble_chains(
    haplotypes: list[WindowHaplotype],
    junctions: list[Junction],
    n_windows: int,
) -> list[Chain]:
    """
    Greedy left-to-right chain assembly.

    Strategy:
      * For each haplotype in window 0, start a chain.
      * At each junction, follow ALL supported links from the current
        upstream haplotype.  Multiple supported branches → emit multiple
        chains (we enumerate, not arbitrate).
      * A chain ends at a window where no supported link exits.
      * Haplotypes never visited as part of any chain are emitted as
        chains-of-one so that no DEVIDER haplotype is silently dropped.

    Abundance (lower bound):
      For a chain that traverses junctions J_a, J_b, J_c with link reads
      r_a, r_b, r_c and total spanning reads T_a, T_b, T_c, the abundance
      is ``min(r_a/T_a, r_b/T_b, r_c/T_c)``.  This is a lower bound
      because any one junction may understate the true population
      proportion of the haplotype.
    """
    # Index haplotypes by (window_index, hap_index) → WindowHaplotype.
    hap_by_pos: dict[tuple[int, int], WindowHaplotype] = {
        (h.window_index, h.hap_index): h for h in haplotypes
    }

    # Index junctions by upstream window.
    junctions_by_upstream: dict[int, Junction] = {j.upstream_window: j for j in junctions}

    # Build adjacency: for each (window, hap), the list of supported downstream haps.
    next_haps: dict[tuple[int, int], list[tuple[int, int, int]]] = defaultdict(list)
    # value tuples: (next_hap_index, link_reads, total_spanning_at_junction)

    for j in junctions:
        for lk in j.links:
            if not lk.supported:
                continue
            next_haps[(lk.upstream_window, lk.upstream_hap)].append(
                (lk.downstream_hap, lk.reads, j.spanning_reads)
            )

    # Walk all chains starting from window 0 haplotypes.
    chains: list[Chain] = []
    visited: set[tuple[int, int]] = set()

    def walk(start_window: int, start_hap: int):
        """Depth-first enumeration of all greedy paths from a starting node."""
        # Each frame: (path_windows, path_haps, per_junction_reads, per_junction_total)
        stack = [([start_window], [start_hap], [], [])]
        while stack:
            windows, haps, link_reads, link_totals = stack.pop()
            cur_w, cur_h = windows[-1], haps[-1]
            visited.add((cur_w, cur_h))
            successors = next_haps.get((cur_w, cur_h), [])
            if not successors or cur_w + 1 >= n_windows:
                # Terminal node — emit chain.
                chains.append(_build_chain(
                    windows, haps, link_reads, link_totals, hap_by_pos
                ))
                continue
            for nh, reads, total in successors:
                stack.append((
                    windows + [cur_w + 1],
                    haps + [nh],
                    link_reads + [reads],
                    link_totals + [total],
                ))

    window_zero_haps = sorted(
        [h.hap_index for h in haplotypes if h.window_index == 0]
    )
    for h in window_zero_haps:
        walk(0, h)

    # Any haplotype not visited by a chain starting at window 0 gets emitted
    # as its own single-window chain.  This preserves DEVIDER's native output
    # for windows that are not reachable from window 0 via supported links
    # (e.g. orphan windows past a junction with no spanning reads).
    for h in haplotypes:
        if (h.window_index, h.hap_index) in visited:
            continue
        # For unvisited haplotypes, try walking forward from them too, so
        # we capture multi-window chains that don't start at window 0.
        if (h.window_index, h.hap_index) not in visited:
            walk(h.window_index, h.hap_index)

    # Dedup: chains can be re-emitted by the multi-start walk above.
    deduped: dict[tuple, Chain] = {}
    for c in chains:
        key = (tuple(c.windows), tuple(c.hap_indices))
        if key not in deduped:
            deduped[key] = c
    chains = list(deduped.values())

    # Sort: longest first, then highest abundance, then lowest leading window.
    chains.sort(
        key=lambda c: (-len(c.windows), -c.abundance_lower_bound, c.windows[0])
    )
    return chains


def _build_chain(
    windows: list[int],
    haps: list[int],
    link_reads: list[int],
    link_totals: list[int],
    hap_by_pos: dict[tuple[int, int], WindowHaplotype],
) -> Chain:
    """Concatenate sequences and compute the abundance lower bound."""
    sequence = _concat_sequences(windows, haps, hap_by_pos)
    if link_totals:
        # min(link_reads / link_totals) — lower bound on abundance.
        fractions = [
            (reads / total) if total > 0 else 0.0
            for reads, total in zip(link_reads, link_totals)
        ]
        abundance = min(fractions) if fractions else 1.0
    else:
        # Single-window chain — use DEVIDER-reported abundance if available.
        reported = hap_by_pos.get((windows[0], haps[0]))
        if reported is not None and reported.reported_abundance is not None:
            abundance = reported.reported_abundance
        else:
            abundance = 1.0
    return Chain(
        windows=list(windows),
        hap_indices=list(haps),
        sequence=sequence,
        spanning_reads_per_junction=list(link_reads),
        abundance_lower_bound=abundance,
    )


def _concat_sequences(
    windows: list[int],
    haps: list[int],
    hap_by_pos: dict[tuple[int, int], WindowHaplotype],
) -> str:
    """
    Concatenate sequences for the windows in the chain.

    Boundary handling at adjacent windows:
      * If windows do not overlap (w_end_i < w_start_{i+1}), there is a gap
        we cannot fill from DEVIDER alone — we emit the upstream sequence
        verbatim and append the downstream sequence (no gap-filling here;
        the consensus-based gap fill described in the prompt body is
        deferred to downstream analysis to keep this script's contract
        minimal and to avoid silently injecting consensus bases into
        haplotypes).  We log a warning in this case.
      * If windows overlap (w_end_i >= w_start_{i+1}), we trim the
        downstream sequence at the midpoint of the overlap to avoid
        emitting duplicated bases at the junction.

    For chains with only one window, we just return its sequence as-is.
    """
    if not windows:
        return ""
    parts: list[str] = []
    prev: Optional[WindowHaplotype] = None
    for w_idx, h_idx in zip(windows, haps):
        cur = hap_by_pos[(w_idx, h_idx)]
        if prev is None:
            parts.append(cur.sequence)
        else:
            # Compute overlap in genomic coordinates.
            overlap = prev.window_end - cur.window_start + 1
            if overlap <= 0:
                # Non-overlapping → just append.
                if overlap < 0:
                    log(
                        f"non-overlapping windows W{prev.window_index} "
                        f"({prev.window_start}-{prev.window_end}) → "
                        f"W{cur.window_index} ({cur.window_start}-{cur.window_end}); "
                        f"emitting without gap fill"
                    )
                parts.append(cur.sequence)
            else:
                # Trim downstream sequence by half the overlap; trim upstream
                # by the other half (split overlap at the midpoint).
                half = overlap // 2
                # Trim the tail of the already-appended upstream part.
                if half > 0 and len(parts[-1]) > half:
                    parts[-1] = parts[-1][:-half]
                trim_downstream = overlap - half
                if trim_downstream >= len(cur.sequence):
                    # Pathological: downstream window is entirely inside
                    # upstream window.  Skip it.
                    log(
                        f"window W{cur.window_index} entirely inside "
                        f"W{prev.window_index}; skipping in concatenation"
                    )
                else:
                    parts.append(cur.sequence[trim_downstream:])
        prev = cur
    return "".join(parts)


# --------------------------------------------------------------------------- #
# Output                                                                      #
# --------------------------------------------------------------------------- #

def write_fasta(
    out_path: Path,
    chains: list[Chain],
    sample_id: str,
    genotype: str,
) -> None:
    """
    Write one record per chain.

    Header: ``>{sample_id}_{genotype}_haplotype_{N}`` (N = chain ordinal).
    The sequence is written in 80-column lines for compatibility with
    downstream FASTA-consuming tools.
    """
    with open(out_path, "w") as fh:
        for n, chain in enumerate(chains):
            fh.write(f">{sample_id}_{genotype}_haplotype_{n}\n")
            seq = chain.sequence
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80])
                fh.write("\n")


def build_report(
    sample_id: str,
    genotype: str,
    haplotypes: list[WindowHaplotype],
    junctions: list[Junction],
    chains: list[Chain],
    unparsed_headers: list[str],
    fallback_used: bool,
    fallback_reason: Optional[str],
    bam_present: bool,
    min_reads: int,
) -> dict[str, Any]:
    """Assemble the stitch report JSON payload."""
    # Per-window haplotype counts (ordered by window_index).
    counts_by_window: dict[int, int] = defaultdict(int)
    for h in haplotypes:
        counts_by_window[h.window_index] += 1
    n_windows = len(counts_by_window)
    haplotypes_per_window = [counts_by_window[i] for i in sorted(counts_by_window)]

    # Identify windows whose downstream junction has zero supported links.
    unlinked: list[str] = []
    for j in junctions:
        if not any(lk.supported for lk in j.links):
            unlinked.append(f"W{j.upstream_window}-W{j.downstream_window}")

    chains_dict = []
    for n, c in enumerate(chains):
        chains_dict.append({
            "id": f"{sample_id}_{genotype}_haplotype_{n}",
            "windows": c.windows,
            "hap_indices": c.hap_indices,
            "length_bp": len(c.sequence),
            "spanning_reads": c.spanning_reads_per_junction,
            "abundance_lower_bound": round(c.abundance_lower_bound, 4),
        })

    report: dict[str, Any] = {
        "sample_id": sample_id,
        "genotype": genotype,
        "min_reads_threshold": min_reads,
        "bam_present": bam_present,
        "windows_found": n_windows,
        "haplotypes_per_window": haplotypes_per_window,
        "n_input_haplotypes": len(haplotypes),
        "junctions": [j.to_dict() for j in junctions],
        "chains_stitched": len(chains),
        "chains": chains_dict,
        "unlinked_windows": unlinked,
        "unparsed_headers": unparsed_headers,
        "fallback_used": fallback_used,
    }
    if fallback_reason:
        report["reason"] = fallback_reason
    return report


def write_report(out_path: Path, report: dict[str, Any]) -> None:
    with open(out_path, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=False)
        fh.write("\n")


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Post-hoc stitching of DEVIDER per-window haplotypes using "
            "read-spanning evidence from a haplotype-tagged BAM."
        ),
    )
    p.add_argument("--devider-dir", required=True, type=Path,
                   help="DEVIDER output directory (will be globbed for FASTAs and BAM).")
    p.add_argument("--output-fasta", required=True, type=Path,
                   help="Output FASTA of stitched haplotypes.")
    p.add_argument("--output-json", required=True, type=Path,
                   help="Output JSON stitching report.")
    p.add_argument("--sample-id", required=True,
                   help="Sample identifier (used in record headers and report).")
    p.add_argument("--genotype", required=True,
                   help="Genotype label (used in record headers and report).")
    p.add_argument("--min-reads", type=int, default=5,
                   help="Minimum spanning reads to support a window-to-window link "
                        "(matches params.stitch_min_reads; default 5).")
    return p.parse_args(argv)


def _emit_empty(args: argparse.Namespace, reason: str) -> None:
    """Emit empty FASTA + JSON when DEVIDER produced nothing usable."""
    log(reason)
    # Empty FASTA — keep the file present so Nextflow output globs resolve.
    args.output_fasta.write_text("")
    report = build_report(
        sample_id=args.sample_id,
        genotype=args.genotype,
        haplotypes=[],
        junctions=[],
        chains=[],
        unparsed_headers=[],
        fallback_used=True,
        fallback_reason=reason,
        bam_present=False,
        min_reads=args.min_reads,
    )
    write_report(args.output_json, report)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)

    # --- Phase 0: discover DEVIDER outputs -------------------------------- #
    fasta_paths, bam_paths, devider_failed = discover_devider_outputs(args.devider_dir)
    if devider_failed:
        _emit_empty(args, "DEVIDER failed (devider.failed marker or missing dir)")
        return 0
    if not fasta_paths:
        _emit_empty(args, "no DEVIDER haplotype FASTAs found")
        return 0

    log(f"discovered {len(fasta_paths)} FASTA file(s), {len(bam_paths)} BAM file(s)")

    # --- Phase 1: parse haplotypes ---------------------------------------- #
    haplotypes, unparsed = load_haplotypes(fasta_paths)
    if not haplotypes:
        _emit_empty(args, "no parseable haplotype sequences in DEVIDER output")
        return 0
    n_windows = len({h.window_index for h in haplotypes})
    log(f"parsed {len(haplotypes)} haplotypes across {n_windows} windows "
        f"({len(unparsed)} unparsed headers)")

    # --- Phase 2: walk haplotype-tagged BAM ------------------------------- #
    # Window coordinates list, indexed by window_index.
    window_coords: dict[int, tuple[int, int]] = {}
    for h in haplotypes:
        if h.window_index not in window_coords:
            window_coords[h.window_index] = (h.window_start, h.window_end)
    n_total_windows = max(window_coords) + 1 if window_coords else 0
    coord_list = [window_coords.get(i, (0, 0)) for i in range(n_total_windows)]

    # If we have a BAM, walk it; otherwise empty assignments (no links).
    bam_present = bool(bam_paths)
    assignments: dict[str, dict[int, int]] = {}
    if bam_present:
        # Use the first BAM found.  DEVIDER v0.0.1 emits exactly one tagged
        # BAM per run; if multiple were found we picked the lexicographically
        # first.  No need to merge — they would carry the same information.
        bam_path = bam_paths[0]
        log(f"using haplotype-tagged BAM: {bam_path}")
        try:
            assignments = collect_window_assignments(bam_path, coord_list)
        except Exception as exc:  # noqa: BLE001 — never fail the sample
            log(f"BAM walking failed ({exc!r}); proceeding without junction evidence")
            assignments = {}
    else:
        log("no haplotype-tagged BAM found; emitting haplotypes as-is "
            "(no junction evidence — every window haplotype becomes its own chain)")

    # --- Phase 3: build junctions and chains ------------------------------ #
    junctions = build_junctions(assignments, n_total_windows, args.min_reads)
    log(f"built {len(junctions)} junction(s); "
        f"{sum(1 for j in junctions for lk in j.links if lk.supported)} supported link(s)")
    chains = assemble_chains(haplotypes, junctions, n_total_windows)
    log(f"assembled {len(chains)} chain(s)")

    # --- Phase 4: write outputs ------------------------------------------- #
    write_fasta(args.output_fasta, chains, args.sample_id, args.genotype)
    report = build_report(
        sample_id=args.sample_id,
        genotype=args.genotype,
        haplotypes=haplotypes,
        junctions=junctions,
        chains=chains,
        unparsed_headers=unparsed,
        fallback_used=(not bam_present),
        fallback_reason=(None if bam_present else "no haplotype-tagged BAM available"),
        bam_present=bam_present,
        min_reads=args.min_reads,
    )
    write_report(args.output_json, report)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 — graceful degradation contract
        # Last-resort guard so we never fail the sample on an unhandled error.
        log(f"unhandled exception: {exc!r}; emitting empty outputs")
        try:
            # Try to write empty outputs if at all possible.
            args = parse_args()
            args.output_fasta.write_text("")
            write_report(
                args.output_json,
                {
                    "sample_id": args.sample_id,
                    "genotype": args.genotype,
                    "fallback_used": True,
                    "reason": f"unhandled exception: {exc!r}",
                },
            )
        except Exception:
            pass
        sys.exit(0)
