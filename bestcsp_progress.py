#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create BEST-CSP progress tables and per-form horizontal progress bars.

The repository is scanned recursively. Every immediate subdirectory of ROOT is
interpreted as one compound. A separate progress bar is produced for every
crystalline form.

The absolute score is:

    score = laboratories_150K + laboratories_room_temperature
            + thermodynamic_data_records

Default target counts are:

* 3 laboratories with crystal-structure measurements at 150 ± 10 K;
* 3 laboratories with crystal-structure measurements near room temperature;
* 5 thermodynamic data records.

Thermodynamic records can include phase-transition temperatures and phase-
transition enthalpies, including melting/fusion and sublimation data. The score
is not capped at the target line and bars may extend beyond it.

The CSV ``Name`` field is used as a laboratory/contributor key after automatic
normalisation. Use ``lab_aliases.csv`` to merge labels belonging to the same
laboratory.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator


# ---------------------------------------------------------------------------
# Defaults that can also be changed through command-line arguments
# ---------------------------------------------------------------------------
DEFAULT_TARGET_150_K = 150.0
DEFAULT_TOLERANCE_150_K = 10.0
DEFAULT_ROOM_TEMPERATURE_K = 298.15
DEFAULT_ROOM_TOLERANCE_K = 10.0
DEFAULT_REQUIRED_150_LABS = 3
DEFAULT_REQUIRED_ROOM_LABS = 3
DEFAULT_REQUIRED_THERMODYNAMIC_DATA = 5
DEFAULT_MINIMUM_VISIBLE_BAR = 1.0

IGNORED_DIRECTORY_NAMES = {
    ".git",
    ".github",
    "__pycache__",
    ".venv",
    "venv",
    "env",
    "progress",
}

STRUCTURE_KEYWORDS = (
    "lattice",
    "unit cell",
    "unit_cell",
    "cell volume",
    "cell_volume",
    "angle alpha",
    "angle beta",
    "angle gamma",
    "thermal expansion",
    "thermal_expansion",
)

MELTING_TEMPERATURE_KEYWORDS = (
    "melting point",
    "melting_point",
)

TRANSITION_TEMPERATURE_KEYWORDS = (
    "temperature of transition",
    "temperature_of_transition",
    "transition temperature",
    "transition_temperature",
    "transition point",
    "transition_point",
    "dsc transition",
    "dsc_transition",
)

THERMODYNAMIC_ENTHALPY_KEYWORDS = (
    "enthalpy of transition",
    "enthalpy_of_transition",
    "transition enthalpy",
    "transition_enthalpy",
    "enthalpy of fusion",
    "enthalpy_of_fusion",
    "enthalpy fusion",
    "enthalpy_fusion",
    "enthalpy melting",
    "enthalpy_melting",
    "enthalpy of melting",
    "enthalpy_of_melting",
    "enthalpy of sublimation",
    "enthalpy_of_sublimation",
    "enthalpy sublimation",
    "enthalpy_sublimation",
)

FORM_WORD_TO_LABEL = {
    "one": "I",
    "two": "II",
    "three": "III",
    "four": "IV",
    "five": "V",
    "six": "VI",
    "seven": "VII",
    "eight": "VIII",
    "nine": "IX",
    "ten": "X",
}

GREEK_FORM_LABELS = {
    "alpha": "α",
    "beta": "β",
    "gamma": "γ",
    "delta": "δ",
    "epsilon": "ε",
}

ROMAN_ORDER = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
    "V": 5,
    "VI": 6,
    "VII": 7,
    "VIII": 8,
    "IX": 9,
    "X": 10,
}

GREEK_ORDER = {"α": 101, "β": 102, "γ": 103, "δ": 104, "ε": 105}


@dataclass(frozen=True)
class Evidence:
    compound: str
    form: str
    category: str
    temperature_k: float | None
    lab_raw: str
    lab_normalized: str
    source_file: str
    source_row: int
    property_text: str
    reference: str
    comment: str

    @property
    def record_id(self) -> str:
        """Stable identifier for one accepted CSV data row."""
        return f"{self.source_file}#row={self.source_row}"


@dataclass(frozen=True)
class SummaryRow:
    compound: str
    form: str
    labs_150k: tuple[str, ...]
    labs_room: tuple[str, ...]
    thermodynamic_records: tuple[str, ...]
    thermodynamic_temperature_records: tuple[str, ...]
    thermodynamic_enthalpy_records: tuple[str, ...]
    thermodynamic_labs: tuple[str, ...]

    @property
    def n_150k(self) -> int:
        return len(self.labs_150k)

    @property
    def n_room(self) -> int:
        return len(self.labs_room)

    @property
    def n_thermodynamic(self) -> int:
        return len(self.thermodynamic_records)

    @property
    def n_thermodynamic_temperature(self) -> int:
        return len(self.thermodynamic_temperature_records)

    @property
    def n_thermodynamic_enthalpy(self) -> int:
        return len(self.thermodynamic_enthalpy_records)

    @property
    def score(self) -> int:
        return self.n_150k + self.n_room + self.n_thermodynamic


# ---------------------------------------------------------------------------
# Text and CSV utilities
# ---------------------------------------------------------------------------
def read_text(path: Path) -> str:
    """Read a text file using common encodings found in the repository."""
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def read_csv_rows(path: Path) -> list[list[str]]:
    """Read a CSV file permissively, including malformed header variants."""
    text = read_text(path)
    return [row for row in csv.reader(text.splitlines()) if any(cell.strip() for cell in row)]


def clean_cell(value: str) -> str:
    value = value.replace("\u00a0", " ").strip().strip('"').strip()
    return re.sub(r"\s+", " ", value)


def first_nonempty_cell(row: Sequence[str]) -> str:
    for cell in row:
        value = clean_cell(cell)
        if value:
            return value
    return ""


def looks_like_data_row(row: Sequence[str]) -> bool:
    if len(row) < 6:
        return False
    first = clean_cell(row[0])
    return bool(re.fullmatch(r"#?\d+(?:\.\d+)?", first))


def standard_data_rows(path: Path) -> Iterable[tuple[int, str, str, str, str]]:
    """Yield row number, property, lab key, reference and comment.

    Most BEST-CSP files use positional columns:
    Identifier, Property, Value, Std, N, Name, Reference, Comment.
    """
    for row_number, row in enumerate(read_csv_rows(path), start=1):
        if not looks_like_data_row(row):
            continue
        padded = list(row) + [""] * max(0, 8 - len(row))
        property_text = clean_cell(padded[1])
        lab = clean_cell(padded[5])
        reference = clean_cell(padded[6])
        comment = clean_cell(", ".join(padded[7:]))
        if not lab:
            lab = reference
        if property_text and lab:
            yield row_number, property_text, lab, reference, comment


def metadata_identifier(path: Path) -> str:
    rows = read_csv_rows(path)
    if not rows:
        return ""
    return first_nonempty_cell(rows[0])


# ---------------------------------------------------------------------------
# Form detection
# ---------------------------------------------------------------------------
def normalise_form_token(token: str) -> str:
    token = token.strip().replace("_", "-")
    low = token.lower()
    if low in FORM_WORD_TO_LABEL:
        return FORM_WORD_TO_LABEL[low]
    if low in GREEK_FORM_LABELS:
        return GREEK_FORM_LABELS[low]
    if re.fullmatch(r"[ivx]+", low):
        return low.upper()
    if re.fullmatch(r"[a-z]", low):
        return low.upper()
    return token


def forms_from_text(text: str) -> tuple[str, ...]:
    """Extract one or more crystalline-form labels from filename or metadata."""
    text_clean = text.replace("–", "-").replace("—", "-")

    match = re.search(
        r"(?i)(?:^|[^A-Za-z0-9])form[\s_-]+"
        r"(one|two|three|four|five|six|seven|eight|nine|ten|[ivx]+|[a-z])"
        r"(?:\s*-\s*([a-z]))?(?=$|[^A-Za-z0-9])",
        text_clean,
    )
    if match:
        first = normalise_form_token(match.group(1))
        second = match.group(2)
        if second:
            return (first, normalise_form_token(second))
        return (first,)

    low = text_clean.lower()
    found: list[str] = []
    for word, label in GREEK_FORM_LABELS.items():
        if re.search(rf"(?:^|[^a-z]){word}(?:$|[^a-z])", low):
            found.append(label)
    return tuple(found)


def build_identifier_form_map(csv_files: Sequence[Path]) -> Mapping[str, set[str]]:
    mapping: DefaultDict[str, set[str]] = defaultdict(set)
    for path in csv_files:
        forms = set(forms_from_text(path.stem))
        if not forms:
            forms.update(forms_from_text(metadata_identifier(path)))
        identifier = metadata_identifier(path).strip().lower()
        if identifier and forms:
            mapping[identifier].update(forms)
    return mapping


def detect_forms(path: Path, identifier_form_map: Mapping[str, set[str]]) -> tuple[str, ...]:
    forms = set(forms_from_text(path.stem))
    metadata = metadata_identifier(path)
    forms.update(forms_from_text(metadata))
    if not forms and metadata:
        forms.update(identifier_form_map.get(metadata.strip().lower(), set()))
    return tuple(sorted(forms, key=form_sort_key))


def form_sort_key(label: str) -> tuple[int, str]:
    if label in ROMAN_ORDER:
        return ROMAN_ORDER[label], label
    if label in GREEK_ORDER:
        return GREEK_ORDER[label], label
    if re.fullmatch(r"[A-Z]", label):
        return 200 + ord(label), label
    return 1000, label


# ---------------------------------------------------------------------------
# Classification and laboratory normalisation
# ---------------------------------------------------------------------------
def contains_any(text: str, keywords: Sequence[str]) -> bool:
    low = text.lower()
    return any(keyword in low for keyword in keywords)


def classify_measurement(
    path: Path,
    property_text: str,
    *,
    include_melting: bool,
    include_glass_transition: bool,
    include_sublimation_enthalpy: bool,
) -> str | None:
    """Classify a CSV row into a score category."""
    combined = f"{path.stem} {property_text}".lower().replace("-", "_")

    if contains_any(combined, STRUCTURE_KEYWORDS):
        return "structure"

    if "glass transition" in combined or "glass_transition" in combined:
        return "thermodynamic_temperature" if include_glass_transition else None

    if contains_any(combined, THERMODYNAMIC_ENTHALPY_KEYWORDS):
        if not include_melting and ("fusion" in combined or "melting" in combined):
            return None
        if not include_sublimation_enthalpy and "sublimation" in combined:
            return None
        return "thermodynamic_enthalpy"

    if contains_any(combined, TRANSITION_TEMPERATURE_KEYWORDS):
        return "thermodynamic_temperature"

    if include_melting and contains_any(combined, MELTING_TEMPERATURE_KEYWORDS):
        return "thermodynamic_temperature"

    return None


def extract_measurement_temperature(path: Path, property_text: str) -> float | None:
    combined = f"{property_text} {path.stem}"
    match = re.search(r"@\s*(-?\d+(?:\.\d+)?)\s*K\b", combined, flags=re.IGNORECASE)
    if not match:
        match = re.search(r"\b(-?\d+(?:\.\d+)?)\s*K\b", combined, flags=re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    if re.search(r"(?:^|[_\s-])RT(?:$|[_\s-])", path.stem, flags=re.IGNORECASE):
        return DEFAULT_ROOM_TEMPERATURE_K
    return None


def automatic_lab_normalisation(raw_lab: str) -> str:
    """Normalise repeated-run labels while preserving contributor identity."""
    lab = clean_cell(raw_lab)
    lab = lab.replace("-", "_")

    lab = re.sub(r"_?\d+(?:\.\d+)?\s*K\s*/?\s*min$", "", lab, flags=re.IGNORECASE)
    lab = re.sub(r"_?\d+(?:\.\d+)?Kmin$", "", lab, flags=re.IGNORECASE)
    lab = re.sub(r"(?<=[A-Za-z])\d+$", "", lab)
    lab = re.sub(r"_(?:19|20)\d{2}(?:\.\d+)?$", "", lab)
    lab = re.sub(r"_+", "_", lab).strip("_ ")
    lab = lab.replace("_", " ")
    return re.sub(r"\s+", " ", lab).strip()


def load_aliases(path: Path | None) -> dict[str, str]:
    """Load alias -> canonical laboratory mappings, case-insensitively."""
    aliases: dict[str, str] = {}
    if path is None or not path.exists():
        return aliases

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"alias", "lab"}
        if not reader.fieldnames or not required.issubset(
            {name.strip() for name in reader.fieldnames}
        ):
            raise ValueError(f"{path}: expected CSV columns 'alias' and 'lab'")
        for row in reader:
            alias = clean_cell(row.get("alias", ""))
            lab = clean_cell(row.get("lab", ""))
            if alias and lab and not alias.startswith("#"):
                aliases[alias.casefold()] = lab
    return aliases


def normalise_lab(raw_lab: str, aliases: Mapping[str, str]) -> str:
    raw = clean_cell(raw_lab)
    auto = automatic_lab_normalisation(raw)
    return aliases.get(raw.casefold(), aliases.get(auto.casefold(), auto))


def status_for_score(score: int, target_total: int) -> str:
    """Return a colour class; percentages are used only internally."""
    if target_total <= 0:
        return "green"
    ratio = score / target_total
    if ratio < 0.25:
        return "red"
    if ratio <= 0.75:
        return "yellow"
    return "green"


# ---------------------------------------------------------------------------
# Repository scan and aggregation
# ---------------------------------------------------------------------------
def immediate_compound_directories(root: Path) -> list[Path]:
    return sorted(
        [
            path
            for path in root.iterdir()
            if path.is_dir()
            and path.name not in IGNORED_DIRECTORY_NAMES
            and not path.name.startswith(".")
        ],
        key=lambda path: path.name.casefold(),
    )


def should_ignore(path: Path, root: Path, include_supplementary: bool) -> bool:
    relative_parts = path.relative_to(root).parts
    if any(part in IGNORED_DIRECTORY_NAMES or part.startswith(".") for part in relative_parts):
        return True
    if not include_supplementary and any(
        part.casefold() == "supplementary_data" for part in relative_parts
    ):
        return True
    return False


def scan_repository(
    root: Path,
    *,
    aliases: Mapping[str, str],
    target_150_k: float,
    tolerance_150_k: float,
    room_temperature_k: float,
    room_tolerance_k: float,
    include_supplementary: bool,
    include_melting: bool,
    include_glass_transition: bool,
    include_sublimation_enthalpy: bool,
) -> tuple[list[Evidence], Mapping[str, set[str]]]:
    evidence: list[Evidence] = []
    discovered_forms: DefaultDict[str, set[str]] = defaultdict(set)

    for compound_dir in immediate_compound_directories(root):
        compound = compound_dir.name
        csv_files = sorted(
            path
            for path in compound_dir.rglob("*.csv")
            if not should_ignore(path, root, include_supplementary)
        )
        identifier_form_map = build_identifier_form_map(csv_files)

        for path in csv_files:
            forms = detect_forms(path, identifier_form_map)
            if forms:
                discovered_forms[compound].update(forms)

            for row_number, property_text, raw_lab, reference, comment in standard_data_rows(path):
                measurement_type = classify_measurement(
                    path,
                    property_text,
                    include_melting=include_melting,
                    include_glass_transition=include_glass_transition,
                    include_sublimation_enthalpy=include_sublimation_enthalpy,
                )
                if measurement_type is None or not forms:
                    continue

                normalised_lab = normalise_lab(raw_lab, aliases)
                if not normalised_lab:
                    continue

                if measurement_type == "structure":
                    temperature = extract_measurement_temperature(path, property_text)
                    if temperature is None:
                        continue
                    if math.isclose(temperature, target_150_k, abs_tol=tolerance_150_k):
                        category = "structure_150K"
                    elif math.isclose(
                        temperature, room_temperature_k, abs_tol=room_tolerance_k
                    ):
                        category = "structure_RT"
                    else:
                        continue
                else:
                    temperature = extract_measurement_temperature(path, property_text)
                    category = measurement_type

                for form in forms:
                    evidence.append(
                        Evidence(
                            compound=compound,
                            form=form,
                            category=category,
                            temperature_k=temperature,
                            lab_raw=raw_lab,
                            lab_normalized=normalised_lab,
                            source_file=str(path.relative_to(root)),
                            source_row=row_number,
                            property_text=property_text,
                            reference=reference,
                            comment=comment,
                        )
                    )

    return evidence, discovered_forms


def build_summary(
    evidence: Sequence[Evidence], discovered_forms: Mapping[str, set[str]]
) -> list[SummaryRow]:
    structure_labs: DefaultDict[tuple[str, str, str], set[str]] = defaultdict(set)
    thermo_records: DefaultDict[tuple[str, str], set[str]] = defaultdict(set)
    thermo_temperature_records: DefaultDict[tuple[str, str], set[str]] = defaultdict(set)
    thermo_enthalpy_records: DefaultDict[tuple[str, str], set[str]] = defaultdict(set)
    thermo_labs: DefaultDict[tuple[str, str], set[str]] = defaultdict(set)

    for item in evidence:
        key = (item.compound, item.form)
        if item.category in {"structure_150K", "structure_RT"}:
            structure_labs[(item.compound, item.form, item.category)].add(item.lab_normalized)
        elif item.category in {"thermodynamic_temperature", "thermodynamic_enthalpy"}:
            thermo_records[key].add(item.record_id)
            thermo_labs[key].add(item.lab_normalized)
            if item.category == "thermodynamic_temperature":
                thermo_temperature_records[key].add(item.record_id)
            else:
                thermo_enthalpy_records[key].add(item.record_id)

    rows: list[SummaryRow] = []
    for compound in sorted(discovered_forms, key=str.casefold):
        for form in sorted(discovered_forms[compound], key=form_sort_key):
            key = (compound, form)
            rows.append(
                SummaryRow(
                    compound=compound,
                    form=form,
                    labs_150k=tuple(
                        sorted(
                            structure_labs[(compound, form, "structure_150K")],
                            key=str.casefold,
                        )
                    ),
                    labs_room=tuple(
                        sorted(
                            structure_labs[(compound, form, "structure_RT")],
                            key=str.casefold,
                        )
                    ),
                    thermodynamic_records=tuple(sorted(thermo_records[key], key=str.casefold)),
                    thermodynamic_temperature_records=tuple(
                        sorted(thermo_temperature_records[key], key=str.casefold)
                    ),
                    thermodynamic_enthalpy_records=tuple(
                        sorted(thermo_enthalpy_records[key], key=str.casefold)
                    ),
                    thermodynamic_labs=tuple(sorted(thermo_labs[key], key=str.casefold)),
                )
            )
    return rows


# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------
def write_evidence_csv(path: Path, evidence: Sequence[Evidence]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "compound",
        "form",
        "category",
        "temperature_K",
        "lab_raw",
        "lab_normalized",
        "source_file",
        "source_row",
        "property",
        "reference",
        "comment",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in sorted(
            evidence,
            key=lambda e: (
                e.compound.casefold(),
                form_sort_key(e.form),
                e.category,
                e.lab_normalized.casefold(),
                e.source_file.casefold(),
                e.source_row,
            ),
        ):
            writer.writerow(
                {
                    "compound": item.compound,
                    "form": item.form,
                    "category": item.category,
                    "temperature_K": ""
                    if item.temperature_k is None
                    else f"{item.temperature_k:g}",
                    "lab_raw": item.lab_raw,
                    "lab_normalized": item.lab_normalized,
                    "source_file": item.source_file,
                    "source_row": item.source_row,
                    "property": item.property_text,
                    "reference": item.reference,
                    "comment": item.comment,
                }
            )


def write_summary_csv(path: Path, summary: Sequence[SummaryRow], target_total: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "compound",
        "form",
        "labs_150K",
        "labs_RT",
        "thermodynamic_data",
        "thermodynamic_temperatures",
        "thermodynamic_enthalpies",
        "score",
        "status",
        "lab_names_150K",
        "lab_names_RT",
        "thermodynamic_lab_names",
        "thermodynamic_record_ids",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary:
            writer.writerow(
                {
                    "compound": row.compound,
                    "form": row.form,
                    "labs_150K": row.n_150k,
                    "labs_RT": row.n_room,
                    "thermodynamic_data": row.n_thermodynamic,
                    "thermodynamic_temperatures": row.n_thermodynamic_temperature,
                    "thermodynamic_enthalpies": row.n_thermodynamic_enthalpy,
                    "score": row.score,
                    "status": status_for_score(row.score, target_total),
                    "lab_names_150K": "; ".join(row.labs_150k),
                    "lab_names_RT": "; ".join(row.labs_room),
                    "thermodynamic_lab_names": "; ".join(row.thermodynamic_labs),
                    "thermodynamic_record_ids": "; ".join(row.thermodynamic_records),
                }
            )


def make_plot(
    summary: Sequence[SummaryRow],
    output_base: Path,
    *,
    target_total: int,
    minimum_visible_bar: float,
    title: str,
) -> None:
    if not summary:
        raise ValueError("No crystalline forms were detected; no plot was produced.")
    if minimum_visible_bar <= 0:
        raise ValueError("Minimum visible bar length must be greater than zero.")

    # Each crystalline form has its own thin coloured progress bar.  All forms
    # belonging to the same compound are placed inside ONE continuous grey
    # background block, so the grey area identifies the compound group rather
    # than duplicating every coloured bar.
    # Compact vertical layout for embedding the figure in a website.  The
    # values below preserve clear separation between compounds while reducing
    # the overall figure height by roughly one third compared with the earlier
    # layout.
    progress_bar_height = 0.25
    row_step = 0.4
    group_gap = 0.5
    group_vertical_padding = 0.24
    heading_gap = 0.01

    grouped_rows: list[tuple[str, list[SummaryRow]]] = []
    for row in summary:
        if not grouped_rows or grouped_rows[-1][0] != row.compound:
            grouped_rows.append((row.compound, [row]))
        else:
            grouped_rows[-1][1].append(row)

    y_positions: list[float] = []
    row_positions: list[tuple[SummaryRow, float]] = []
    group_geometry: list[tuple[str, float, float, float]] = []
    y = 0.0

    for compound, rows in grouped_rows:
        positions = [y + index * row_step for index in range(len(rows))]
        first_y = positions[0]
        last_y = positions[-1]
        group_top = first_y - group_vertical_padding
        group_bottom = last_y + group_vertical_padding
        group_center = (group_top + group_bottom) / 2.0
        group_height = group_bottom - group_top
        heading_y = group_top - heading_gap
        group_geometry.append((compound, group_center, group_height, heading_y))

        for row, ypos in zip(rows, positions):
            y_positions.append(ypos)
            row_positions.append((row, ypos))

        y = last_y + row_step + group_gap

    max_displayed_score = max(
        max(float(row.score), minimum_visible_bar) for row in summary
    )

    # The common grey width is a visual reference based on the largest current
    # score. It is not a ceiling, and coloured bars may extend beyond the
    # vertical target line. The x axis is linear up to the target and changes
    # to a base-2 logarithmic scale afterwards.
    background_width = max(float(target_total), max_displayed_score)

    # Extend the axis to the next target × 2^n tick, leaving a small margin.
    # This keeps the logarithmic part compact even if one compound has much
    # more data than all the others.
    post_target_ratio = max(1.0, background_width / float(target_total))
    highest_power = max(1, math.ceil(math.log(post_target_ratio, 2)))
    highest_post_tick = float(target_total) * (2**highest_power)
    x_max = highest_post_tick * 1.12

    figure_height = max(
        4.6,
        0.34 * len(summary) + 0.42 * len(grouped_rows) + 1.55,
    )
    fig, ax = plt.subplots(figsize=(12.5, figure_height))

    status_colours = {
        "red": "#D73027",
        "yellow": "#F2C94C",
        "green": "#2E9D55",
    }

    # One continuous grey background block per compound.
    for _compound, group_center, group_height, _heading_y in group_geometry:
        ax.barh(
            group_center,
            background_width,
            height=group_height,
            color="#E3E3E3",
            edgecolor="none",
            zorder=1,
        )

    # One thin coloured bar per crystalline form.
    for row, ypos in row_positions:
        status = status_for_score(row.score, target_total)
        displayed_score = max(float(row.score), minimum_visible_bar)
        ax.barh(
            ypos,
            displayed_score,
            height=progress_bar_height,
            color=status_colours[status],
            edgecolor="none",
            zorder=3,
        )

    # Crystalline-form names are placed before the bars, never inside them.
    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [row.form for row, _ypos in row_positions],
        fontsize=13,
        fontweight="bold",
    )
    ax.tick_params(axis="y", length=0, pad=10)

    # Compound names are larger and positioned above their single grey block.
    for compound, _group_center, _group_height, heading_y in group_geometry:
        ax.text(
            0.0,
            heading_y,
            compound.replace("_", " "),
            transform=ax.get_yaxis_transform(),
            va="bottom",
            ha="left",
            fontsize=13,
            fontweight="bold",
            clip_on=False,
            zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4},
        )

    # Clear, unlabelled target line. Bars remain free to extend beyond it.
    ax.axvline(
        target_total,
        color="black",
        linewidth=1.7,
        linestyle="--",
        zorder=4,
    )

    top_limit = min(heading_y for *_rest, heading_y in group_geometry) - 0.28
    last_group_bottom = max(
        center + height / 2.0 for _compound, center, height, _heading in group_geometry
    )
    # Keep the important 0–target region linear and compress only values above
    # the dashed target line. A base of 2 makes every post-target major tick a
    # doubling of the previous one.
    ax.set_xscale(
        "symlog",
        linthresh=float(target_total),
        linscale=1.0,
        base=2,
    )
    ax.set_xlim(0, x_max)
    ax.set_ylim(last_group_bottom + 0.42, top_limit)

    # Use explicit ticks so the change from the linear region to doubling is
    # immediately readable. The four pre-target ticks divide the linear region
    # into quarters; post-target ticks are target × 2, × 4, × 8, ...
    if float(target_total).is_integer():
        integer_target = int(target_total)
        linear_step = max(1, round(integer_target / 4))
        linear_ticks = [float(value) for value in range(0, integer_target, linear_step)]
        if not linear_ticks or not math.isclose(linear_ticks[-1], float(target_total)):
            linear_ticks.append(float(target_total))
    else:
        linear_ticks = [
            0.0,
            float(target_total) * 0.25,
            float(target_total) * 0.50,
            float(target_total) * 0.75,
            float(target_total),
        ]
    logarithmic_ticks = [
        float(target_total) * (2**power)
        for power in range(1, highest_power + 1)
    ]
    major_ticks = linear_ticks + logarithmic_ticks

    def tick_label(value: float) -> str:
        if math.isclose(value, round(value), abs_tol=1e-9):
            return str(int(round(value)))
        return f"{value:g}"

    ax.xaxis.set_major_locator(FixedLocator(major_ticks))
    ax.xaxis.set_major_formatter(
        FixedFormatter([tick_label(value) for value in major_ticks])
    )

    ax.set_xlabel(
        "Absolute data-coverage score " "\n"
        "Linear to the dashed line; beyond it, each major tick doubles (×2).",
        fontsize=14,
        labelpad=8,
    )
    #ax.set_title(title, pad=12, fontsize=12.5, fontweight="bold")

    # Standard double slash at the transition point makes the axis break clear.
    ax.text(
        float(target_total),
        0.006,
        "//",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        clip_on=False,
        zorder=8,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.05},
    )

    ax.grid(axis="x", linestyle=":", linewidth=0.7, alpha=0.55, zorder=0)
    ax.spines[["top", "right", "left"]].set_visible(False)

    # Reserve enough room for form labels and long compound names.
    fig.subplots_adjust(left=0.22, right=0.98, top=0.94, bottom=0.14)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

def create_alias_template(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["alias", "lab"])
        writer.writerow(["# Example: D_Braun", "# Example: Innsbruck"])
        writer.writerow(["# Example: W_Wood", "# Example: Innsbruck"])


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------
def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate BEST-CSP per-form progress CSV files and progress bars."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Repository root containing one directory per compound (default: current directory).",
    )
    parser.add_argument(
        "--output-dir",
        default="progress",
        help="Output directory, relative to ROOT unless absolute (default: progress).",
    )
    parser.add_argument(
        "--aliases",
        default="lab_aliases.csv",
        help="Alias CSV with columns alias,lab; relative to ROOT unless absolute.",
    )
    parser.add_argument("--target-150", type=float, default=DEFAULT_TARGET_150_K)
    parser.add_argument("--tolerance-150", type=float, default=DEFAULT_TOLERANCE_150_K)
    parser.add_argument("--room-temperature", type=float, default=DEFAULT_ROOM_TEMPERATURE_K)
    parser.add_argument("--room-tolerance", type=float, default=DEFAULT_ROOM_TOLERANCE_K)
    parser.add_argument(
        "--required-150-labs",
        type=int,
        default=DEFAULT_REQUIRED_150_LABS,
        help="Target number of laboratories near 150 K (default: 3).",
    )
    parser.add_argument(
        "--required-room-labs",
        type=int,
        default=DEFAULT_REQUIRED_ROOM_LABS,
        help="Target number of laboratories near room temperature (default: 3).",
    )
    parser.add_argument(
        "--required-thermodynamic-data",
        type=int,
        default=DEFAULT_REQUIRED_THERMODYNAMIC_DATA,
        help="Target number of thermodynamic data records (default: 5).",
    )
    parser.add_argument(
        "--minimum-visible-bar",
        type=float,
        default=DEFAULT_MINIMUM_VISIBLE_BAR,
        help="Smallest displayed bar length; the CSV score remains unchanged (default: 1).",
    )
    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Also scan directories named Supplementary_data.",
    )
    parser.add_argument(
        "--exclude-melting",
        action="store_true",
        help="Exclude melting points and fusion/melting enthalpies.",
    )
    parser.add_argument(
        "--exclude-sublimation-enthalpy",
        action="store_true",
        help="Exclude sublimation enthalpy from thermodynamic-data counts.",
    )
    parser.add_argument(
        "--include-glass-transition",
        action="store_true",
        help="Count glass-transition measurements when assigned to a form.",
    )
    parser.add_argument(
        "--title",
        default="BEST-CSP experimental-data progress",
        help="Plot title.",
    )
    return parser.parse_args(argv)


def resolve_against_root(root: Path, path_text: str) -> Path:
    path = Path(path_text).expanduser()
    return path if path.is_absolute() else root / path


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    root = Path(args.root).expanduser().resolve()
    if not root.is_dir():
        print(f"ERROR: repository root does not exist or is not a directory: {root}", file=sys.stderr)
        return 2

    for name, value in (
        ("required 150 K laboratories", args.required_150_labs),
        ("required room-temperature laboratories", args.required_room_labs),
        ("required thermodynamic data", args.required_thermodynamic_data),
    ):
        if value < 0:
            print(f"ERROR: {name} cannot be negative.", file=sys.stderr)
            return 2

    output_dir = resolve_against_root(root, args.output_dir)
    aliases_path = resolve_against_root(root, args.aliases)
    create_alias_template(aliases_path)

    try:
        aliases = load_aliases(aliases_path)
        evidence, discovered_forms = scan_repository(
            root,
            aliases=aliases,
            target_150_k=args.target_150,
            tolerance_150_k=args.tolerance_150,
            room_temperature_k=args.room_temperature,
            room_tolerance_k=args.room_tolerance,
            include_supplementary=args.include_supplementary,
            include_melting=not args.exclude_melting,
            include_glass_transition=args.include_glass_transition,
            include_sublimation_enthalpy=not args.exclude_sublimation_enthalpy,
        )
        summary = build_summary(evidence, discovered_forms)
        target_total = (
            args.required_150_labs
            + args.required_room_labs
            + args.required_thermodynamic_data
        )

        write_evidence_csv(output_dir / "progress_evidence.csv", evidence)
        write_summary_csv(output_dir / "progress_summary.csv", summary, target_total)
        make_plot(
            summary,
            output_dir / "progress_bars",
            target_total=target_total,
            minimum_visible_bar=args.minimum_visible_bar,
            title=args.title,
        )
    except (OSError, ValueError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Scanned repository: {root}")
    print(f"Detected forms: {len(summary)}")
    print(f"Accepted evidence rows: {len(evidence)}")
    print(f"Target line: {target_total}")
    print(f"Summary: {output_dir / 'progress_summary.csv'}")
    print(f"Evidence: {output_dir / 'progress_evidence.csv'}")
    print(f"Plot: {output_dir / 'progress_bars.svg'}")
    print(f"Plot: {output_dir / 'progress_bars.png'}")
    print(f"Laboratory aliases: {aliases_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
