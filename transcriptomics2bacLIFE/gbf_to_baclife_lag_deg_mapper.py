#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Map GBFF gene identifiers to bacLIFE combined_proteins identifiers via "
            "DIAMOND, then join the resulting cluster mapping with bacLIFE LAGs and DEGs."
        )
    )
    parser.add_argument("--baclife-LAGs", dest="baclife_lags", required=True,
                        help="Path to plant_pathogen_LAGs.xlsx")
    parser.add_argument("--megamatrix", required=True,
                        help="Path to MEGAMATRIX_renamed.txt")
    parser.add_argument("--combined-proteins", required=True,
                        help="Path to combined_proteins fasta")
    parser.add_argument("--DEGs", dest="degs", required=True,
                        help="Path to DEGs.csv")
    parser.add_argument("--gbf", required=True,
                        help="Input GenBank/GBFF file")
    parser.add_argument("--gbf-tag", required=True,
                        help="CDS qualifier to use as FASTA header, e.g. locus_tag, old_locus_tag, gene, protein_id")
    parser.add_argument("--output", required=True,
                        help="Output directory")
    parser.add_argument("--identity-threshold", type=float, default=90.0,
                        help="Minimum percent identity to keep a DIAMOND hit [default: 90.0]")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for DIAMOND [default: 1]")
    parser.add_argument("--keep-temp", action="store_true",
                        help="Keep intermediate files")
    return parser.parse_args()


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def check_executable(name: str):
    try:
        subprocess.run([name, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    except FileNotFoundError:
        sys.exit(f"ERROR: Required executable not found in PATH: {name}")


def parse_megamatrix(megamatrix_path: Path) -> Dict[str, List[str]]:
    gene_to_clusters: Dict[str, List[str]] = {}
    with megamatrix_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=' ', quotechar='"', skipinitialspace=True)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Empty megamatrix file: {megamatrix_path}")

        try:
            cluster_idx = header.index("clusters")
            desc_idx = header.index("descriptions")
        except ValueError as exc:
            raise ValueError("The megamatrix file must contain 'clusters' and 'descriptions' columns.") from exc

        for line_number, fields in enumerate(reader, start=2):
            if not fields:
                continue
            if len(fields) <= max(cluster_idx, desc_idx):
                sys.stderr.write(
                    f"WARNING: Skipping malformed megamatrix line {line_number}: "
                    f"expected >= {max(cluster_idx, desc_idx) + 1} fields, got {len(fields)}\n"
                )
                continue

            cluster_id = fields[cluster_idx].strip()
            descriptions = fields[desc_idx].strip()
            if not descriptions or descriptions.upper() == "NA":
                continue

            genes = [gene.strip() for gene in descriptions.split(",") if gene.strip()]
            for gene in genes:
                gene_to_clusters.setdefault(gene, []).append(cluster_id)

    return gene_to_clusters


def extract_proteins_from_gbf(gbf_path: Path, header_tag: str, output_faa: Path) -> int:
    written = 0
    records: List[SeqRecord] = []

    for record in SeqIO.parse(str(gbf_path), "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            qualifiers = feature.qualifiers
            if "translation" not in qualifiers:
                continue
            if header_tag not in qualifiers:
                continue

            header_value = qualifiers[header_tag][0].strip()
            protein_seq = qualifiers["translation"][0].strip()
            if not header_value or not protein_seq:
                continue

            records.append(SeqRecord(Seq(protein_seq), id=header_value, name=header_value, description=""))
            written += 1

    if written == 0:
        raise ValueError(f"No CDS proteins were extracted from {gbf_path} using qualifier '{header_tag}'.")
    SeqIO.write(records, str(output_faa), "fasta")
    return written


def run_command(cmd: List[str], description: str):
    sys.stderr.write(f"[RUN] {description}\n")
    sys.stderr.write("[CMD] " + " ".join(cmd) + "\n")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def make_diamond_db(reference_fasta: Path, db_prefix: Path):
    run_command(
        ["diamond", "makedb", "--in", str(reference_fasta), "--db", str(db_prefix)],
        "Building DIAMOND database from combined_proteins"
    )


def run_diamond(query_faa: Path, db_prefix: Path, out_tsv: Path, threads: int):
    run_command(
        [
            "diamond", "blastp",
            "--query", str(query_faa),
            "--db", str(db_prefix),
            "--out", str(out_tsv),
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
            "--max-target-seqs", "25",
            "--threads", str(threads),
            "--quiet",
        ],
        "Running DIAMOND blastp"
    )


def read_best_hits(diamond_tsv: Path, identity_threshold: float) -> Dict[str, Tuple[str, float, int, float, float]]:
    best_hits: Dict[str, Tuple[str, float, int, float, float]] = {}
    with diamond_tsv.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                sys.stderr.write(f"WARNING: Skipping malformed DIAMOND line {line_number}\n")
                continue

            qseqid, sseqid, pident, length, evalue, bitscore = parts
            pident_f = float(pident)
            if pident_f < identity_threshold:
                continue

            candidate = (sseqid, pident_f, int(length), float(evalue), float(bitscore))
            current = best_hits.get(qseqid)

            if current is None:
                best_hits[qseqid] = candidate
                continue

            replace = False
            if candidate[4] > current[4]:
                replace = True
            elif candidate[4] == current[4]:
                if candidate[3] < current[3]:
                    replace = True
                elif candidate[3] == current[3]:
                    if candidate[1] > current[1]:
                        replace = True
                    elif candidate[1] == current[1]:
                        if candidate[2] > current[2]:
                            replace = True
                        elif candidate[2] == current[2] and candidate[0] < current[0]:
                            replace = True

            if replace:
                best_hits[qseqid] = candidate
    return best_hits


def write_mapping_table(
    output_tsv: Path,
    best_hits: Dict[str, Tuple[str, float, int, float, float]],
    gene_to_clusters: Dict[str, List[str]],
):
    rows: List[Dict[str, str]] = []
    missing_cluster = 0

    with output_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["clusters", "baclife_gene", "gene_name"])

        for gene_name in sorted(best_hits):
            baclife_gene = best_hits[gene_name][0]
            clusters = gene_to_clusters.get(baclife_gene)

            if not clusters:
                missing_cluster += 1
                continue

            for cluster_id in sorted(set(clusters)):
                writer.writerow([cluster_id, baclife_gene, gene_name])
                rows.append(
                    {"clusters": cluster_id, "baclife_gene": baclife_gene, "gene_name": gene_name}
                )

    return rows, missing_cluster


def col_to_index(cell_ref: str) -> int:
    letters = ''.join(ch for ch in cell_ref if ch.isalpha())
    index = 0
    for ch in letters:
        index = index * 26 + (ord(ch.upper()) - ord('A') + 1)
    return index - 1


def xlsx_cell_value(cell, ns, shared_strings):
    cell_type = cell.attrib.get("t")
    v = cell.find("a:v", ns)
    inline = cell.find("a:is", ns)

    if cell_type == "s" and v is not None:
        return shared_strings[int(v.text)]
    if cell_type == "inlineStr" and inline is not None:
        return ''.join((t.text or '') for t in inline.iterfind(".//a:t", ns))
    if v is not None:
        return v.text or ""
    return ""


def read_xlsx_first_sheet(xlsx_path: Path) -> List[Dict[str, str]]:
    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(xlsx_path) as zf:
        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", ns):
                shared_strings.append(''.join((t.text or '') for t in si.iterfind(".//a:t", ns)))

        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        sheets = workbook.find("a:sheets", ns)
        if sheets is None or len(list(sheets)) == 0:
            raise ValueError(f"No worksheets found in {xlsx_path}")

        first_sheet_xml = "xl/worksheets/sheet1.xml"
        sheet = ET.fromstring(zf.read(first_sheet_xml))
        sheet_data = sheet.find("a:sheetData", ns)
        if sheet_data is None:
            return []

        sparse_rows: List[List[str]] = []
        for row in sheet_data.findall("a:row", ns):
            values: List[str] = []
            for cell in row.findall("a:c", ns):
                idx = col_to_index(cell.attrib.get("r", "A1"))
                while len(values) <= idx:
                    values.append("")
                values[idx] = xlsx_cell_value(cell, ns, shared_strings).strip()
            sparse_rows.append(values)

    # detect header row: first row containing "id"
    header_idx = None
    for i, row in enumerate(sparse_rows):
        if "id" in row:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Could not find header row with 'id' in {xlsx_path}")

    header = sparse_rows[header_idx]
    records: List[Dict[str, str]] = []
    for row in sparse_rows[header_idx + 1:]:
        if not any(value != "" for value in row):
            continue
        padded = row + [""] * (len(header) - len(row))
        record = {header[j]: padded[j] for j in range(len(header)) if header[j] != ""}
        records.append(record)
    return records


def read_deg_csv(csv_path: Path) -> List[Dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def join_annotations(
    mapping_rows: List[Dict[str, str]],
    lags_rows: List[Dict[str, str]],
    deg_rows: List[Dict[str, str]],
) -> List[Dict[str, str]]:
    lags_by_cluster: Dict[str, List[Dict[str, str]]] = {}
    for row in lags_rows:
        cluster = row.get("id", "").strip()
        if cluster:
            lags_by_cluster.setdefault(cluster, []).append(row)

    deg_by_id: Dict[str, List[Dict[str, str]]] = {}
    for row in deg_rows:
        gene_id = row.get("ID", "").strip()
        if gene_id:
            deg_by_id.setdefault(gene_id, []).append(row)

    joined_rows: List[Dict[str, str]] = []
    excluded_lag_columns = {"id", "gene", "completeness"}

    for mapping in mapping_rows:
        cluster = mapping["clusters"]
        gene_name = mapping["gene_name"]

        matching_lags = lags_by_cluster.get(cluster)
        matching_degs = deg_by_id.get(gene_name)

        # Inner join: skip this mapping row entirely if either side is absent
        if not matching_lags or not matching_degs:
            continue

        for lag in matching_lags:
            for deg in matching_degs:
                merged = {}
                merged.update(mapping)

                for key, value in lag.items():
                    if key in excluded_lag_columns:
                        continue
                    merged[f"LAGs_{key}"] = value

                for key, value in deg.items():
                    merged[f"DEGs_{key}"] = value

                joined_rows.append(merged)

    return joined_rows


def write_table(rows: List[Dict[str, str]], output_path: Path):
    fieldnames: List[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def cleanup_tempdir(tempdir: Path):
    if not tempdir.exists():
        return
    for path in tempdir.iterdir():
        path.unlink()
    tempdir.rmdir()


def main():
    args = parse_args()
    check_executable("diamond")

    outdir = Path(args.output).resolve()
    ensure_dir(outdir)

    megamatrix = Path(args.megamatrix).resolve()
    combined_proteins = Path(args.combined_proteins).resolve()
    gbf = Path(args.gbf).resolve()
    lags_xlsx = Path(args.baclife_lags).resolve()
    degs_csv = Path(args.degs).resolve()

    tempdir = outdir / "tmp_baclife_gene_mapping"
    ensure_dir(tempdir)

    query_faa = tempdir / f"{gbf.stem}.{args.gbf_tag}.proteins.faa"
    diamond_db_prefix = tempdir / "combined_proteins_db"
    diamond_out = tempdir / "gbf_vs_combined_proteins.diamond.tsv"
    mapping_tsv = outdir / "baclife_gene_cluster_mapping.tsv"
    joined_tsv = outdir / "baclife_LAGs_DEGs_merged.tsv"

    sys.stderr.write("Parsing megamatrix...\n")
    gene_to_clusters = parse_megamatrix(megamatrix)
    sys.stderr.write(f"Loaded cluster annotations for {len(gene_to_clusters)} bacLIFE genes.\n")

    sys.stderr.write("Extracting proteins from GBF...\n")
    n_queries = extract_proteins_from_gbf(gbf, args.gbf_tag, query_faa)
    sys.stderr.write(f"Extracted {n_queries} protein sequences.\n")

    make_diamond_db(combined_proteins, diamond_db_prefix)
    run_diamond(query_faa, diamond_db_prefix, diamond_out, args.threads)

    sys.stderr.write("Selecting best DIAMOND hits...\n")
    best_hits = read_best_hits(diamond_out, args.identity_threshold)
    sys.stderr.write(f"Retained {len(best_hits)} query proteins with identity >= {args.identity_threshold:.1f}%.\n")

    sys.stderr.write("Writing mapping table...\n")
    mapping_rows, n_missing_cluster = write_mapping_table(mapping_tsv, best_hits, gene_to_clusters)
    sys.stderr.write(f"Wrote {len(mapping_rows)} mapping rows to {mapping_tsv}\n")
    if n_missing_cluster:
        sys.stderr.write(f"NOTE: {n_missing_cluster} mapped baclife genes had no cluster in megamatrix descriptions.\n")

    sys.stderr.write("Reading bacLIFE LAGs Excel file...\n")
    lags_rows = read_xlsx_first_sheet(lags_xlsx)
    sys.stderr.write(f"Loaded {len(lags_rows)} LAG rows.\n")

    sys.stderr.write("Reading DEGs CSV file...\n")
    deg_rows = read_deg_csv(degs_csv)
    sys.stderr.write(f"Loaded {len(deg_rows)} DEG rows.\n")

    sys.stderr.write("Joining mapping + LAGs + DEGs...\n")
    joined_rows = join_annotations(mapping_rows, lags_rows, deg_rows)
    write_table(joined_rows, joined_tsv)
    sys.stderr.write(f"Wrote {len(joined_rows)} joined rows to {joined_tsv}\n")

    if not args.keep_temp:
        cleanup_tempdir(tempdir)

    sys.stderr.write("Done.\n")


if __name__ == "__main__":
    main()
