#!/usr/bin/env python3
"""Protein model evaluation with explicit paper/tool and supplementary fields.

Experimental reference structures are evaluation-only inputs. This module does
not create model predictions, dock models, fit to density, or infer a structure.

US-align fields reproduce the named metric definitions with a pinned, explicit
complex-alignment convention; MICA's original US-align command is unavailable.
The spatial fields are added safeguards, not a reimplementation of Phenix.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
if (HERE / "vendor").exists():
    sys.path.insert(0, str(HERE / "vendor"))
import gemmi
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def resolve_recorded_path(recorded: str) -> Path:
    """Rebase preserved absolute provenance after unpacking the project archive."""
    path = Path(recorded)
    parts = path.parts
    if "mica_benchmark" in parts:
        local = HERE.parent.joinpath(*parts[parts.index("mica_benchmark") + 1:])
        if local.exists():
            return local
    return path


@dataclass(frozen=True)
class CA:
    chain: str
    auth_chain: str
    label_chain: str
    seqid: str
    label_seq: int | None
    resname: str
    aa: str
    xyz: tuple[float, float, float]
    occupancy: float
    altloc: str


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_ca(path: str | Path, polymer_only: bool = False) -> tuple[list[CA], dict]:
    """Read first coordinate model; retain one occupied CA per amino-acid residue.

    Gemmi handles mmCIF quoting, long chain IDs, insertion codes, modified amino
    acids and PDB files. Chain keys preserve both author and label chain identity.
    Alternate locations use highest occupancy, then blank, A, alphabetical order.
    Zero-occupancy/non-finite atoms and non-amino-acid residues are excluded.
    No missing residues or biological copies are invented from metadata.
    polymer_only (used for the deposited reference, 2026-10-13): residues of non-polymer entities are excluded, so a
    free amino-acid ligand (e.g. the four GLU ligands of 9P9D) is not counted as a protein residue.
    """
    path = Path(path)
    structure = gemmi.read_structure(str(path))
    if not len(structure):
        raise ValueError(f"No coordinate models in {path}")
    sites: list[CA] = []
    excluded = Counter()
    alt_choices = 0
    seen = set()
    if polymer_only:
        structure.setup_entities()
    for chain in structure[0]:
        for residue in chain:
            if polymer_only and residue.entity_type != gemmi.EntityType.Polymer:
                excluded["non_polymer_entity_residues"] += 1
                continue
            info = gemmi.find_tabulated_residue(residue.name)
            if not info.is_amino_acid():
                excluded["non_amino_acid_residues"] += 1
                continue
            candidates = [a for a in residue if a.name.strip() == "CA"]
            candidates = [a for a in candidates if a.occ > 0 and all(
                np.isfinite([a.pos.x, a.pos.y, a.pos.z]))]
            if not candidates:
                excluded["amino_acid_residues_without_occupied_finite_CA"] += 1
                continue
            if len(candidates) > 1:
                alt_choices += 1
            def alt_rank(atom):
                alt = atom.altloc.replace("\x00", "").strip()
                return (-float(atom.occ), 0 if alt == "" else 1 if alt == "A" else 2, alt)
            atom = sorted(candidates, key=alt_rank)[0]
            label = residue.subchain or ""
            chain_key = f"auth={chain.name}|label={label}"
            key = (chain_key, str(residue.seqid))
            if key in seen:
                raise ValueError(f"Duplicate residue identity in {path}: {key}; resolve microheterogeneity explicitly")
            seen.add(key)
            aa = info.one_letter_code.upper()
            if aa not in CANONICAL_AA:
                aa = "X"
            sites.append(CA(chain_key, chain.name, label, str(residue.seqid),
                            residue.label_seq, residue.name, aa,
                            (atom.pos.x, atom.pos.y, atom.pos.z),
                            float(atom.occ), atom.altloc.replace("\x00", "").strip()))
    if not sites:
        raise ValueError(f"No occupied protein CA atoms in {path}")
    xyz = np.array([a.xyz for a in sites])
    return sites, {
        "path": str(path.resolve()), "sha256": sha256(path),
        "coordinate_model_index": 0, "coordinate_models_in_file": len(structure),
        "protein_ca_count": len(sites),
        "chain_ca_counts": dict(Counter(a.chain for a in sites)),
        "alternate_location_choices": alt_choices,
        "unknown_amino_acid_count": sum(a.aa == "X" for a in sites),
        "excluded": dict(excluded),
        "coordinate_bbox_angstrom": [xyz.min(axis=0).tolist(), xyz.max(axis=0).tolist()],
    }


def write_ca_cif(sites: list[CA], path: str | Path) -> None:
    """Canonical CA-only mmCIF for consistent inputs to US-align and controls.

    Sequential chain IDs and label residue numbers are local evaluation IDs;
    original identifiers remain in the JSON metadata of the input files.
    """
    document = gemmi.cif.Document()
    block = document.add_new_block("evaluation_CA")
    tags = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
            "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
            "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
            "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
            "auth_atom_id", "pdbx_PDB_model_num"]
    loop = block.init_loop("_atom_site.", tags)
    chain_map = {name: f"C{i+1}" for i, name in enumerate(dict.fromkeys(a.chain for a in sites))}
    chain_counts = Counter()
    for i, site in enumerate(sites, 1):
        chain_counts[site.chain] += 1
        resnum = str(chain_counts[site.chain])
        ch = chain_map[site.chain]
        # Standardise recognised modified amino acids to parent residue names.
        names = {"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY",
                 "H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN",
                 "P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL",
                 "W":"TRP","Y":"TYR","X":"UNK"}
        name = names[site.aa]
        loop.add_row(["ATOM", str(i), "C", "CA", ".", name, ch,
                      str(list(chain_map).index(site.chain)+1), resnum, "?",
                      *[f"{x:.6f}" for x in site.xyz], "1.0", "0.0", resnum,
                      name, ch, "CA", "1"])
    document.write_file(str(path))


def spatial_one_to_one(predicted: list[CA], reference: list[CA], threshold: float) -> dict:
    """Maximum-cardinality, then minimum-distance injective spatial matching.

    This is deliberately sequence-independent and ignores chain labels, testing
    coordinate recovery only. Chain order/topology is assessed separately by
    US-align. A sparse assignment with one dummy column per prediction avoids
    a dense all-by-all distance matrix and penalises duplicate predictions.
    """
    px = np.array([a.xyz for a in predicted], dtype=float)
    rx = np.array([a.xyz for a in reference], dtype=float)
    npred, nref = len(px), len(rx)
    near = cKDTree(px).sparse_distance_matrix(cKDTree(rx), threshold,
                                             output_type="coo_matrix")
    # Every real edge costs between 1 and 2. Losing one match costs more than
    # all possible improvements in real-edge costs combined.
    dummy_penalty = 3.0 * (min(npred, nref) + 1)
    rows = np.concatenate([near.row, np.arange(npred)])
    cols = np.concatenate([near.col, nref + np.arange(npred)])
    costs = np.concatenate([1.0 + near.data / threshold,
                            np.full(npred, dummy_penalty)])
    graph = coo_matrix((costs, (rows, cols)), shape=(npred, nref + npred)).tocsr()
    ri, ci = min_weight_full_bipartite_matching(graph)
    keep = ci < nref
    ri, ci = ri[keep], ci[keep]
    d = np.linalg.norm(px[ri] - rx[ci], axis=1)
    matched = len(ri)
    identical = sum(predicted[i].aa == reference[j].aa and predicted[i].aa != "X"
                    for i, j in zip(ri, ci))
    p = matched / npred
    r = matched / nref
    return {
        "threshold_angstrom_inclusive": threshold,
        "matching": "maximum cardinality then minimum total distance; one-to-one; sequence-independent",
        "matched_ca": matched, "predicted_ca": npred, "reference_ca": nref,
        "precision": p, "recall": r, "f1": 2*p*r/(p+r) if p+r else 0.0,
        "identical_amino_acid_pairs": int(identical),
        "sequence_identity_among_spatial_matches": identical / matched if matched else None,
        "sequence_correct_reference_coverage": identical / nref,
        "matched_pair_rmsd_angstrom": float(np.sqrt(np.mean(d*d))) if matched else None,
        "matched_distance_max_angstrom": float(np.max(d)) if matched else None,
        "excess_prediction_fraction": max(0.0, (npred-nref)/nref),
        "candidate_edges": int(near.nnz),
    }


def run_usalign(predicted: list[CA], reference: list[CA], executable: Path,
                output: Path, timeout: int = 180) -> tuple[dict, list[CA] | None]:
    # Time limit only (no metric is affected): large extension complexes (24 reference chains) exceed the 180 s default.
    timeout = int(os.environ.get("MICA_EVAL_TIMEOUT_S", timeout))
    output.mkdir(parents=True, exist_ok=True)
    pred_path = output / "prediction_CA.cif"
    ref_path = output / "reference_CA.cif"
    write_ca_cif(predicted, pred_path)
    write_ca_cif(reference, ref_path)
    matrix_path = output / "usalign_matrix.txt"
    cmd = [str(executable.resolve()), str(pred_path.resolve()), str(ref_path.resolve()),
           "-mol", "prot", "-mm", "1", "-ter", "1", "-m", str(matrix_path.resolve())]
    (output / "usalign_command.json").write_text(json.dumps(cmd, indent=2)+"\n")
    started = time.monotonic()
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        (output / "usalign_stdout.txt").write_text(str(error.stdout or ""))
        return {"status": "timeout", "timeout_seconds": timeout, "command": cmd}, None
    (output / "usalign_stdout.txt").write_text(process.stdout)
    (output / "usalign_stderr.txt").write_text(process.stderr)
    result = {"status": "ok" if process.returncode == 0 else "failed",
              "returncode": process.returncode, "wall_seconds": time.monotonic()-started,
              "command": cmd, "executable_sha256": sha256(executable),
              "normalization": "Structure_2: full retained experimental reference CA length",
              "alignment": "US-align -mol prot -mm 1 -ter 1; one global complex transform"}
    if process.returncode:
        return result, None
    aligned = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+),\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)", process.stdout)
    tm = re.search(r"TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_2:\s*L=(\d+)", process.stdout)
    pred_length = re.search(r"Length of Structure_1:\s*(\d+)", process.stdout)
    if not aligned or not tm or not pred_length:
        result["status"] = "parse_failed"
        return result, None
    result.update({"tm_score_reference_normalized": float(tm[1]),
                   "aligned_ca_length": int(aligned[1]),
                   "sequence_identity": float(aligned[3]),
                   "aligned_rmsd_angstrom": float(aligned[2]),
                   "reference_ca_count_used": int(tm[2]),
                   "predicted_ca_count_used": int(pred_length[1])})
    if int(tm[2]) != len(reference) or int(pred_length[1]) != len(predicted):
        result["status"] = "length_mismatch"
        result["warning"] = "US-align discarded atoms; resolve before using score"
    rows = []
    for line in matrix_path.read_text().splitlines():
        fields = line.split()
        if len(fields) == 5 and fields[0] in ["0", "1", "2"]:
            rows.append([float(x) for x in fields[1:]])
    if len(rows) != 3:
        result["transform_status"] = "parse_failed"
        return result, None
    matrix = np.array(rows)
    rotation, translation = matrix[:, 1:], matrix[:, 0]
    result["rotation_determinant"] = float(np.linalg.det(rotation))
    result["transform_status"] = "ok"
    xyz = np.array([a.xyz for a in predicted]) @ rotation.T + translation
    moved = [replace(a, xyz=tuple(x)) for a, x in zip(predicted, xyz)]
    return result, moved


def run_chain_comparison(predicted: list[CA], reference: list[CA], output: Path) -> dict:
    """Run the official cctbx module installed in an isolated optional directory.

    The package supplies Phenix's chain_comparison implementation without
    supplying or claiming to install the full Phenix suite.
    """
    vendor = HERE / "cctbx_vendor"
    if not (vendor / "mmtbx").exists():
        return {"status": "chain_comparison_unavailable", "ca_match_percent": None,
                "ca_quality_fraction": None, "sequence_match_percent": None,
                "original_mica_command_recovered": False,
                "reason": "Install cctbx-base==2025.11 in cctbx_vendor or run phenix_adapter.py in Phenix"}
    output.mkdir(parents=True, exist_ok=True)
    ref_path, pred_path = output / "reference_CA.cif", output / "prediction_CA.cif"
    write_ca_cif(reference, ref_path)
    write_ca_cif(predicted, pred_path)
    result_path = output / "chain_comparison.json"
    cmd = [sys.executable, str(HERE / "phenix_adapter.py"), str(ref_path.resolve()),
           str(pred_path.resolve()), str(result_path.resolve())]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(vendor)
    env["LIBTBX_BUILD"] = str(vendor / "libtbx" / "core" / "share" / "cctbx")
    (output / "command.json").write_text(json.dumps(cmd, indent=2)+"\n")
    try:
        cc_timeout = int(os.environ.get("MICA_EVAL_TIMEOUT_S", 180))
        process = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=cc_timeout)
    except subprocess.TimeoutExpired:
        return {"status": "chain_comparison_timeout", "command": cmd, "timeout_seconds": cc_timeout}
    (output / "stdout.txt").write_text(process.stdout)
    (output / "stderr.txt").write_text(process.stderr)
    if process.returncode:
        return {"status": "chain_comparison_failed", "command": cmd, "returncode": process.returncode}
    return json.loads(result_path.read_text())


def evaluate(pred_path: Path, ref_path: Path, output: Path,
             thresholds=(1.0, 2.0, 3.0, 5.0), executable: Path | None = None) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    predicted, pmeta = read_ca(pred_path)
    reference, rmeta = read_ca(ref_path, polymer_only=True)
    result = {
        "schema_version": "1.0", "experiment_type": "model_evaluation",
        "prediction": pmeta, "reference": rmeta,
        "definitions": {
            "map_frame": "Original Cartesian coordinates; no centering, scaling, symmetry expansion, or rigid fit",
            "superposed": "Single global US-align rigid transform; report separately from map-frame metrics",
            "protein_subset": "First model, observed occupied amino-acid CA only; fixed original stoichiometry",
            "sequence_identity_unknowns": "Supplementary metrics count X as non-identical; US-align identity follows tool",
        },
        "supplementary_map_frame": [spatial_one_to_one(predicted, reference, t) for t in thresholds],
        "phenix_paper_metrics": run_chain_comparison(predicted, reference, output / "chain_comparison"),
    }
    if executable and executable.exists():
        us, aligned = run_usalign(predicted, reference, executable, output / "usalign")
        result["paper_named_usalign_metrics"] = us
        if aligned is not None:
            result["supplementary_after_global_superposition"] = [
                spatial_one_to_one(aligned, reference, t) for t in thresholds]
    else:
        result["paper_named_usalign_metrics"] = {"status": "executable_unavailable"}
    (output / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--usalign", type=Path, default=HERE / "USalign" / "USalign")
    parser.add_argument("--thresholds", type=float, nargs="+", default=[1, 2, 3, 5])
    args = parser.parse_args()
    if any(t <= 0 for t in args.thresholds):
        parser.error("Distance thresholds must be positive")
    result = evaluate(args.prediction, args.reference, args.output, args.thresholds, args.usalign)
    print(json.dumps({"metrics": str((args.output / 'metrics.json').resolve()),
                      "usalign_status": result['paper_named_usalign_metrics']['status'],
                      "predicted_ca": result['prediction']['protein_ca_count'],
                      "reference_ca": result['reference']['protein_ca_count']}))


if __name__ == "__main__":
    main()
