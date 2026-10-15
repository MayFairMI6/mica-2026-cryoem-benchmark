#!/usr/bin/env python
"""Run with `phenix.python phenix_adapter.py reference_CA.cif model_CA.cif out.json`.

Uses the official cctbx chain_comparison implementation and explicit denominators.
Executed here with cctbx-base 2025.11 through evaluate.run_chain_comparison.
It can also run with phenix.python, with that release recorded and validated.
This is a declared convention, not a recovered original MICA evaluation command.
"""
from __future__ import print_function
import hashlib
import json
import sys
import cctbx
from mmtbx.validation import chain_comparison


def main():
    reference, prediction, output = sys.argv[1:4]
    args = [reference, prediction,
            "unique_part_of_target_only=False",
            "test_unique_part_of_target_only=False",
            "unique_query_only=False",
            "use_crystal_symmetry=False",
            "target_length_from_matching_chains=False",
            "chain_type=PROTEIN", "max_dist=3.0", "remove_alt_conf=True"]
    with open(output + ".raw.log", "w") as log:
        params = chain_comparison.get_params(args, out=log)
        rv = chain_comparison.run(params=params, out=log)
        rv.show_summary(out=log)
    rmsd, nclose = rv.get_values("close")
    _, nfar = rv.get_values("far_away")
    nclose, nfar = int(nclose or 0), int(nfar or 0)
    nref = int(rv.used_target)
    npred = int(rv.used_chain)
    precision = nclose / float(npred) if npred else 0.0
    result = {
        "status": "executed_official_cctbx_chain_comparison",
        "implementation_version": cctbx.__version__,
        "implementation_source_sha256": hashlib.sha256(open(chain_comparison.__file__, "rb").read()).hexdigest(),
        "implementation_note": "Official cctbx mmtbx.validation.chain_comparison module; full Phenix suite is not required for this calculation",
        "original_mica_command_recovered": False,
        "arguments": args, "reference_ca_used": nref, "predicted_ca_used": npred,
        "matched_ca": nclose, "far_ca": nfar,
        "ca_match_percent_paper_prose_denominator": 100.0*precision,
        "ca_quality_fraction_paper_formula": precision*npred/nref if nref else None,
        "sequence_match_percent": rv.get_match_percent("close"),
        "phenix_percent_found_reference_denominator": rv.get_close_to_target_percent("close"),
        "matched_ca_rmsd_angstrom": rmsd,
        "count_reconciliation": nclose+nfar==npred,
        "note": "Phenix printed %Found uses reference length; CAlphaMatch above uses prediction length as specified in MICA prose."
    }
    with open(output, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
