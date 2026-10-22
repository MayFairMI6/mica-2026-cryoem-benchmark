"""Check the frozen release without running model inference.

AI-tool assistance with project scripts is disclosed in README.md.
This check validates archived records; it does not re-query release dates
or establish independence from a model's complete training corpus.
"""
import csv
import json
import statistics
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    with (root / "dataset/case_rounds.csv").open() as stream:
        rounds = list(csv.DictReader(stream))
    cases = {row["case_id"] for row in rounds}
    if len(rounds) != 24 or len(cases) != 24:
        raise SystemExit("Expected 24 unique benchmark cases")
    with (root / "dataset/primary24_scope_screened_candidate.csv").open() as stream:
        selected = [row for row in csv.DictReader(stream) if row["case_id"] in cases]
    if len(selected) != 24 or any(
        not row[key].startswith("2026-")
        for row in selected
        for key in ("pdb_initial_release_date", "emdb_initial_map_release_date")
    ):
        raise SystemExit("Missing cases or unexpected recorded release dates")
    expected = {"MICA": .857, "EModelX(+AF)": .838,
                "ModelAngelo": .683, "CryoAtom": .721, "EMProt": .575}
    means = {}
    for method, published in expected.items():
        results = json.loads((root / "scoring/results" / (method + ".json")).read_text())
        scores = [results[case]["metrics"]["tm_score"] for case in sorted(cases)]
        if not all(isinstance(value, (int, float)) and 0 <= value <= 1 for value in scores):
            raise SystemExit("Invalid TM-score: " + method)
        means[method] = statistics.mean(scores)
        if round(means[method], 3) != published:
            raise SystemExit("Reported mean differs: " + method)
    print(json.dumps({"cases": len(cases), "primary_scores": 120,
                      "recorded_release_year": 2026, "mean_tm_score": means}, indent=2))


if __name__ == "__main__":
    main()
