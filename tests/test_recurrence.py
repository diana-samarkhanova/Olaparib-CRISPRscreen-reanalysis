"""Synthetic checks for risks relevant to manuscript recurrence reporting."""
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "collect_hits.py"
SPEC = importlib.util.spec_from_file_location("collect_hits", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class RecurrenceChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.rows = []
        for cid, unit, study in [("t1", "model_a", "study_a"), ("t2", "model_a", "study_a"), ("t3", "model_b", "study_b")]:
            self.rows.append({"comparison_id": cid, "study_id": study, "model_id": unit,
                "gene_summary_file": f"{cid}.tsv", "include_in_recurrence": "yes",
                "recurrence_unit_id": unit, "review_status": "verified", "decision_note": "Synthetic fixture"})
            values = [("POS", ".01", "1", ".9"), ("NEG_ONLY", ".8", "-2", ".001"),
                      ("EDGE", ".05", "2", ".9"), ("DISCORDANT", ".01", "-1", ".9")]
            if cid in {"t1", "t2"}:
                values.append(("ONE_UNIT", ".01", "1", ".9"))
            with (self.root / f"{cid}.tsv").open("w", newline="") as handle:
                w = csv.writer(handle, delimiter="\t")
                w.writerow(["id", "pos|fdr", "pos|lfc", "neg|fdr"])
                w.writerows(values)

    def manifest(self):
        path = self.root / "manifest.tsv"
        with path.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(self.rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(self.rows)
        return path

    def run_analysis(self, policy="report-only"):
        return MODULE.collect(self.manifest(), self.root, self.root / "out", .05, 2, policy)

    def test_positive_selection_and_unique_units(self):
        report = self.run_analysis()
        _, genes = MODULE.read_table(self.root / "out/recurrent_genes.tsv")
        self.assertEqual({x["gene_id"] for x in genes}, {"POS", "DISCORDANT"})
        self.assertEqual(report["n_recurrent_genes"], 2)
        pos = next(x for x in genes if x["gene_id"] == "POS")
        self.assertEqual(pos["n_recurrence_units"], "2")
        self.assertEqual(pos["n_supporting_comparisons"], "3")

    def test_explicit_lfc_policy_changes_selection(self):
        self.run_analysis("require-positive")
        self.assertEqual((self.root / "out/recurrent_genes.txt").read_text(), "POS\n")

    def test_pending_decision_stops_before_outputs(self):
        self.rows[0]["include_in_recurrence"] = "PENDING"
        with self.assertRaisesRegex(ValueError, "Pending"):
            self.run_analysis()
        self.assertFalse((self.root / "out").exists())

    def test_missing_dataset_does_not_silently_disappear(self):
        (self.root / "t3.tsv").unlink()
        with self.assertRaises(FileNotFoundError):
            self.run_analysis()
        self.assertFalse((self.root / "out").exists())

    def test_duplicate_result_file_cannot_count_twice(self):
        self.rows[2]["gene_summary_file"] = "t1.tsv"
        with self.assertRaisesRegex(ValueError, "same gene-summary"):
            self.run_analysis()
        self.assertFalse((self.root / "out").exists())


if __name__ == "__main__":
    unittest.main()
