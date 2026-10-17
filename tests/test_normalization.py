"""Regression checks for the MAGeCK 0.5.9.5 median-to-total fallback."""

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    path = ROOT / "scripts" / "mle" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCRIPTS = [load_script(name) for name in (
    "run_lau_mle.py", "run_lin_mle_joint.py",
    "run_tuano_mle.py", "run_clements_mle.py",
)]


def qc(zero_count, selected=True):
    return {"sample": "OL_1", "guides": 100, "zero_count_guides": zero_count,
            "zero_percent": zero_count, "selected": selected}


class MedianNormalizationTests(unittest.TestCase):
    def test_fallback_boundary_and_unused_columns(self):
        for module in SCRIPTS:
            with self.subTest(script=module.__name__):
                # At 45% the engine warns but retains median; above it switches.
                module.require_median_counts([qc(45), qc(90, selected=False)])
                with self.assertRaisesRegex(ValueError, "switches to total"):
                    module.require_median_counts([qc(46)])

    def test_log_verification_distinguishes_warning_from_fallback(self):
        cases = (
            ("INFO: Final size factor: 1.1 0.9\n", True),
            ("WARNING: Consider switching to total read count normalization.\n"
             "INFO: Final size factor: 1.1 0.9\n", True),
            ("WARNING: Switch to total read count normalization.\n"
             "INFO: Final size factor: 1.1 0.9\n", False),
            ("INFO: Parameters only; no normalization completed.\n", False),
        )
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder) / "fit"
            for module in SCRIPTS:
                for log_text, accepted in cases:
                    with self.subTest(script=module.__name__, accepted=accepted,
                                      log=log_text):
                        Path(str(prefix) + ".log").write_text(log_text)
                        if accepted:
                            module.verify_median_log(prefix)
                        else:
                            with self.assertRaises(RuntimeError):
                                module.verify_median_log(prefix)
                        report = json.loads(Path(str(prefix) + ".normalization.json").read_text())
                        self.assertEqual(report["median_verified_from_log"], accepted)

    def test_missing_log_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as folder:
            for module in SCRIPTS:
                with self.subTest(script=module.__name__):
                    with self.assertRaises(FileNotFoundError):
                        module.verify_median_log(Path(folder) / "missing")

    def test_lau_sparse_dry_run_stops_and_preserves_inputs(self):
        module = SCRIPTS[0]
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            counts = root / "lau_counts.count.txt"
            original = "sgRNA\tGene\tControl\tOL_1\tOL_2\n" + "".join(
                f"g{i}\tG{i // 5}\t100\t{0 if i < 82 else 100}\t100\n"
                for i in range(100)
            )
            counts.write_text(original)
            with patch.object(module, "FOLDER", root), patch("sys.argv", ["run_lau_mle.py", "--dry-run"]), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, "No MLE fit was run"):
                    module.main()
            self.assertEqual(counts.read_text(), original)
            selected = list(root.glob("MAGeCK_MLE_results/*/selected.count.txt"))
            self.assertEqual(len(selected), 1)
            self.assertEqual(len(selected[0].read_text().splitlines()), 101)
            self.assertFalse(list(root.rglob("*.gene_summary.txt")))

    def test_rra_log_guard(self):
        helper = ROOT / "scripts" / "mageck" / "median_normalization.sh"
        cases = (
            ("INFO: Final size factor: 1 1\n", True),
            ("Consider switching to total read count normalization.\nFinal size factor: 1 1\n", True),
            ("Switch to total read count normalization.\nFinal size factor: 1 1\n", False),
            ("Parameters only\n", False),
            (None, False),
        )
        with tempfile.TemporaryDirectory() as folder:
            for i, (content, accepted) in enumerate(cases):
                log = Path(folder) / f"case_{i}.log"
                if content is not None:
                    log.write_text(content)
                result = subprocess.run(
                    ["bash", "-c", 'source "$1"; require_median_log "$2"',
                     "test", str(helper), str(log)], capture_output=True, text=True
                )
                self.assertEqual(result.returncode == 0, accepted, result.stderr)


if __name__ == "__main__":
    unittest.main()
