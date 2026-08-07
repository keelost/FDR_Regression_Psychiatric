import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ProjectConfigTests(unittest.TestCase):
    def test_target_table_has_16_unique_targets(self):
        with (ROOT / "config" / "targets.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(16, len(rows))
        self.assertEqual(16, len({row["target"] for row in rows}))
        self.assertEqual({"target", "traits_no", "traits_with"}, set(rows[0]))

    def test_overlap_groups_do_not_intersect(self):
        with (ROOT / "config" / "targets.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            no = {x for x in row["traits_no"].split(",") if x}
            with_overlap = {x for x in row["traits_with"].split(",") if x}
            self.assertFalse(no & with_overlap, row["target"])

    def test_version_directories_are_separate(self):
        with (ROOT / "config" / "versions.tsv").open(encoding="utf-8", newline="") as handle:
            rows = {row["version"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual({"v7", "v8"}, set(rows))
        for column in ("metaxcan_input", "metaxcan_fdrreg", "smultixcan_input", "smultixcan_fdrreg"):
            self.assertNotEqual(rows["v7"][column], rows["v8"][column])
            self.assertTrue(rows["v7"][column].endswith("_v7"))
            self.assertFalse(rows["v8"][column].endswith("_v7"))

    def test_regions_has_13_unique_entries(self):
        regions = [x.strip() for x in (ROOT / "config" / "regions.txt").read_text().splitlines() if x.strip()]
        self.assertEqual(13, len(regions))
        self.assertEqual(13, len(set(regions)))

    def test_v7_metaxcan_script_uses_v7_directories(self):
        script = (ROOT / "scripts" / "real" / "v7" / "metaxcan_fdrreg.R").read_text()
        self.assertIn('METAXCAN_SUBDIR <- "06.metaxcan_v7"', script)
        self.assertIn('"07.metaxcan_fdrreg_v7"', script)

    def test_imputation_reference_table_covers_configured_traits(self):
        with (ROOT / "config" / "targets.tsv").open(encoding="utf-8", newline="") as handle:
            target_rows = list(csv.DictReader(handle, delimiter="\t"))
        configured = set()
        for row in target_rows:
            configured.add(row["target"])
            configured.update(x for x in row["traits_no"].split(",") if x)
            configured.update(x for x in row["traits_with"].split(",") if x)

        with (ROOT / "config" / "imputation_refs.tsv").open(encoding="utf-8", newline="") as handle:
            reference_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(configured, {row["trait"] for row in reference_rows})
        self.assertEqual({"EUR", "EAS"}, {row["panel"] for row in reference_rows})

    def test_frozen_annotation_headers(self):
        annotation_dir = ROOT / "data" / "annotations"
        expected = {
            "magma-library-all.csv": {"ID", "ENSEMBL_GENE_ID", "GeneName"},
            "magma-library-uniq-entrez.csv": {"ID", "GeneName"},
            "magma-library-uniq-ensembl.csv": {"ENSEMBL_GENE_ID"},
        }
        for name, required in expected.items():
            with (annotation_dir / name).open(encoding="utf-8-sig", newline="") as handle:
                header = set(next(csv.reader(handle)))
            self.assertTrue(required <= header, name)


if __name__ == "__main__":
    unittest.main()
