import pandas as pd
import pytest

from conftest import load_script

shared = load_script("02_find_shared_variants.py")
annot = load_script("03_annotate_variants.py")
freq = load_script("04_freq_comparison.py")


def sites(*positions):
    return pd.DataFrame({"CHROM": [10] * len(positions), "POS": list(positions)})


def test_intersect_keeps_only_sites_present_in_every_population():
    result = shared.intersect_sites([sites(1, 2, 3, 4), sites(2, 3, 4), sites(3, 4, 9), sites(4, 3)])
    assert sorted(result["POS"]) == [3, 4]


def test_sample_overlap_check_reports_double_assigned_samples():
    overlaps = shared.sample_overlaps({"GIH": {"a", "b"}, "PJL": {"c"}, "BEB": {"b", "d"}})
    assert overlaps == {("GIH", "BEB"): {"b"}}


def test_max_freq_spread_ignores_missing_populations():
    df = pd.DataFrame({"ID": ["v1", "v2"],
                       "GIH_ALT_FREQ": [0.296, 0.10], "PJL_ALT_FREQ": [0.234, None],
                       "BEB_ALT_FREQ": [0.238, 0.30], "ITU_ALT_FREQ": [0.201, 0.15]})
    result = freq.add_spread(df).set_index("ID")["max_freq_spread"]
    assert result["v1"] == pytest.approx(0.095)
    assert result["v2"] == pytest.approx(0.20)


def test_max_freq_spread_sorts_most_divergent_first():
    df = pd.DataFrame({"ID": ["small", "large"], "GIH_ALT_FREQ": [0.5, 0.1], "BEB_ALT_FREQ": [0.52, 0.4]})
    assert list(freq.add_spread(df)["ID"]) == ["large", "small"]


def write(path, text):
    path.write_text(text)
    return path


def test_detects_plink_and_vcftools_formats(tmp_path):
    assert annot.detect_filetype(write(tmp_path / "a.afreq", "#CHROM\tID\n")) == "plink"
    assert annot.detect_filetype(write(tmp_path / "a.frq", "CHROM\tPOS\n")) == "vcftools"
    with pytest.raises(ValueError):
        annot.detect_filetype(write(tmp_path / "x.txt", "something else\n"))


def test_vcftools_parser_handles_multiallelic_rows(tmp_path):
    path = write(tmp_path / "x.frq",
                 "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n"
                 "10\t100\t2\t200\tC:0.9\tA:0.1\n"
                 "10\t200\t3\t200\tA:0.7\tAT:0.2\tATT:0.1\n")
    df = annot.load_vcftools_freq(path).set_index("POS")
    assert df.loc[100, "ALT_ALLELE"] == "A" and df.loc[100, "ALT_FREQ"] == pytest.approx(0.1)
    assert df.loc[200, "N_ALLELES"] == 3
    assert df.loc[200, "ALL_ALLELES_RAW"] == "A:0.7;AT:0.2;ATT:0.1"


def test_annotation_matches_known_variants_by_position(tmp_path):
    afreq = write(tmp_path / "x.afreq",
                  "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
                  "10\t10:96541616:G:A\tG\tA\t0.35\t200\n"
                  "10\t10:96500000:C:T\tC\tT\t0.10\t200\n")
    lookup = write(tmp_path / "lookup.csv",
                   "CHROM,POS,rsID,Star_Allele,Gene,Function,Source\n"
                   "10,96541616,rs4244285,CYP2C19*2,CYP2C19,Loss-of-function,ClinVar\n")
    merged = annot.annotate(afreq, lookup, tmp_path / "out.tsv").set_index("POS")
    assert merged.loc[96541616, "Star_Allele"] == "CYP2C19*2"
    assert pd.isna(merged.loc[96500000, "rsID"])
