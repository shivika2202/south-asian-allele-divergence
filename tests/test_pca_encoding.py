import numpy as np

from conftest import load_script

pca = load_script("06_pca_scikit_allel.py")

# Genotypes are (variants x samples x ploidy); -1 means missing or padding.
# Sample 0 is female (diploid), sample 1 is male (hemizygous on chrX).
GT = np.array([
    [[0, 1], [1, -1]],   # variant 0: female het, male carries ALT
    [[1, 1], [0, -1]],   # variant 1: female hom ALT, male carries REF
    [[-1, -1], [1, -1]],  # variant 2: female missing
])


def test_autosomal_encoding_sums_both_alleles():
    counts = pca.to_allele_counts(GT[:2, :1, :], is_x_linked=False)
    assert counts[:, 0].tolist() == [1, 2]


def test_male_x_calls_are_doubled_onto_the_diploid_scale():
    counts = pca.to_allele_counts(GT, is_x_linked=True, is_male=np.array([False, True]))
    assert counts[:, 1].tolist() == [2, 0, 2]   # one ALT copy in a male = "2", same scale as hom ALT
    assert counts[:, 0].tolist() == [1, 2, -1]  # female unchanged; missing stays -1


def test_haploid_padding_is_not_mistaken_for_missing_in_males():
    # Without sex information, a male's padded second allele (-1) would read as missing.
    counts = pca.to_allele_counts(GT, is_x_linked=True, is_male=None)
    assert counts[:, 1].tolist() == [-1, -1, -1]


def test_drop_unusable_variants_removes_missing_and_pooled_monomorphic_sites():
    ac = np.array([
        [0, 1, 2, 1],    # usable
        [2, 2, 2, 2],    # monomorphic once pooled: would break Patterson scaling
        [0, -1, 1, 2],   # has a missing call
        [1, 0, 0, 1],    # usable
    ])
    kept, n_missing, n_mono = pca.drop_unusable_variants(ac)
    assert kept.tolist() == [[0, 1, 2, 1], [1, 0, 0, 1]]
    assert (n_missing, n_mono) == (1, 1)
