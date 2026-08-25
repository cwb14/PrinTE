"""PrinTE BED: chrom, start, end, feature_ID, TSD, strand. 0-based, half-open."""

import pytest

pytestmark = pytest.mark.slow


def chrom_lengths(fasta):
    lengths, name = {}, None
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        else:
            lengths[name] += len(line.strip())
    return lengths


def bed_records(bed):
    for line in bed.read_text().splitlines():
        if line.strip():
            f = line.split("\t")
            assert len(f) == 6, f"expected 6 columns, got {len(f)}: {f}"
            yield f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]


def test_records_lie_inside_their_chromosome(burnin_run):
    lengths = chrom_lengths(burnin_run / "burnin.fasta")
    seen = 0
    for chrom, start, end, fid, _tsd, _strand in bed_records(burnin_run / "burnin.bed"):
        assert chrom in lengths, f"{chrom} absent from the genome"
        assert 0 <= start < end <= lengths[chrom], f"{fid}: {start}-{end} outside {chrom}"
        seen += 1
    assert seen > 0


def test_strand_and_tsd_are_well_formed(burnin_run):
    for _c, _s, _e, _fid, tsd, strand in bed_records(burnin_run / "burnin.bed"):
        assert strand in ("+", "-")
        assert tsd == "NA" or set(tsd.upper()) <= set("ACGTN"), tsd


def test_genes_and_tes_are_distinguishable(burnin_run):
    genes = tes = 0
    for _c, _s, _e, fid, _t, _st in bed_records(burnin_run / "burnin.bed"):
        if "#" in fid:
            tes += 1
        else:
            assert fid.startswith("gene"), fid
            genes += 1
    assert genes > 0 and tes > 0


def test_intact_ltr_elements_carry_an_ltrlen_tag(burnin_run):
    for _c, _s, _e, fid, _t, _st in bed_records(burnin_run / "burnin.bed"):
        if "#LTR/" in fid and "_FRAG" not in fid and "_SOLO" not in fid:
            assert "~LTRlen:" in fid, f"intact LTR-RT with no LTRlen: {fid}"


def test_features_do_not_overlap(burnin_run):
    by_chrom = {}
    for chrom, start, end, fid, _t, _s in bed_records(burnin_run / "burnin.bed"):
        by_chrom.setdefault(chrom, []).append((start, end, fid))
    for chrom, feats in by_chrom.items():
        feats.sort()
        for (s1, e1, f1), (s2, _e2, f2) in zip(feats, feats[1:]):
            assert e1 <= s2, f"{chrom}: {f1} ({s1}-{e1}) overlaps {f2} (starts {s2})"
