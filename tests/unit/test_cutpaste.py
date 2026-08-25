from printe import cutpaste_common as cc


def test_cutpaste_families_read_from_the_transposition_column(tmp_path):
    r = tmp_path / "r.tsv"
    r.write_text(
        "# class\tsuperfamily\tintact\tfrag\ttransposition\n"
        "DNA\tDTA\t0.5\t0.5\tcutpaste\n"
        "LTR\tCopia\t0.3\t0.3\tcopypaste\n"
        "LTR\tGypsy\t0.2\t0.2\n"
    )
    fams = cc.load_cutpaste_set(str(r))
    assert ("DNA", "DTA") in fams
    assert ("LTR", "Copia") not in fams
    assert ("LTR", "Gypsy") not in fams  # a missing column means copy-and-paste


def test_cutpaste_flag_is_case_insensitive(tmp_path):
    r = tmp_path / "r.tsv"
    r.write_text("MITE\tDTH\t0.5\t0.5\tCutPaste\n")
    assert ("MITE", "DTH") in cc.load_cutpaste_set(str(r))


def test_missing_ratios_file_means_everything_is_copypaste():
    assert cc.load_cutpaste_set("/nonexistent/ratios.tsv") == set()
    assert cc.load_cutpaste_set(None) == set()


def test_parse_family_ignores_ltrlen_frag_and_solo_suffixes():
    assert cc.parse_family("elem1#LTR/Copia") == ("LTR", "Copia")
    assert cc.parse_family("elem1#LTR/Copia~LTRlen:430") == ("LTR", "Copia")
    assert cc.parse_family("elem1#LTR/Copia~LTRlen:430_SOLO") == ("LTR", "Copia")
    assert cc.parse_family("TE_00016444_FRAG#MITE/DTA") == ("MITE", "DTA")


def test_parse_family_on_a_gene_returns_unknown():
    assert cc.parse_family("gene1") == ("unknown", "unknown")


def test_nested_annotation_does_not_confuse_the_family():
    fid = "tuteh_AC183372_584#LTR/unknown~LTRlen:126;CUT_BY:Os2721#DNAnona/hAT"
    assert cc.parse_family(fid) == ("LTR", "unknown")


def test_debt_round_trips(tmp_path):
    p = tmp_path / "debt.tsv"
    cc.write_debt(str(p), {("DNA", "DTA"): 3, ("MITE", "DTH"): 1})
    assert cc.read_debt(str(p)) == {("DNA", "DTA"): 3, ("MITE", "DTH"): 1}
