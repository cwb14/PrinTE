from printe import ltr_dens

# Two rows as legacy Kmer2LTR writes them (its README's example): 16 tab-separated
# columns, no header. The last four (trims and LTR coordinates) came later than the
# twelve ltr_dens names.
LEGACY_ROWS = (
    "Gypsy1#LTR_Ty3\t584\t574\t51\t40\t11\t0.088850\t1480833\t0.094572\t1576200"
    "\t0.096055\t1600917\t5\t5\t584\t8210\n"
    "Copia3#LTR_Ty1\t301\t298\t19\t15\t4\t0.063758\t1062633\t0.066666\t1111100"
    "\t0.067383\t1123050\t4\t6\t301\t4502\n"
)


def test_legacy_kmer2ltr_rows_are_read_by_position(tmp_path):
    f = tmp_path / "gen100_LTR.tsv"
    f.write_text(LEGACY_ROWS)
    df = ltr_dens.read_ltr_tsv(str(f))
    assert list(df["LTR_RT_name"]) == ["Gypsy1#LTR_Ty3", "Copia3#LTR_Ty1"]
    assert list(df["K2P_d"]) == [0.096055, 0.067383]
    assert list(df["K2P_T"]) == [1600917, 1123050]
