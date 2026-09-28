# Manuscript figures

`fig1.pdf` … `fig7.pdf` are the paper's seven main figures; `scripts/` drew them. Source
Data are published with the article.

| script | makes |
|---|---|
| `build_fig1.py` | Fig. 1, from the PowerPoint schematic (`fig1_hypercube.py` draws the hypercube) |
| `make_fig2.py` … `make_fig7.py` | Figs 2-7 (`make_fig3_trees.R` and `fig7_contmap.R` draw the trees in 3d-e and 7a) |
| `lowercase_panels.py` | last step: lower-case panel letters, fonts subset |
| `ncomms_style.py`, `srcdata.py` | shared journal style; per-panel Source Data tables |

The scripts read raw simulation output from absolute paths on our lab server, which is not
in this repository, so they record how each figure was made rather than rebuild it
elsewhere. Run in the directory that holds `scripts/`:

```bash
python scripts/build_fig1.py            # -> fig1/fig1.pdf; likewise make_fig2.py ... make_fig7.py
python scripts/lowercase_panels.py --scale 1.1 fig1/fig1.pdf        # -> fig1/fig1_fixed.pdf = fig1.pdf here
python scripts/lowercase_panels.py --scale 1.3 --dy 2 fig2/fig2.pdf
python scripts/lowercase_panels.py --scale 1.3 fig3/fig3.pdf fig4/fig4.pdf fig5/fig5.pdf fig6/fig6.pdf fig7/fig7.pdf
```

Python 3.10 (matplotlib 3.10, pandas, scipy, scikit-learn, scikit-image, seaborn, statsmodels,
Pillow, PyMuPDF 1.25, fontTools, openpyxl), R 4.4 (ape, phytools, ggtree, ggplot2, optparse),
and the Arial fonts.
