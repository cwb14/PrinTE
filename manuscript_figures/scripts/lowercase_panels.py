#!/usr/bin/env python3
"""Rewrite the panel letters of a finished figure PDF in lower case (Nature Communications style).

Finds every single-letter bold text span (the panel letters), deletes it and re-typesets the
lower-case letter in Arial Bold at the same baseline, size and colour. Writes <stem><suffix>.pdf.
"""
import argparse
import sys
from pathlib import Path

import fitz

ARIAL_BOLD = "/home/chris/mamba/pkgs/mscorefonts-0.0.1-3/fonts/arialbd.ttf"


def panel_letters(page, min_size=6.0):
    """Single upper-case letters set in a bold font = the panel labels."""
    out = []
    for b in page.get_text("dict")["blocks"]:
        for ln in b.get("lines", []):
            for s in ln["spans"]:
                t = s["text"].strip()
                if len(t) == 1 and t.isalpha() and t.isupper() and "bold" in s["font"].lower() \
                        and s["size"] >= min_size:
                    out.append(s)
    return out


def convert(path, suffix, scale=1.4, dy=0.0, verbose=False):
    doc = fitz.open(path)
    if doc.page_count != 1:
        sys.exit(f"ERROR: {path} has {doc.page_count} pages")
    page = doc[0]
    spans = panel_letters(page)
    if not spans:
        sys.exit(f"ERROR: no panel letters found in {path}")
    # delete the old letters first: a redaction would also remove text inserted beforehand
    for s in spans:
        page.add_redact_annot(fitz.Rect(s["bbox"]), fill=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    font = fitz.Font(fontfile=ARIAL_BOLD)
    for s in spans:
        rgb = tuple(((s["color"] >> k) & 255) / 255 for k in (16, 8, 0))
        letter = s["text"].strip().lower()
        size = s["size"] * scale  # lower case x-height < upper case cap height; scale to match
        x, y = s["origin"]
        y -= dy  # lift the whole set (e.g. so a descender clears the panel below)
        # keep the enlarged letter on the page (ink top of a lower-case letter is its x-height)
        top = y - size * font.ascender * 0.75
        if top < 1:
            y += 1 - top
        page.insert_text(fitz.Point(x, y), letter, fontsize=size, fontfile=ARIAL_BOLD,
                         fontname="ArialBd", color=rgb, render_mode=0)
        if verbose:
            print(f"  {s['text'].strip()} -> {letter} at ({x:.1f}, {y:.1f}) "
                  f"{s['size']:.2f} -> {size:.2f} pt", file=sys.stderr)
    out = Path(path).with_name(Path(path).stem + suffix + ".pdf")
    doc.subset_fonts()  # keep only the glyphs used; the full Arial Bold adds ~150 KB per figure
    doc.save(out, garbage=4, deflate=True)

    check = fitz.open(out)[0]
    got = sorted(s["text"].strip() for s in panel_letters(check))
    lower = sorted(s["text"].strip() for b in check.get_text("dict")["blocks"] for ln in b.get("lines", [])
                   for s in ln["spans"] if len(s["text"].strip()) == 1 and s["text"].strip().islower()
                   and "bold" in s["font"].lower())
    want = sorted(s["text"].strip().lower() for s in spans)
    if got or lower != want:
        sys.exit(f"ERROR: {out}: expected {want}, found lower={lower}, upper left={got}")
    print(f"wrote {out}: {len(spans)} panel letters -> {''.join(want)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdfs", nargs="+", help="figure PDFs")
    ap.add_argument("--suffix", default="_fixed", help="output name suffix (default: %(default)s)")
    ap.add_argument("--scale", type=float, default=1.4,
                    help="panel-letter size relative to the original upper-case letter (default: %(default)s)")
    ap.add_argument("--dy", type=float, default=0.0,
                    help="raise every panel letter by this many points (default: %(default)s)")
    ap.add_argument("-v", "--verbose", action="store_true", help="list every letter replaced")
    a = ap.parse_args()
    for p in a.pdfs:
        convert(p, a.suffix, a.scale, a.dy, a.verbose)


if __name__ == "__main__":
    main()
