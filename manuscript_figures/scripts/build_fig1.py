#!/usr/bin/env python3
"""Build final Fig. 1 (schematic) for Nature Communications.

Starts from the macOS/PowerPoint vector export (Figure_1_v2.1.pdf) and
  1. swaps the 9 rasterised icons for their original SVG vectors (from the PPTX),
  2. redraws the hypercube as vector lines fitted to its native 1024 px bitmap
     (fig1_hypercube.py; --hypercube bitmap keeps the native-resolution bitmap instead),
  3. replaces the one Cambria Math glyph with a vector arrow,
  4. scales the page to 180 mm width (text 5.0 pt -> 6.7 pt).
"""
import argparse
import math
import re
import sys
import zipfile
from pathlib import Path

import fitz  # PyMuPDF

MM = 72 / 25.4
# PDF image xref -> PPTX media SVG (same order as the <p:pic> elements on slide 3)
SVG_FOR_XREF = {10: "image2.svg", 13: "image4.svg", 14: "image6.svg", 19: "image8.svg",
                24: "image10.svg", 25: "image12.svg", 26: "image14.svg",
                27: "image16.svg", 30: "image18.svg"}
HYPERCUBE_XREF, HYPERCUBE_MEDIA = 31, "image19.png"


def drop_resource(doc, page_xref, kind, name):
    """Delete /name from the page's Resources/<kind> dictionary (no null placeholders)."""
    rtyp, rval = doc.xref_get_key(page_xref, "Resources")
    # Resources may be indirect (original file) or inline (after PyMuPDF edits)
    holder, key = (int(rval.split()[0]), kind) if rtyp == "xref" else (page_xref, f"Resources/{kind}")
    typ, val = doc.xref_get_key(holder, key)
    pat = rf"/{name}\s+\d+\s+\d+\s+R"
    if typ == "xref":  # indirect sub-dictionary
        dxref = int(val.split()[0])
        doc.update_object(dxref, re.sub(pat, "", doc.xref_object(dxref, compressed=True)))
    elif typ == "dict":
        doc.xref_set_key(holder, key, re.sub(pat, "", val))
    else:
        sys.exit(f"ERROR: unexpected Resources/{kind} type {typ}")
    if doc.xref_get_key(holder, key)[0] not in ("dict", "xref"):
        sys.exit(f"ERROR: Resources/{kind} corrupted while removing /{name}")


def replace_cambria_arrow(page, verbose=False):
    """Remove the single '↻' set in CambriaMath and redraw it as a vector path."""
    spans = [s for b in page.get_text("dict")["blocks"] for ln in b.get("lines", [])
             for s in ln["spans"] if "Cambria" in s["font"]]
    if len(spans) != 1 or spans[0]["text"].strip() != "↻":
        sys.exit(f"ERROR: expected one CambriaMath '↻' span, found {[s['text'] for s in spans]}")
    s = spans[0]
    bb = fitz.Rect(s["bbox"])
    rgb = tuple(((s["color"] >> k) & 255) / 255 for k in (16, 8, 0))
    page.add_redact_annot(bb + (0.15, 0.9, -0.15, -0.2), fill=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    # glyph geometry measured from a 1600 dpi render of the original
    cx, cy, r, lw = bb.x0 + 0.47 * bb.width, s["origin"][1] - 1.63, 1.42, 0.36
    phis = [math.radians(-50 + i * 270 / 48) for i in range(49)]  # y-down: increasing = clockwise
    pts = [fitz.Point(cx + r * math.cos(p), cy + r * math.sin(p)) for p in phis]
    sh = page.new_shape()
    sh.draw_polyline(pts)
    end, phi = pts[-1], phis[-1]
    tx, ty = -math.sin(phi), math.cos(phi)  # travel direction at the arrow tip
    for rot in (math.radians(150), math.radians(-150)):
        bx = tx * math.cos(rot) - ty * math.sin(rot)
        by = tx * math.sin(rot) + ty * math.cos(rot)
        sh.draw_line(end, fitz.Point(end.x + 0.85 * bx, end.y + 0.85 * by))
    sh.finish(color=rgb, width=lw, lineCap=1, lineJoin=1, closePath=False)
    sh.commit()
    # drop the now-unused font resource so CambriaMath is no longer embedded
    for f in page.get_fonts(full=True):
        if "Cambria" in f[3]:
            drop_resource(page.parent, page.xref, "Font", f[4])
    if verbose:
        print(f"  replaced CambriaMath '↻' at {bb} with vector arrow", file=sys.stderr)


def draw_hypercube(page, box, png, verbose=False):
    """Draw the fitted hypercube (image-pixel geometry) into the bitmap's placement box."""
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fig1_hypercube

    g = fig1_hypercube.geometry(png)
    sx, sy = box.width / g["size"], box.height / g["size"]
    s = (sx + sy) / 2

    def P(xy):
        return fitz.Point(box.x0 + xy[0] * sx, box.y0 + xy[1] * sy)

    sh = page.new_shape()
    for p, q in g["solid"]:
        sh.draw_line(P(p), P(q))
    sh.finish(color=g["color"], width=g["width"] * s, lineCap=1, lineJoin=1, closePath=False)
    for p, q, pitch, phase in g["dotted"]:  # zero-length dashes + round caps = dots
        u = (q - p) / math.hypot(*(q - p))
        sh.draw_line(P(p + phase * u), P(q))
        sh.finish(color=g["color"], width=g["dot_d"] * s, lineCap=1, dashes=f"[0 {pitch * s:.4f}] 0",
                  closePath=False)
    sh.commit()
    if verbose:
        print(f"  hypercube: {len(g['solid'])} solid + {len(g['dotted'])} dotted edges, "
              f"stroke {g['width'] * s:.3f} pt, dots {g['dot_d'] * s:.3f} pt", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", default="/data2/chris/manuscript/fig1/Figure_1_v2.1.pdf", help="vector PDF export of slide 3")
    ap.add_argument("--pptx", default="/data2/chris/manuscript/fig1/PrinTE_Figure_1_v0_THedited_SO.pptx", help="source PPTX (for SVG/PNG media)")
    ap.add_argument("--width-mm", type=float, default=180.0, help="final page width in mm (default 180)")
    ap.add_argument("--hypercube", choices=["vector", "bitmap"], default="vector",
                    help="vector = lines fitted to the 1024 px bitmap (default); bitmap = native 1024 px image")
    ap.add_argument("-o", "--out", default="fig1/fig1.pdf", help="output PDF")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    a = ap.parse_args()

    src = fitz.open(a.pdf)
    page = src[0]
    zf = zipfile.ZipFile(a.pptx)
    boxes = {i["xref"]: fitz.Rect(i["bbox"]) for i in page.get_image_info(xrefs=True)}
    missing = (set(SVG_FOR_XREF) | {HYPERCUBE_XREF}) - set(boxes)
    if missing:
        sys.exit(f"ERROR: expected image xrefs not found in {a.pdf}: {sorted(missing)}")

    # 1. raster icons -> vector SVG: drop each "/ImN Do" paint operator and its
    #    resource entry (delete_image() would leave 1x1 placeholders), then overlay the SVG.
    names = {im[0]: im[7] for im in page.get_images(full=True)}
    cxref = page.get_contents()[0]
    content = src.xref_stream(cxref)
    to_remove = list(SVG_FOR_XREF) + ([HYPERCUBE_XREF] if a.hypercube == "vector" else [])
    for xref in to_remove:
        op = re.compile(rb"/" + names[xref].encode() + rb"\s+Do\b")
        content, n = op.subn(b"", content)
        if n != 1:
            sys.exit(f"ERROR: expected exactly one '/{names[xref]} Do' in page content, found {n}")
        drop_resource(src, page.xref, "XObject", names[xref])
    src.update_stream(cxref, content)
    for xref, svg in SVG_FOR_XREF.items():
        svgdoc = fitz.open("svg", zf.read(f"ppt/media/{svg}"))
        vec = fitz.open("pdf", svgdoc.convert_to_pdf())
        page.show_pdf_page(boxes[xref], vec, 0, keep_proportion=False)
        if a.verbose:
            print(f"  xref {xref} -> {svg} at {boxes[xref]}", file=sys.stderr)

    # 2. hypercube
    png = zf.read(f"ppt/media/{HYPERCUBE_MEDIA}")
    if a.hypercube == "vector":
        draw_hypercube(page, boxes[HYPERCUBE_XREF], png, a.verbose)
    else:  # native-resolution bitmap (placement + clip are kept by replace_image)
        page.replace_image(HYPERCUBE_XREF, stream=png)

    # 3. the lone Cambria Math glyph (serif font) -> vector clockwise open-circle arrow
    replace_cambria_arrow(page, a.verbose)

    # 4. scale to final width
    w = a.width_mm * MM
    h = page.rect.height * w / page.rect.width
    out = fitz.open()
    out.new_page(width=w, height=h).show_pdf_page(fitz.Rect(0, 0, w, h), src, 0)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.save(a.out, garbage=4, deflate=True)
    print(f"wrote {a.out}: {a.width_mm:.0f} x {h / MM:.1f} mm, scale {w / src[0].rect.width:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
