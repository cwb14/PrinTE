"""Vector reconstruction of the Fig. 1B hypercube icon from its 1024 px bitmap (PPTX media image19.png).

Solid edges: total-least-squares line fits to the skeleton of the solid stroke network; vertices are
the intersections of the fitted lines. Dotted edges: line fits to the centroids of the isolated dots,
redrawn as a regular dotted stroke (pitch/phase fitted to the detected dots). Coordinates are image
pixels (x right, y down).
"""
import io

import numpy as np
from PIL import Image
from scipy import ndimage
from skimage import measure, morphology

# approximate edges (seeds, image px) -> refined from the pixels
SOLID_SEEDS = {
    "obTop": ((383, 176), (826, 176)), "ofTop": ((208, 301), (654, 301)), "ibTop": ((461, 370), (700, 370)),
    "ifTop": ((383, 414), (608, 414)), "ifBot": ((383, 603), (577, 603)), "ofBot": ((208, 739), (654, 739)),
    "obRight": ((826, 176), (826, 605)), "ibRight": ((700, 370), (703, 545)), "ofRight": ((654, 301), (654, 739)),
    "x608": ((608, 414), (608, 589)), "ifRight": ((577, 414), (577, 603)), "ifLeft": ((383, 414), (383, 603)),
    "ofLeft": ((208, 301), (208, 739)), "d1": ((208, 301), (383, 414)), "d2": ((577, 603), (654, 739)),
    "d3": ((577, 414), (654, 301)), "d4": ((700, 370), (826, 176)), "d5": ((654, 739), (826, 605)),
    "d6": ((654, 301), (826, 176)), "d7": ((208, 301), (383, 176)), "d8": ((383, 414), (461, 370)),
    "d9": ((608, 414), (700, 370)), "d10": ((577, 603), (703, 545)),
}
DOTTED_SEEDS = {
    "obLeft": ((383, 176), (383, 414)), "k1": ((383, 176), (461, 370)), "ibLeft": ((461, 370), (461, 540)),
    "ibBot": ((461, 540), (703, 545)), "k2": ((208, 739), (383, 603)), "k3": ((461, 540), (383, 603)),
    "k4": ((703, 545), (826, 605)), "obBot": ((577, 603), (826, 605)),
}
# vertex = intersection of two fitted lines
VERTICES = {
    "obTL": ("obTop", "d7"), "obTR": ("obTop", "obRight"), "obBR": ("obRight", "d5"),
    "ofTL": ("ofTop", "ofLeft"), "ofTR": ("ofTop", "ofRight"), "ofBR": ("ofBot", "ofRight"),
    "ofBL": ("ofBot", "ofLeft"), "ibTL": ("ibTop", "d8"), "ibTR": ("ibTop", "ibRight"),
    "ibBR": ("ibRight", "d10"), "ifTL": ("ifTop", "ifLeft"), "ifTR": ("ifTop", "ifRight"),
    "ifTR2": ("ifTop", "x608"), "ifBR": ("ifBot", "ifRight"), "ifBL": ("ifBot", "ifLeft"),
    "X": ("x608", "d10"), "ibBL": ("ibLeft", "ibBot"),
    "ofBRd5": ("d5", "ofBot"),  # the back diagonal leaves the bottom edge ~9 px right of the corner
}
SOLID = {"obTop": ("obTL", "obTR"), "ofTop": ("ofTL", "ofTR"), "ibTop": ("ibTL", "ibTR"),
         "ifTop": ("ifTL", "ifTR2"), "ifBot": ("ifBL", "ifBR"), "ofBot": ("ofBL", "ofBRd5"),
         "obRight": ("obTR", "obBR"), "ibRight": ("ibTR", "ibBR"), "ofRight": ("ofTR", "ofBR"),
         "x608": ("ifTR2", "X"), "ifRight": ("ifTR", "ifBR"), "ifLeft": ("ifTL", "ifBL"),
         "ofLeft": ("ofTL", "ofBL"), "d1": ("ofTL", "ifTL"), "d2": ("ifBR", "ofBR"), "d3": ("ifTR", "ofTR"),
         "d4": ("ibTR", "obTR"), "d5": ("ofBRd5", "obBR"), "d6": ("ofTR", "obTR"), "d7": ("ofTL", "obTL"),
         "d8": ("ifTL", "ibTL"), "d9": ("ifTR2", "ibTR"), "d10": ("ifBR", "ibBR")}
DOTTED = {"obLeft": ("obTL", "ifTL"), "k1": ("obTL", "ibTL"), "ibLeft": ("ibTL", "ibBL"),
          "ibBot": ("ibBL", "ibBR"), "k2": ("ofBL", "ifBL"), "k3": ("ibBL", "ifBL"),
          "k4": ("ibBR", "obBR"), "obBot": ("ifBR", "obBR")}


def _fit(points):
    """Total-least-squares line: (centroid, unit direction)."""
    c = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - c)
    return c, vt[0]


def _near(points, p, q, tol, trim):
    p, q = np.asarray(p, float), np.asarray(q, float)
    L = np.linalg.norm(q - p)
    u = (q - p) / L
    n = np.array([-u[1], u[0]])
    t, d = (points - p) @ u, np.abs((points - p) @ n)
    return points[(d < tol) & (t > trim) & (t < L - trim)]


def _intersect(l1, l2):
    (c1, u1), (c2, u2) = l1, l2
    A = np.column_stack([u1, -u2])
    s = np.linalg.solve(A, c2 - c1)
    return c1 + s[0] * u1


def geometry(png_bytes):
    im = np.array(Image.open(io.BytesIO(png_bytes)).convert("RGBA"))
    alpha = im[..., 3] / 255.0
    mask = alpha > 0.5
    lab = measure.label(mask, connectivity=2)
    props = measure.regionprops(lab)
    solid_lbl = max(props, key=lambda r: r.area).label
    solid = lab == solid_lbl
    dots = np.array([r.centroid[::-1] for r in props if r.label != solid_lbl])  # (x, y)
    dot_d = float(np.median([r.equivalent_diameter for r in props if r.label != solid_lbl]))
    sk = np.argwhere(morphology.skeletonize(solid))[:, ::-1].astype(float)  # (x, y)
    # stroke width: 2 x median distance-transform at the skeleton overestimates by ~0.5 px (pixel centres);
    # 9.5 px matches the bitmap ink area to within 1%
    width = float(2 * np.median(ndimage.distance_transform_edt(solid)[morphology.skeletonize(solid)])) - 0.5
    color = tuple(np.median(im[..., :3][alpha > 0.9], axis=0) / 255.0)

    lines = {}
    for k, (p, q) in SOLID_SEEDS.items():
        # pass 1: wide capture around the approximate seed; pass 2: tight refit around the pass-1 line
        c, u = _fit(_near(sk, p, q, tol=9, trim=20))
        L = np.linalg.norm(np.subtract(q, p))
        t0 = (np.asarray(p, float) - c) @ u
        pts = _near(sk, c + t0 * u, c + (t0 + np.sign((np.asarray(q, float) - c) @ u - t0) * L) * u, tol=3, trim=16)
        if len(pts) < 10:
            raise ValueError(f"hypercube: too few skeleton pixels for edge {k}")
        lines[k] = _fit(pts)
    dotfit = {}
    for k, (p, q) in DOTTED_SEEDS.items():
        pts = _near(dots, p, q, tol=6, trim=-6)
        if len(pts) < 4:
            raise ValueError(f"hypercube: too few dots for edge {k}")
        lines[k] = _fit(pts)
        dotfit[k] = pts
    V = {v: _intersect(lines[a], lines[b]) for v, (a, b) in VERTICES.items()}

    dotted = []
    for k, (a, b) in DOTTED.items():
        # keep the dots on their own fitted line; only the extent comes from the solid-line vertices
        c, u = lines[k]
        p, q = c + ((V[a] - c) @ u) * u, c + ((V[b] - c) @ u) * u
        u = (q - p) / np.linalg.norm(q - p)
        t = np.sort((dotfit[k] - p) @ u)
        g = np.diff(t)
        p0 = float(np.median(g[(g > 10) & (g < 20)]))  # single-pitch gaps (dots merged into lines leave gaps)
        idx = np.round((t - t[0]) / p0)
        pitch, phase = np.polyfit(idx, t, 1)  # t_i = phase + pitch * k_i
        phase = float(np.mod(phase, pitch))
        dotted.append((p, q, float(pitch), phase))
    solid_segs = [(V[a], V[b]) for a, b in SOLID.values()]
    return dict(size=im.shape[1], color=color, width=width, dot_d=dot_d, solid=solid_segs, dotted=dotted,
                vertices=V)


def render_mask(g, n=1024):
    """Rasterise the vector geometry (for QC against the bitmap)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(n / 100, n / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, n)
    ax.set_ylim(n, 0)
    ax.axis("off")
    lw = g["width"] * 72 / 100  # px -> pt at dpi 100
    for p, q in g["solid"]:
        ax.plot([p[0], q[0]], [p[1], q[1]], color="k", lw=lw, solid_capstyle="round")
    for p, q, pitch, phase in g["dotted"]:
        L = np.linalg.norm(q - p)
        u = (q - p) / L
        ts = np.arange(phase, L + 1e-6, pitch)
        xy = p + np.outer(ts, u)
        ax.scatter(xy[:, 0], xy[:, 1], s=(g["dot_d"] * 72 / 100) ** 2, color="k", lw=0)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., 0] < 128
    plt.close(fig)
    return img
