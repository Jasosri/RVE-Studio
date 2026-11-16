"""Generate the dashboard illustrations in rvestudio/gui/assets (run from the repository root).

Every picture is drawn from the software itself: the ply stack of configuration 1, the hexagonal RVE, the
built-in periodic shell FE of configuration 1 under in-plane shear, a seeded parametric database with its
Pareto front, and a wall buckling mode.
"""
from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "rvestudio" / "gui" / "assets"
OUT.mkdir(parents=True, exist_ok=True)
W, H, DPI = 4.8, 3.0, 150          # 720 x 450 px
BG = "#FFFFFF"
PALETTE = ["#2F80ED", "#F2994A", "#1F9D72", "#7A5AF8", "#D64545", "#56CCF2"]


def save(fig, name):
    fig.savefig(OUT / name, dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)
    print("wrote", OUT / name)


# ------------------------------------------------------------------------------------------- laminate
def laminate():
    fig = plt.figure(figsize=(W, H), facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 10); ax.set_ylim(0, 6.25); ax.axis("off")
    angles = [45, -45, 0, -45, 45]
    cols = ["#2F80ED", "#56CCF2", "#F2C94C", "#56CCF2", "#2F80ED"]
    for k, (ang, col) in enumerate(zip(angles, cols)):
        y0 = 0.6 + k * 1.05
        dx, dy = 1.6, 0.55                          # oblique depth of the ply
        P = np.array([[1.2, y0], [7.2, y0], [7.2 + dx, y0 + dy], [1.2 + dx, y0 + dy]])
        ax.add_patch(Polygon(P, closed=True, fc=col, ec="#1B2D49", lw=0.8, alpha=0.95))
        ax.add_patch(Polygon([[1.2, y0 - 0.16], [7.2, y0 - 0.16], [7.2, y0], [1.2, y0]], fc="#1B2D49", alpha=0.35, lw=0))
        # fibres: lines in the ply plane at the ply angle, clipped to the parallelogram
        th = math.radians(ang)
        d = np.array([math.cos(th), math.sin(th)])
        for s in np.linspace(-8, 8, 34):
            pts = []
            for t in np.linspace(-10, 10, 400):
                u, v = 4.0 + s * (-d[1]) + t * d[0], 0.5 + s * d[0] * 0.18 + t * d[1] * 0.18
                if 0 <= v <= 1 and 0 <= u - 1.2 - v * dx <= 6.0 + 1e-9 and 1.2 + v * dx <= u <= 7.2 + v * dx:
                    pts.append((u, y0 + v * dy))
            if len(pts) > 2:
                pts = np.array(pts); ax.plot(pts[:, 0], pts[:, 1], color="white", lw=0.55, alpha=0.75)
        ax.text(9.15, y0 + 0.25, f"{ang:+d}°" if ang else "0°", fontsize=9, fontweight="bold", color="#172033", va="center")
    save(fig, "laminate.png")


# ------------------------------------------------------------------------------------------- 3D RVE
def rve3d():
    from rvestudio.examples_data import paper_hexagonal_case
    m, _, _ = paper_hexagonal_case(1)
    X = np.asarray(m.nodes, float); a1 = np.asarray(m.a1, float); a2 = np.asarray(m.a2, float); h = 60.0
    fig = plt.figure(figsize=(W, H), facecolor=BG)
    ax = fig.add_axes([-0.22, -0.32, 1.44, 1.42], projection="3d"); ax.set_facecolor(BG)
    polys, cols = [], []
    for i in range(-2, 3):
        for j in range(-1, 2):
            sh = i * a1 + j * a2
            for e in m.elements:
                p, q = X[e.node_i] + sh, X[e.node_j] + sh
                polys.append([[p[0], p[1], 0], [q[0], q[1], 0], [q[0], q[1], h], [p[0], p[1], h]])
                main = (i == 0 and j == 0)
                base = matplotlib.colors.to_rgb(PALETTE[list(m.elements).index(e) % 3])
                cols.append((*base, 0.95) if main else (*(0.55 * np.array(base) + 0.45), 0.55))
    pc = Poly3DCollection(polys, facecolors=cols, edgecolors=(0.1, 0.18, 0.29, 0.6), linewidths=0.4)
    ax.add_collection3d(pc)
    allp = np.array([p for poly in polys for p in poly])
    ax.set_xlim(allp[:, 0].min(), allp[:, 0].max()); ax.set_ylim(allp[:, 1].min(), allp[:, 1].max()); ax.set_zlim(0, h)
    ax.set_box_aspect((np.ptp(allp[:, 0]), np.ptp(allp[:, 1]), h * 1.4))
    ax.view_init(elev=52, azim=-75); ax.set_axis_off()
    save(fig, "rve3d.png")


# ------------------------------------------------------------------------------------------- deformed
def deformed():
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.shell_fe import run_pbc_shell_check
    m, mats, lams = paper_hexagonal_case(1)
    r = run_pbc_shell_check(m, mats, lams, mesh_size=10.0, return_fields=True)
    f = r["fields"]; Xc = np.asarray(f["coords"], float); E = np.asarray(f["elements"], int); U = np.asarray(f["u"]["12"], float)
    span = float(np.max(np.ptp(Xc[:, :2], axis=0))); scale = 0.07 * span / max(np.abs(U).max(), 1e-30)
    Xd = Xc + scale * U
    mag = np.linalg.norm(U, axis=1); fv = mag[E].mean(axis=1); fv = (fv - fv.min()) / max(np.ptp(fv), 1e-30)
    fig = plt.figure(figsize=(W, H), facecolor=BG)
    ax = fig.add_axes([0.0, -0.08, 1.0, 0.95], projection="3d")
    ax.add_collection3d(Poly3DCollection(Xc[E], facecolors=(0.6, 0.65, 0.72, 0.10), edgecolors=(0.4, 0.45, 0.5, 0.25), linewidths=0.2))
    ax.add_collection3d(Poly3DCollection(Xd[E], facecolors=matplotlib.colormaps["turbo"](0.08 + 0.9 * fv), edgecolors=(0, 0, 0, 0.12), linewidths=0.15))
    P = np.vstack([Xc, Xd]); ax.set_xlim(P[:, 0].min(), P[:, 0].max()); ax.set_ylim(P[:, 1].min(), P[:, 1].max()); ax.set_zlim(P[:, 2].min(), P[:, 2].max())
    ax.set_box_aspect((np.ptp(P[:, 0]), np.ptp(P[:, 1]), np.ptp(P[:, 2]) * 0.75)); ax.view_init(elev=30, azim=-60); ax.set_axis_off()
    save(fig, "deformed.png")


# ------------------------------------------------------------------------------------------- database, Pareto
def database_and_pareto():
    from rvestudio.database import generate_parametric_database, run_batch
    from rvestudio.postprocess import pareto_mask_min_x_max_y
    with tempfile.TemporaryDirectory() as d:
        generate_parametric_database(Path(d) / "in.csv", 160, seed=7, wall_shear="timoshenko")
        df = run_batch(Path(d) / "in.csv", Path(d) / "out.csv")
    df = df[df.Status == "OK"].copy()
    df["Qmin"] = np.minimum(df.Q11, df.Q22)
    geo = df.GeometryID.astype(str).str.upper()
    gcols = {g: PALETTE[k % len(PALETTE)] for k, g in enumerate(sorted(geo.unique()))}
    # database (Ashby-type) plot
    fig, ax = plt.subplots(figsize=(W, H), facecolor=BG)
    for g, c in gcols.items():
        s = df[geo == g]; ax.scatter(s.AverageDensityKgM3, s.E3, s=16, color=c, alpha=0.8, edgecolor="white", lw=0.4, label=g.title())
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("core density (kg/m³)", fontsize=8); ax.set_ylabel("E3 (MPa)", fontsize=8)
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(ScalarFormatter()); axis.set_minor_formatter(NullFormatter())
    ax.set_yticks([500, 1000, 2000, 5000]); ax.set_xticks([50, 100, 200, 300])
    ax.tick_params(labelsize=7); [ax.spines[k].set_visible(False) for k in ("top", "right")]
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="lower right")
    fig.tight_layout(); save(fig, "database.png")
    # multi-objective / Pareto
    x, y = df.AverageDensityKgM3.to_numpy(), df.Qmin.to_numpy()
    pm = pareto_mask_min_x_max_y(x, y)
    fig, ax = plt.subplots(figsize=(W, H), facecolor=BG)
    ax.scatter(x[~pm], y[~pm], s=14, color="#AFC3DB", edgecolor="white", lw=0.3, label="designs")
    o = np.argsort(x[pm]); ax.plot(x[pm][o], y[pm][o], "-", color="#D64545", lw=1.6)
    ax.scatter(x[pm], y[pm], s=30, color="#F2994A", edgecolor="#D64545", lw=0.8, zorder=3, label="Pareto front")
    best = np.argmax(y[pm] / x[pm]); ax.scatter([x[pm][best]], [y[pm][best]], marker="*", s=160, color="#7A5AF8", zorder=4, label="compromise")
    ax.set_yscale("log"); ax.set_xlabel("core density (kg/m³)  → minimise", fontsize=8); ax.set_ylabel("min(Q11, Q22) (MPa)  → maximise", fontsize=8)
    ax.tick_params(labelsize=7); [ax.spines[k].set_visible(False) for k in ("top", "right")]; ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout(); save(fig, "pareto.png")


# ------------------------------------------------------------------------------------------- buckling mode
def buckling():
    a, b = 50.0, 150.0
    x = np.linspace(0, a, 60); z = np.linspace(0, b, 150); Xg, Zg = np.meshgrid(x, z)
    # shear-type mode: two skewed half-waves (superposition of sine terms, as in the Ritz screen)
    wv = np.sin(np.pi * Xg / a) * (np.sin(2 * np.pi * Zg / b) + 0.55 * np.sin(3 * np.pi * Zg / b + 0.8 * np.pi * Xg / a))
    wv /= np.abs(wv).max()
    fig = plt.figure(figsize=(W, H), facecolor=BG)
    ax = fig.add_axes([-0.05, -0.12, 1.1, 1.05], projection="3d")
    ax.plot_surface(Xg, Zg, 9 * wv, cmap="coolwarm", rstride=2, cstride=2, linewidth=0.1, edgecolor=(0, 0, 0, 0.15), antialiased=True)
    ax.set_box_aspect((a, b, 30)); ax.view_init(elev=28, azim=-50); ax.set_axis_off()
    save(fig, "buckling.png")


# ------------------------------------------------------------------------------------------- banner
def banner():
    """Light honeycomb on a transparent background, for the dark dashboard banner."""
    from rvestudio.examples_data import paper_hexagonal_case
    m, _, _ = paper_hexagonal_case(1)
    X = np.asarray(m.nodes, float); a1 = np.asarray(m.a1, float); a2 = np.asarray(m.a2, float); h = 70.0
    fig = plt.figure(figsize=(4.8, 2.4))
    ax = fig.add_axes([-0.25, -0.45, 1.5, 1.7], projection="3d")
    ax.patch.set_alpha(0)
    polys, cols = [], []
    tones = [(0.34, 0.80, 0.95), (0.80, 0.92, 1.0), (0.56, 0.86, 0.98)]
    for i in range(-2, 3):
        for j in range(-1, 2):
            sh = i * a1 + j * a2
            for k, e in enumerate(m.elements):
                p, q = X[e.node_i] + sh, X[e.node_j] + sh
                polys.append([[p[0], p[1], 0], [q[0], q[1], 0], [q[0], q[1], h], [p[0], p[1], h]])
                cols.append((0.95, 0.60, 0.29, 0.95) if (i, j) == (0, 0) else (*tones[k % 3], 0.22))
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cols, edgecolors=(1, 1, 1, 0.55), linewidths=0.5))
    allp = np.array([p for poly in polys for p in poly])
    ax.set_xlim(allp[:, 0].min(), allp[:, 0].max()); ax.set_ylim(allp[:, 1].min(), allp[:, 1].max()); ax.set_zlim(0, h)
    ax.set_box_aspect((np.ptp(allp[:, 0]), np.ptp(allp[:, 1]), h * 1.3)); ax.view_init(elev=40, azim=-70); ax.set_axis_off()
    fig.savefig(OUT / "rve3d_banner.png", dpi=DPI, transparent=True); plt.close(fig); print("wrote", OUT / "rve3d_banner.png")


if __name__ == "__main__":
    banner(); laminate(); rve3d(); deformed(); database_and_pareto(); buckling()
