"""Location of the RVE Studio user guide and theory guide (PDF), and how to cite the software."""
from __future__ import annotations

import sys
from pathlib import Path

GUIDES = {
    "user": ("RVE_Studio_User_Guide.pdf", "USER_GUIDE.md", "User Guide"),
    "theory": ("RVE_Studio_Theory_Guide.pdf", "THEORY.md", "Theory and Modelling Guidance"),
}

# The software builds on these papers; please cite them together (most recent first).
CITATIONS = (
    ("Sriharan J, et al. RVE Studio: open-source strain-energy homogenisation of laminated-composite cellular cores "
     "- verification of the full stiffness tensor, a bend-twist reliability limit and design tools. In preparation.", None),
    ("Sriharan J, Dias M, Adhikari S, Fernando D. A semi-analytical model for predicting the shear buckling of laminated "
     "composite honeycomb cores in sandwich panels. Composite Structures 351 (2025) 118629.",
     "https://doi.org/10.1016/j.compstruct.2024.118629"),
    ("Sriharan J, Dias M, Fernando D, Adhikari S. A unified approach for the prediction of the effective properties of "
     "laminated composite cellular core. Engineering Structures 299 (2024) 117106.",
     "https://doi.org/10.1016/j.engstruct.2023.117106"),
)


def _roots() -> list[Path]:
    here = Path(__file__).resolve().parent
    roots = [here / "docs", here.parent / "docs" / "pdf", here.parent / "docs"]
    base = getattr(sys, "_MEIPASS", None)            # PyInstaller bundle
    if base:
        roots += [Path(base) / "rvestudio" / "docs", Path(base) / "docs" / "pdf"]
    return roots


def guide_path(kind: str) -> Path | None:
    """PDF of the 'user' or 'theory' guide, or its Markdown source when no PDF is present; None if neither exists."""
    pdf, md, _ = GUIDES[kind]
    for name in (pdf, md):
        for root in _roots():
            p = root / name
            if p.is_file():
                return p
    return None
