"""Reliability of the one-dimensional wall model for a given design (theory guide, Section 9).

The analytical model treats every wall as a beam in cylindrical bending (twist suppressed) placed on its neutral
axis. Three measures say how far that idealisation can be trusted for the effective properties:

* the bend-twist index  Pi = eta psi^2, with the bend-twist ratio psi = D16 / sqrt(D11 D66) of the wall laminate and
  the twist release eta = (2.0 + 2.15 psi^2) lambda^0.73, lambda = sqrt(D66/A44)/L the twist length of the wall.
  Pi is the fraction by which the cylindrical wall over-estimates the wall bending stiffness; the bending-governed
  in-plane properties (E1, E2, G12, nu13) are over-predicted by about 100 Pi %;
* the membrane-bending coupling  chi_B = B11^2 / (A11 D11)  of an unsymmetric wall. The stretch-governed terms
  (C1111, C2222, C1122, Q11, Q12) are over-predicted by about chi_B;
* the wall shear parameter  beta = 12 D_b / (L^2 S_eq): with the Euler-Bernoulli wall the bending-governed terms
  are over-predicted by about 100*beta %, in addition to the bend-twist error, so the Euler-Bernoulli cylindrical
  wall is graded on Pi + beta; with the shear-flexible wall the shear model itself is accurate to about
  2 % for beta <= 0.10 and to about 5 % up to beta = 0.5.

The out-of-plane properties (E3, G13, G23, nu31) depend on none of the bending measures and are reliable for every
symmetric wall.

The 'twist-corrected' bending option replaces D11 by D11 (1 - eta psi^2) with the calibrated twist release
eta = (2.0 + 2.15 psi^2) lambda^0.73; with it the in-plane properties of every tested wall with beta <= 0.10
(|psi| up to 0.85, 0.010 <= lambda <= 0.100) are within 2 % of the converged shell FE (limits below).
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np

from .laminate import TWIST_LENGTH_RANGE, TWIST_RELEASE_DEFAULT, twist_length, twist_release

# Cylindrical wall (original model), shear-flexible beam: limits on the bend-twist index Pi = eta psi^2.
PI_GOOD = 0.018                   # in-plane error within about 2 %
PI_FAIR = 0.045                   # in-plane error within about 5 %
# Twist-corrected wall: limits on |psi| and on the calibrated twist-length range.
PSI_TC_GOOD = 0.85                # within 2 % for beta <= 0.10 over the whole tested range
PSI_TC_FAIR = 0.85                # beyond the tested range: outside
LAMBDA_RANGE = TWIST_LENGTH_RANGE  # twist lengths covered by the calibration
# Shear-flexible wall: accuracy of the wall shear model itself.
BETA_TIM_GOOD = 0.10              # shear-flexible wall: about 2 % for beta <= 0.10, about 5 % up to 0.5
BETA_TIM_FAIR = 0.50
# Unsymmetric walls: stretch-governed terms.
CHI_B_GOOD = 0.02
CHI_B_FAIR = 0.05
# Euler-Bernoulli walls: wall shear parameter.
BETA_GOOD = 0.02
BETA_FAIR = 0.05

LEVELS = ("reliable", "use with care", "outside the verified range")


@dataclass
class WallReliability:
    laminate: int
    psi: float
    twist_length: float
    twist_release: float
    bend_twist_index: float
    chi_B: float
    beta_max: float
    level: str
    notes: list = field(default_factory=list)


@dataclass
class Reliability:
    level: str
    in_plane: str
    out_of_plane: str
    option: str
    beam: str
    walls: list
    messages: list

    def as_dict(self) -> dict:
        return {"level": self.level, "in_plane": self.in_plane, "out_of_plane": self.out_of_plane,
                "option": self.option, "beam": self.beam,
                "walls": [asdict(w) for w in self.walls], "messages": list(self.messages),
                "limits": {"pi_good": PI_GOOD, "pi_fair": PI_FAIR,
                           "psi_tc_good": PSI_TC_GOOD, "psi_tc_fair": PSI_TC_FAIR, "lambda_range": list(LAMBDA_RANGE),
                           "chi_B_good": CHI_B_GOOD, "chi_B_fair": CHI_B_FAIR,
                           "beta_good": BETA_GOOD, "beta_fair": BETA_FAIR, "beta_shear_flexible_good": BETA_TIM_GOOD,
                           "beta_shear_flexible_fair": BETA_TIM_FAIR}}

    def short(self) -> str:
        return f"in-plane: {self.in_plane}; out-of-plane: {self.out_of_plane}"


def _rank(*levels: str) -> str:
    return LEVELS[max(LEVELS.index(x) for x in levels)] if levels else LEVELS[0]


def _grade(value: float, good: float, fair: float) -> str:
    v = abs(value)
    return LEVELS[0] if v <= good else (LEVELS[1] if v <= fair else LEVELS[2])


def section_measures(sec) -> dict:
    """psi and chi_B of a LaminateSection (D about the mid-plane)."""
    D = np.asarray(getattr(sec, "Doriginal", sec.D), float)
    A = np.asarray(sec.A, float)
    B = np.asarray(sec.B, float)
    psi = float(D[0, 2] / math.sqrt(D[0, 0] * D[2, 2])) if D[0, 0] > 0 and D[2, 2] > 0 else 0.0
    chi = float(B[0, 0] ** 2 / (A[0, 0] * D[0, 0])) if A[0, 0] > 0 and D[0, 0] > 0 else 0.0
    return {"psi": psi, "chi_B": chi}


def _infer(result):
    meta = result.element_data[0].get("strip_meta", {}) if result.element_data else {}
    beam = meta.get("beam", "euler-bernoulli")
    mode = str(meta.get("mode", "paper")).lower().replace("_", "")
    if mode in ("paper", "cylindrical", "thesis"):
        option = "cylindrical"
    elif mode.startswith("twist"):
        option = "twist-corrected"
    else:
        option = "free moment" if mode.startswith("free") else "anisotropic strip"
    return beam, option


def _user_coeffs(result):
    """Twist-release coefficients of the solve when they differ from the calibrated defaults, else None."""
    meta = result.element_data[0].get("strip_meta", {}) if result.element_data else {}
    c = meta.get("twist_coeffs")
    if c is None or all(abs(float(x) - float(y)) <= 1e-12 for x, y in zip(c, TWIST_RELEASE_DEFAULT)):
        return None
    return tuple(float(x) for x in c)


def assess(result, wall_shear: str | None = None, wall_option: str | None = None) -> Reliability:
    """Reliability of a solved design (HomogenizationResult) from its wall laminates and wall lengths.

    wall_shear / wall_option default to the ones the result was solved with."""
    from .walls import is_timoshenko
    beam, option = _infer(result)
    wall_shear = wall_shear or beam
    wall_option = wall_option or option
    shear = is_timoshenko(wall_shear)
    lengths: dict[int, float] = {}
    for e in result.element_data:
        k = int(e["laminate"])
        lengths[k] = min(lengths.get(k, math.inf), float(e["L"]))
    walls, msgs = [], []
    bound = wall_option in ("free moment", "anisotropic strip")
    if bound:
        msgs.append(f"the '{wall_option}' option is a bound on the wall bending stiffness, not the verified model; "
                    "the limits apply to the cylindrical and twist-corrected options")
    for k, Lmin in sorted(lengths.items()):
        sec = result.sections[k]
        m = section_measures(sec)
        lam = twist_length(sec, Lmin) if math.isfinite(Lmin) else 0.0
        eta = twist_release(m["psi"], lam)
        pi_bt = eta * m["psi"] ** 2
        Db = float(np.asarray(sec.D, float)[0, 0])
        beta = 12 * Db / (Lmin * Lmin * sec.Seq) if sec.Seq > 0 and math.isfinite(Lmin) else 0.0
        notes = []
        if wall_option == "twist-corrected":
            lv = [_grade(m["psi"], PSI_TC_GOOD, PSI_TC_FAIR)]
            user = _user_coeffs(result)
            if user is not None and abs(m["psi"]) > 0.05:
                lv.append(LEVELS[1])
                notes.append("user-defined twist-release coefficients (a, b, n) = "
                             f"({user[0]:g}, {user[1]:g}, {user[2]:g}): the limits were verified for the calibrated "
                             f"values {TWIST_RELEASE_DEFAULT}")
            if abs(m["psi"]) > PSI_TC_GOOD:
                notes.append(f"psi = {m['psi']:.3f} is beyond the calibrated range of the twist-corrected wall (|psi| <= 0.85); "
                             "verify with the FE check")
            if abs(m["psi"]) > 0.05 and not (LAMBDA_RANGE[0] <= lam <= LAMBDA_RANGE[1]):
                lv.append(LEVELS[1])
                notes.append(f"twist length lambda = {lam:.4f} is outside the calibrated range "
                             f"{LAMBDA_RANGE[0]:.3f}-{LAMBDA_RANGE[1]:.3f}")
        else:
            lv = [_grade(pi_bt, PI_GOOD, PI_FAIR)]
            if pi_bt > PI_GOOD:
                notes.append(f"bend-twist index Pi = {pi_bt:.4f} (psi = {m['psi']:.3f}, lambda = {lam:.4f}): "
                             f"E1, E2, G12 and nu13 are over-predicted by about {100 * pi_bt:.0f} %"
                             + "; the twist-corrected option extends the range")
        if shear and beta > BETA_TIM_GOOD:
            lv.append(LEVELS[1] if beta <= BETA_TIM_FAIR else LEVELS[2])
            notes.append(f"thick wall, beta = {beta:.3f}: the shear-flexible wall model is accurate to about 5 % "
                         "rather than 2 % in the bending-governed terms" + ("" if beta <= BETA_TIM_FAIR else
                                                                        "; beta is beyond the tested range (0.5)"))
        lv.append(_grade(m["chi_B"], CHI_B_GOOD, CHI_B_FAIR))
        if m["chi_B"] > CHI_B_GOOD:
            notes.append(f"unsymmetric wall, chi_B = {m['chi_B']:.3f}: C1111, C2222, C1122 (Q11, Q12) "
                         f"are over-predicted by about {100 * m['chi_B']:.0f} %")
        if not shear:
            # The wall-shear error (about 100 beta %) and the bend-twist error of the cylindrical wall (about 100 Pi %)
            # act on the same bending-governed terms and add: grade the Euler-Bernoulli cylindrical wall on Pi + beta
            # (320-wall study: largest in-plane error 2.0 % for Pi + beta <= 0.02 and 5.2 % for Pi + beta <= 0.05).
            idx = beta + (0.0 if wall_option == "twist-corrected" else pi_bt)
            lv.append(_grade(idx, BETA_GOOD, BETA_FAIR))
            if idx > BETA_GOOD:
                what = (f"beta = {beta:.3f}" if wall_option == "twist-corrected" else
                        f"beta + Pi = {beta:.3f} + {pi_bt:.3f}")
                notes.append(f"Euler-Bernoulli wall with {what}: bending-governed terms are about "
                             f"{100 * idx:.0f} % too stiff; use the shear-flexible (Timoshenko) wall")
        walls.append(WallReliability(k, m["psi"], lam, eta, pi_bt, m["chi_B"], beta, _rank(*lv), notes))
    in_plane = LEVELS[2] if bound else _rank(*[w.level for w in walls])
    out_of_plane = _rank(*[_grade(w.chi_B, CHI_B_GOOD, CHI_B_FAIR) for w in walls])
    for w in walls:
        for n in w.notes:
            msgs.append(f"laminate {w.laminate}: {n}")
    return Reliability(_rank(in_plane, out_of_plane), in_plane, out_of_plane, wall_option,
                       "timoshenko" if shear else "euler-bernoulli", walls, msgs)


def assess_model(model, result) -> Reliability:
    """Same as assess(result); kept for callers that hold the model."""
    return assess(result)
