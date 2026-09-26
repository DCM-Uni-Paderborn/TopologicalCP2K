"""Numerical gap bounds on finite-localizer parameter boxes.

Parameters are (E, kappa, x, y), in hartree and bohr. The Hamiltonian,
orthonormalized positions and reference origin remain fixed. This is a
Weyl/norm estimate with a numerical margin, not interval arithmetic.
Thomas D. Kuehne, tkuehne@cp2k.org
"""

import numpy as np


def variation_bound(lower, upper, position_norm, reference_xy):
    """Bound ||L(p)-L(midpoint)|| throughout an axis-aligned parameter box.

    Delta L has diagonal blocks -dE I and dE I. Its off-diagonal
    block is dk*((X-xc*I)-i*(Y-yc*I)) - k*(dx-i*dy)*I. Squaring Delta L
    therefore gives ||Delta L||=hypot(dE,||offdiag||).
    """
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    origin = np.asarray(reference_xy, float)
    if (lower.shape != (4,) or upper.shape != (4,) or origin.shape != (2,)
            or not np.isfinite([*lower, *upper, *origin, position_norm]).all()
            or np.any(lower > upper) or lower[1] <= 0 or position_norm < 0):
        raise ValueError("Invalid localizer box or norm")
    center = lower + (upper - lower) / 2
    radius = (upper - lower) / 2
    norm_at_center = position_norm + np.linalg.norm(center[2:] - origin)
    position_variation = upper[1] * np.linalg.norm(radius[2:])
    scale_variation = radius[1] * norm_at_center
    bound = float(np.hypot(radius[0], scale_variation + position_variation))
    contributions = np.r_[radius[0], scale_variation, upper[1] * radius[2:]]
    if not np.isfinite([bound, *contributions]).all():
        raise ValueError("Unresolved finite parameter variation")
    return bound, contributions


def cover_region(lower, upper, position_norm, reference_xy, evaluate,
                 margin=1e-10, max_nodes=1023):
    """Subdivide until every leaf has a positive gap bound, or retain failures.

    The box is connected; a valid continuous class-AII localizer then has
    one Pfaffian index throughout a resolved cover, determined separately
    at an anchor. This function itself never assigns an index.
    """
    variation_bound(lower, upper, position_norm, reference_xy)
    if not np.isfinite(margin) or margin <= 0 or type(max_nodes) is not int or max_nodes < 1:
        raise ValueError("Invalid margin or evaluation budget")
    pending = [(np.asarray(lower, float), np.asarray(upper, float))]
    covered, unresolved, evaluations = [], [], []
    while pending:
        left, right = pending.pop()
        if len(evaluations) >= max_nodes:
            unresolved.append(dict(lower=left.tolist(), upper=right.tolist(), reason="budget"))
            continue
        middle = left + (right - left) / 2
        variation, effects = variation_bound(left, right, position_norm, reference_xy)
        gap = float(evaluate(middle))
        if not np.isfinite(gap) or gap < 0:
            raise ValueError("Invalid localizer gap")
        row = dict(lower=left.tolist(), upper=right.tolist(), midpoint=middle.tolist(),
                   gap=gap, variation_bound=variation, lower_bound=gap - variation - margin)
        evaluations.append(row)
        if row["lower_bound"] > 0:
            covered.append(row)
        else:
            split = int(np.argmax(effects))
            if middle[split] in (left[split], right[split]) or effects[split] == 0:
                unresolved.append(dict(lower=left.tolist(), upper=right.tolist(), reason="resolution"))
                continue
            lower_right, upper_left = left.copy(), right.copy()
            lower_right[split] = upper_left[split] = middle[split]
            pending.extend([(lower_right, right), (left, upper_left)])
    return dict(covered=covered, unresolved=unresolved, evaluations=evaluations,
                resolved=not unresolved, roundoff_margin_hartree=margin)
