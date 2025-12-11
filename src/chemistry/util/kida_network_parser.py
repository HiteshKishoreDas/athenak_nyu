# -*- coding: utf-8 -*-
"""Minimal functional helpers for parsing KIDA species and reaction files.

This replaces the unused KidaNetwork class with two simple functions that mirror
the current functional loader style.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import math


def parse_species_file(path: Path, element_order: List[str]) -> Dict[str, Dict]:
    """Parse a KIDA species file into a dict keyed by species name."""
    species = {}
    try:
        lines = path.read_text().splitlines()
    except FileNotFoundError:
        return species

    for line in lines:
        line = line.strip()
        if not line or line.startswith("!"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name = parts[0]
        try:
            charge = int(parts[1])
        except Exception:
            charge = 0
        composition = {}
        counts = parts[2:]
        for idx, count_str in enumerate(counts):
            try:
                val = int(count_str)
                if val > 0 and idx < len(element_order):
                    composition[element_order[idx]] = val
            except Exception:
                continue
        elements = {e: 0 for e in element_order}
        for e, v in composition.items():
            elements[e] = v
        species[name] = {"charge": charge, "composition": composition, "elements": elements}
    return species


def parse_reactions_file(path: Path) -> List[Dict]:
    """Parse a KIDA reactions file into a list of reaction dicts."""
    reactions: List[Dict] = []
    try:
        lines = path.read_text().splitlines()
    except FileNotFoundError:
        return reactions

    for line in lines:
        if line.startswith("!") or len(line) < 80:
            continue
        try:
            r1 = line[0:11].strip()
            r2 = line[11:22].strip()
            r3 = line[22:33].strip()
            reactants = [r for r in (r1, r2, r3) if r]
            products = []
            for i in range(5):
                p = line[34 + (i * 11) : 34 + ((i + 1) * 11)].strip()
                if p:
                    products.append(p)
            alpha = float(line[90:100].strip() or 0.0)
            beta = float(line[101:111].strip() or 0.0)
            gamma = float(line[112:122].strip() or 0.0)
            f_val = float(line[123:131].strip() or 0.0)
            g_val = float(line[132:140].strip() or 0.0)
            remainder = line[140:].split()
            try:
                itype = int(remainder[0]) if len(remainder) > 0 else 0
            except Exception:
                itype = 0
            try:
                t_min = int(remainder[1]) if len(remainder) > 1 else -9999
                t_max = int(remainder[2]) if len(remainder) > 2 else 9999
            except Exception:
                t_min, t_max = -9999, 9999
            reactions.append(
                {
                    "reactants": reactants,
                    "products": products,
                    "alpha": alpha,
                    "beta": beta,
                    "gamma": gamma,
                    "f": f_val,
                    "g": g_val,
                    "itype": itype,
                    "t_min": t_min,
                    "t_max": t_max,
                }
            )
        except Exception:
            continue
    return reactions
