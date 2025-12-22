# -*- coding: utf-8 -*-
"""Minimal KIDA network loader (functional, no classes)."""

from __future__ import annotations

import math
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set


def _default_elements() -> List[str]:
    return [
        "H",
        "He",
        "C",
        "N",
        "O",
        "Si",
        "S",
        "Fe",
        "Na",
        "Mg",
        "Cl",
        "P",
        "F",
        "D",
        "Al",
        "Li",
        "Be",
        "B",
        "Ca",
        "Br",
        "Ar",
        "Ne",
    ]


def _read_lines(path: Path) -> List[str]:
    try:
        return path.read_text().splitlines()
    except FileNotFoundError:
        return []


def _parse_metadata(path: Path) -> Dict:
    meta = {}
    lines = _read_lines(path)
    for line in lines:
        if "Number of reactions" in line:
            try:
                meta["num_reactions"] = int(line.split(":")[1].strip())
            except Exception:
                pass
        if "Number of species" in line:
            try:
                meta["num_species"] = int(line.split(":")[1].strip())
            except Exception:
                pass
        if "Format of the file" in line:
            meta["format_desc"] = line.split(":")[1].strip()
        if "Elements are the following:" in line:
            elems = line.split("following:")[1].strip()
            meta["elements"] = [e.strip() for e in elems.split(",")]
    return meta


def _load_species_file(path: Path, element_order: List[str]) -> Dict[str, Dict]:
    species = {}
    for line in _read_lines(path):
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
        species[name] = {
            "charge": charge,
            "composition": composition,
            "elements": elements,
        }
    return species


def _load_reactions_file(path: Path) -> List[Dict]:
    rxns = []
    for line in _read_lines(path):
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
            rxns.append(
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
    return rxns


def load_network(directory: Path) -> Dict:
    directory = Path(directory)
    meta = {}
    element_order: List[str] = []
    for f in sorted(directory.glob("kida_spec_readme_*.txt")):
        m = _parse_metadata(f)
        meta.update(m)
        if "elements" in m and not element_order:
            element_order = m["elements"]
    for f in sorted(directory.glob("kida_readme_*.txt")):
        m = _parse_metadata(f)
        meta.update(m)
        if "elements" in m and not element_order:
            element_order = m["elements"]
    if not element_order:
        element_order = _default_elements()
    species = {}
    for f in sorted(directory.glob("kida_spec_*.dat")):
        species.update(_load_species_file(f, element_order))
    reactions: List[Dict] = []
    for f in sorted(directory.glob("kida_reac_*.dat")):
        reactions.extend(_load_reactions_file(f))
    surface = {}
    for f in sorted(directory.glob("kida_surface_*.dat")):
        surface[f.name] = [
            ln for ln in _read_lines(f) if ln.strip() and not ln.startswith("!")
        ]
    return {
        "metadata": meta,
        "element_order": element_order,
        "species": species,
        "reactions": reactions,
        "surface": surface,
    }


def visualize_network(
    network: Dict,
    out_path: Optional[Path | str] = None,
    max_reactions: int = 50,
    temperature: float = 100.0,
) -> str:
    """Build a GraphViz DOT string for the reaction network.

    The graph is bipartite: species (ellipses) connect to reaction nodes
    (rounded boxes) which then connect to product species. External drivers
    like CR or Photon are omitted to keep the view compact. If `out_path`
    is provided, the DOT text is written to that file.

    Parameters:
        network: loaded network dict from `load_network`
        out_path: optional path to save DOT text
        max_reactions: cap the number of reactions to include (to avoid giant graphs)
        temperature: temperature used to annotate k values (Arrhenius coefficient)
    Returns:
        DOT-format string describing the graph
    """
    reactions = network.get("reactions") or []
    species = network.get("species") or {}
    if not reactions or not species:
        raise ValueError("Network is empty; load a network before visualizing.")

    def esc(name: str) -> str:
        # basic quote escaping for DOT labels
        return name.replace("\\", "\\\\").replace('"', '\\"')

    lines: List[str] = [
        "digraph ChemicalNetwork {",
        "  rankdir=LR;",
        '  node [shape=ellipse, style=filled, fillcolor="#eef4ff", color="#33415c"];',
        '  edge [color="#5a6e8c"];',
    ]

    used_species: Set[str] = set()

    for idx, rxn in enumerate(reactions[:max_reactions]):
        reactants = [r for r in rxn.get("reactants", []) if r not in {"CR", "CRP", "Photon"}]
        products = [p for p in rxn.get("products", []) if p not in {"CR", "CRP", "Photon"}]
        rid = f"R{idx}"
        rate = calculate_rate_coefficient(rxn, temperature=temperature)
        left = " + ".join(reactants) if reactants else "(none)"
        right = " + ".join(products) if products else "(none)"
        eqn = f"{left} -> {right}"
        label = f"{rid}: {eqn}\\n k={rate:.2e}"
        lines.append(
            f'  "{rid}" [shape=box, style="rounded,filled", fillcolor="#dde7ff", label="{esc(label)}"];'
        )

        for r in reactants:
            if r in species:
                used_species.add(r)
                lines.append(f'  "{esc(r)}" -> "{rid}";')
        for p in products:
            if p in species:
                used_species.add(p)
                lines.append(f'  "{rid}" -> "{esc(p)}";')

    for s in sorted(used_species):
        charge = species[s].get("charge", 0)
        label = f"{s} (q={charge})" if charge else s
        lines.append(f'  "{esc(s)}" [label="{esc(label)}"];')

    lines.append("}")
    dot = "\n".join(lines)

    if out_path:
        out_path = Path(out_path)
        out_path.write_text(dot)

    return dot


def calculate_rate_coefficient(
    rxn: Dict, temperature: float = 300.0, av: float = 10.0, cr_zeta: float = 1.3e-17
) -> float:
    """Simple Arrhenius-style rate coefficient."""
    alpha = rxn.get("alpha", 0.0)
    beta = rxn.get("beta", 0.0)
    gamma = rxn.get("gamma", 0.0)
    if temperature <= 0.0:
        temperature = 1.0
    try:
        return alpha * (temperature / 300.0) ** beta * math.exp(-gamma / temperature)
    except OverflowError:
        return 0.0


if __name__ == "__main__":
    from pprint import pprint

    # Default to the repository custom directory next to this file
    base_dir = Path(__file__).resolve().parent / "custom"
    net = load_network(base_dir)
    if not net.get("species") or not net.get("reactions"):
        print(f"Failed to load network from {base_dir}")
        raise SystemExit(1)

    print(f"Loaded network from {base_dir}")
    print(f"  Species: {len(net['species'])}")
    print(f"  Reactions: {len(net['reactions'])}")
    print(f"  Elements: {net.get('element_order', [])}")
    print("\nSample species (up to 5):")
    for name in list(net["species"].keys())[:5]:
        pprint({name: net["species"][name]})
    print("\nSample reactions (up to 5):")
    for rxn in net["reactions"][:5]:
        pprint(rxn)

    # Generate a DOT visualization and render to PNG under test_output/
    try:
        # Prefer existing pluralized directory; fall back to singular if missing
        base_dir = Path(__file__).resolve().parent
        candidate = base_dir / "test_outputs"
        out_dir = candidate if candidate.exists() else base_dir / "test_output"
        out_dir.mkdir(parents=True, exist_ok=True)

        dot_path = out_dir / "network.dot"
        png_path = out_dir / "network.png"

        visualize_network(net, out_path=dot_path, max_reactions=25, temperature=100.0)
        print(f"\nVisualization DOT written to {dot_path}")

        dot_cmd = shutil.which("dot")
        if dot_cmd:
            subprocess.run([dot_cmd, "-Tpng", str(dot_path), "-o", str(png_path)], check=True)
            print(f"PNG rendered to {png_path}")
        else:
            print("GraphViz `dot` not found; install graphviz to render the PNG.")
            print(f"You can run: dot -Tpng {dot_path} -o {png_path}")
    except Exception as exc:
        print(f"\nFailed to generate visualization DOT: {exc}")
