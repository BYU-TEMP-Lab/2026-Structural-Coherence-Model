"""
batch_builder_from_table.py
Read batch_table.csv and generate a Python file batch_plot_configs.py containing
CONFIGS = [ ... ] dicts suitable for direct copy/paste into batch_plot.py.

This script DOES NOT execute any plotting. It only writes the configs file and
prints where it wrote it.
"""

import pandas as pd
from typing import List
import json


def parse_methods(row) -> List[str]:
    methods: List[str] = []
    if str(row.get('Gheribi-KT24 at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Gheribi-KT24')
    if str(row.get('Gheribi-KT24, Mix Data at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Gheribi-KT24, Mix Data')
    if str(row.get('Zhao-PGM at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Zhao-PGM')
    if str(row.get('Zhao-PGM, Mix Data at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Zhao-PGM, Mix Data')
    if str(row.get('Present Model at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Present Model')
    if str(row.get('Present Model, Mix Data at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Present Model, Mix Data')
    if str(row.get('Ideal at Melt Temp (W/m-K)', 0)) == '1':
        methods.append('Ideal')
    return methods


def parse_sources(s) -> List[str]:
    """Split references on commas that are OUTSIDE parentheses, preserving commas inside '(...)'."""
    if not isinstance(s, str) or not s.strip():
        return []
    parts: List[str] = []
    buf: List[str] = []
    depth = 0
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
        if ch == ',' and depth == 0:
            token = ''.join(buf).strip()
            if token:
                parts.append(token)
            buf = []
        else:
            buf.append(ch)
    # last token
    token = ''.join(buf).strip()
    if token:
        parts.append(token)
    return parts


def to_python_literal(v):
    """Render Python code literal for scalars/lists/tuples/strings without json's null/true/false."""
    if isinstance(v, str):
        # Escape backslashes and quotes for safety
        return repr(v)
    if isinstance(v, tuple):
        return f"({to_python_literal(v[0])}, {to_python_literal(v[1])})"
    if isinstance(v, list):
        return '[' + ', '.join(to_python_literal(x) for x in v) + ']'
    if isinstance(v, (int, float)):
        return repr(v)
    if v is None:
        return 'None'
    # Fallback via json then patch tokens
    s = json.dumps(v)
    return s.replace('null', 'None').replace('true', 'True').replace('false', 'False')


def build_configs_from_csv(csv_path: str) -> List[dict]:
    df = pd.read_csv(csv_path)
    configs: List[dict] = []
    for _, r in df.iterrows():
        comp = str(r.get('Salt Composition', '')).strip()
        src = str(r.get('PDF Source', '')).strip()
        # Skip blank or NaN-like rows
        if (not comp or comp.lower() == 'nan') or (not src or src.lower() == 'nan'):
            # Skip blank/separator rows
            continue
        try:
            melt = float(r.get('Melting Temp (K)', ''))
        except Exception:
            continue
        methods = parse_methods(r)
        cfg = {
            'composition': comp,
            'temp_range': (melt, 1500.0),
            'methods': methods,
            'measurement_sources': parse_sources(r.get('Exp. Reference', '')),
            'scl_composition_with_source': f"{comp} ({src})",
            'save_results_csv': True,
            'show_plot': False,
        }
        configs.append(cfg)
    return configs


def write_configs_py(configs: List[dict], out_path: str) -> None:
    lines: List[str] = []
    lines.append("CONFIGS = [")
    for cfg in configs:
        lines.append("    {")
        for k in ['composition', 'temp_range', 'methods', 'measurement_sources', 'scl_composition_with_source', 'save_results_csv', 'show_plot']:
            v = cfg.get(k)
            if k in ('save_results_csv', 'show_plot'):
                lit = 'True' if v else 'False'
            else:
                lit = to_python_literal(v)
            lines.append(f"        '{k}': {lit},")
        lines.append("    },")
    lines.append("]\n")
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


if __name__ == "__main__":
    configs = build_configs_from_csv('batch_table.csv')
    out_file = 'batch_plot_configs.py'
    write_configs_py(configs, out_file)
    print(f"Wrote {len(configs)} configs to {out_file}. Import CONFIGS from that file or copy/paste into batch_plot.py.")