from typing import List, Dict, Tuple


def _methods_from_flags(
    g_kt24: int,
    g_kt24_mix: int,
    z_pgm: int,
    z_pgm_mix: int,
    present: int,
    present_mix: int,
    ideal: int,
) -> List[str]:
    methods: List[str] = []
    if g_kt24: methods.append("Gheribi-KT24")
    if g_kt24_mix: methods.append("Gheribi-KT24, Mix Data")
    if z_pgm: methods.append("Zhao-PGM")
    if z_pgm_mix: methods.append("Zhao-PGM, Mix Data")
    if present: methods.append("Present Model")
    if present_mix: methods.append("Present Model, Mix Data")
    if ideal: methods.append("Ideal")
    return methods


def _sources_list(exp_refs: str) -> List[str]:
    if not exp_refs:
        return []
    return [s.strip() for s in exp_refs.split(',') if s.strip()]


def _cfg(comp: str, pdf_source: str, melt_k: float,
         g_kt24: int, g_kt24_mix: int, z_pgm: int, z_pgm_mix: int, present: int, present_mix: int, ideal: int,
         exp_refs: str) -> Dict:
    return dict(
        composition=comp,
        temp_range=(float(melt_k), 1500.0),
        methods=_methods_from_flags(g_kt24, g_kt24_mix, z_pgm, z_pgm_mix, present, present_mix, ideal),
        measurement_sources=_sources_list(exp_refs),
        scl_composition_with_source=f"{comp} ({pdf_source})",
        save_results_csv=True,
        show_plot=False,
    )


def get_configs() -> List[Dict]:
    cfgs: List[Dict] = []

    # 1.0LiF
    cfgs.append(_cfg(
        "1.0LiF", "Walz, 2019", 1121.2,
        1, 1, 1, 1, 1, 1, 1,
        "LiF (Golyshev, 1992),LiF (Smirnov, 1987),LiF (Khlebnikov, 1981)"
    ))

    # 1.0NaF
    cfgs.append(_cfg(
        "1.0NaF", "Walz, 2019", 1268,
        1, 1, 1, 1, 1, 1, 1,
        "NaF (Smirnov, 1987),NaF (Polyakov, 1975)"
    ))

    # 1.0KF
    cfgs.append(_cfg(
        "1.0KF", "Walz, 2019", 1131.2,
        1, 1, 1, 1, 1, 1, 1,
        "KF (Smirnov, 1987)"
    ))

    # 1.0LiCl
    cfgs.append(_cfg(
        "1.0LiCl", "Walz, 2019", 883,
        1, 1, 1, 1, 1, 1, 1,
        "LiCl (Nagasaka, 1992)"
    ))

    # 1.0NaCl (Walz, 2019)
    cfgs.append(_cfg(
        "1.0NaCl", "Walz, 2019", 1073.8,
        1, 1, 1, 1, 1, 1, 1,
        "NaCl (Harada, 1992),NaCl (Nagasaka 1992)"
    ))

    # 1.0NaCl (Lu, 2021)
    cfgs.append(_cfg(
        "1.0NaCl", "Lu, 2021", 1073.8,
        1, 1, 1, 1, 1, 1, 1,
        "NaCl (Harada, 1992),NaCl (Nagasaka 1992)"
    ))

    # 1.0KCl
    cfgs.append(_cfg(
        "1.0KCl", "Walz, 2019", 1042.7,
        1, 1, 1, 1, 1, 1, 1,
        "KCl (Harada, 1992),KCl (Nagasaka, 1992)"
    ))

    # 1.0MgCl2 (Roy, 2021)
    cfgs.append(_cfg(
        "1.0MgCl2", "Roy, 2021", 987,
        1, 1, 1, 1, 1, 1, 1,
        "MgCl2 (Filatov, 2005)"
    ))

    # 1.0MgCl2 (McGreevy, 1987)
    cfgs.append(_cfg(
        "1.0MgCl2", "McGreevy, 1987", 987,
        1, 1, 1, 1, 1, 1, 1,
        "MgCl2 (Filatov, 2005)"
    ))

    # 1.0MgCl2 (Lu, 2021)
    cfgs.append(_cfg(
        "1.0MgCl2", "Lu, 2021", 987,
        1, 1, 1, 1, 1, 1, 1,
        "MgCl2 (Filatov, 2005)"
    ))

    # 1.0CaCl2 (McGreevy, 1987)
    cfgs.append(_cfg(
        "1.0CaCl2", "McGreevy, 1987", 1045,
        1, 1, 1, 1, 1, 1, 1,
        "CaCl2 (Wei, 2022; MD)"
    ))

    # 1.0CaCl2 (Bu, 2021)
    cfgs.append(_cfg(
        "1.0CaCl2", "Bu, 2021", 1045,
        1, 1, 1, 1, 1, 1, 1,
        "CaCl2 (Wei, 2022; MD)"
    ))

    # 1.0SrCl2 (no experimental)
    cfgs.append(_cfg(
        "1.0SrCl2", "McGreevy, 1987", 1146,
        1, 1, 1, 1, 1, 0, 0,
        ""
    ))

    # 0.6LiF-0.4NaF
    cfgs.append(_cfg(
        "0.6LiF-0.4NaF", "Grizzi, 2024", 927.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.6LiF-0.4NaF (BYU)"
    ))

    # 0.5BeF2-0.5LiF
    cfgs.append(_cfg(
        "0.5BeF2-0.5LiF", "Sun, 2024", 648.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.66LiF-0.34BeF2 (Kato, 1983)"
    ))

    # 0.34BeF2-0.66LiF (Fayfar, 2023) entry 1
    cfgs.append(_cfg(
        "0.34BeF2-0.66LiF", "Fayfar, 2023", 713.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.66LiF-0.34BeF2 (Redkin, 2022),0.66LiF-0.34BeF2 (Kato, 1983),0.66LiF-0.34BeF2 (Bobrova, 2023)"
    ))

    # 0.34BeF2-0.66LiF (Fayfar, 2023) entry 2
    cfgs.append(_cfg(
        "0.34BeF2-0.66LiF", "Fayfar, 2023", 713.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.66LiF-0.34BeF2 (Schorne-Pinto, 2024)"
    ))

    # 0.42KF-0.465LiF-0.115NaF
    cfgs.append(_cfg(
        "0.42KF-0.465LiF-0.115NaF", "Frandsen, 2020", 735,
        1, 1, 1, 1, 1, 1, 1,
        "0.465LiF-0.115NaF-0.42KF (Merritt, 2022),0.465LiF-0.115NaF-0.42KF (Gallagher, 2022),0.465LiF-0.115NaF-0.42KF (Rudenko, 2022)"
    ))

    # 0.59KF-0.065MgF2-0.345NaF
    cfgs.append(_cfg(
        "0.59KF-0.065MgF2-0.345NaF", "Solano, 2021", 958,
        1, 1, 1, 1, 1, 1, 1,
        "0.345NaF-0.59KF-0.065MgF2 (Rudenko, 2024)"
    ))

    # separator row in table (ignored)

    # 0.5KCl-0.5NaCl
    cfgs.append(_cfg(
        "0.5KCl-0.5NaCl", "Manga, 2013", 927.15,
        1, 1, 1, 1, 1, 1, 1,
        ""
    ))

    # 0.3CaCl2-0.7LiCl
    cfgs.append(_cfg(
        "0.3CaCl2-0.7LiCl", "Liang, 2024", 771.15,
        1, 1, 1, 1, 1, 1, 1,
        ""
    ))

    # 0.5097CaCl2-0.4903NaCl (entry 1)
    cfgs.append(_cfg(
        "0.5097CaCl2-0.4903NaCl", "Wei, 2022", 779.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.52NaCl-0.48CaCl2 (Tian, 2021)"
    ))

    # 0.5097CaCl2-0.4903NaCl (entry 2)
    cfgs.append(_cfg(
        "0.5097CaCl2-0.4903NaCl", "Wei, 2022", 779.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.4903NaCl-0.5097CaCl2 (Wei, 2022; MD),0.52NaCl-0.48CaCl2 (Tian, 2021)"
    ))

    # 0.282CaCl2-0.718KCl (two entries, both with no exp refs)
    cfgs.append(_cfg(
        "0.282CaCl2-0.718KCl", "Wei, 2022", 921.15,
        1, 1, 1, 1, 1, 1, 1,
        ""
    ))
    cfgs.append(_cfg(
        "0.282CaCl2-0.718KCl", "Wei, 2022", 921.15,
        1, 1, 1, 1, 1, 1, 1,
        ""
    ))

    # 0.22KCl-0.45MgCl2-0.33NaCl
    cfgs.append(_cfg(
        "0.22KCl-0.45MgCl2-0.33NaCl", "Jiang, 2024", 658.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.1511NaCl-0.3891KCl-0.4598MgCl2 (Wang, 2021)"
    ))

    # 0.41KCl-0.38MgCl2-0.21NaCl
    cfgs.append(_cfg(
        "0.41KCl-0.38MgCl2-0.21NaCl", "Jiang, 2024", 660.15,
        1, 1, 1, 1, 1, 1, 1,
        "0.1511NaCl-0.3891KCl-0.4598MgCl2 (Wang, 2021)"
    ))

    # 0.525CaCl2-0.058KCl-0.417NaCl
    cfgs.append(_cfg(
        "0.525CaCl2-0.058KCl-0.417NaCl", "Wei, 2022", 769.45,
        1, 1, 1, 1, 1, 1, 1,
        ""
    ))

    # 0.15CaCl2-0.315MgCl2-0.535NaCl
    cfgs.append(_cfg(
        "0.15CaCl2-0.315MgCl2-0.535NaCl", "Wei, 2022", 693.98,
        1, 1, 1, 1, 1, 1, 1,
        "0.535NaCl-0.15CaCl2-0.315MgCl2 (Rong, 2020)"
    ))

    # 0.64NaCl-0.36UCl3
    cfgs.append(_cfg(
        "0.64NaCl-0.36UCl3", "Andersson, 2022", 796.15,
        1, 1, 1, 1, 1, 1, 1,
        "63NaCl-37UCl3 (Termini, 2024),65.8NaCl-34.2UCl3 (Rose, 2023)"
    ))

    # 0.5454LiF-0.3636NaF-0.091UF4
    cfgs.append(_cfg(
        "0.5454LiF-0.3636NaF-0.091UF4", "Grizzi, 2024", 923,
        1, 1, 1, 1, 1, 1, 1,
        ""
    ))

    return cfgs
