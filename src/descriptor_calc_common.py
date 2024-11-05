#!/usr/bin/env python

import os
from pathlib import Path
import yaml
import pandas
import re
import pymol


def get_decoys(file_path):
    with open(file_path) as f:
        decoys = [line.strip().split(" ")[0] for line in f.readlines()]
        return decoys


def get_coord_df(ppdb):
    df = pandas.concat([ppdb.df["ATOM"], ppdb.df["HETATM"]])
    return df


def get_chain_id(residue):
    chain_id = residue.split(":")[0]
    return chain_id


def get_wt_res_name(residue: str) -> str:
    """残基書式から Wild Type 残基名を取得する関数
    残基書式：
    X-DEC: [CHAIN]:[WT-RESN][RESI]
      ex): H:F27
    M-DEC: [CHAIN]:[WT-RESN][RESI][MT-RESN]
      ex): H:F27A

    Args:
        residue (str): X-dec, M-dec の残基入力書式
                       * 1残基分の入力書式とする。
    Returns:
        str: 残基書式から取得された Wild Type 残基名
    """
    return residue.split(":")[1][0]


def get_mt_res_name(residue: str) -> str:
    """残基書式から Mutant Type 残基名を取得する関数
    残基書式：
    X-DEC: [CHAIN]:[WT-RESN][RESI]
      ex): H:F27
    M-DEC: [CHAIN]:[WT-RESN][RESI][MT-RESN]
      ex): H:F27A

    Args:
        residue (str): X-dec, M-dec の残基入力書式
                       * 1残基分の入力書式とする。
    Returns:
        str: 残基書式から取得された Mutant Type 残基名
    """
    if re.match(r"[A-Z]", residue[-1]):
        return residue[-1]
    else:
        return residue.split(":")[1][0]


def get_res_no(residue: str) -> int:
    """残基書式から残基番号を取得する関数
    残基書式：
    X-DEC: [CHAIN]:[WT-RESN][RESI]
      ex): H:F27
    M-DEC: [CHAIN]:[WT-RESN][RESI][MT-RESN]
      ex): H:F27A

    Args:
        residue (str): X-dec, M-dec の残基入力書式
                       * 1残基分の入力書式とする。

    Returns:
        int: 残基書式から取得された残基番号
    """
    if re.match(r"[A-Z]", residue[-1]):
        return int(int(residue.split(":")[1][1:-1]))
    else:
        return int(residue.split(":")[1][1:])


def get_res_bone_type(residue: str) -> str:
    """残基書式から骨格情報を取得する関数
    残基書式：
    X-DEC: [CHAIN]:[WT-RESN][RESI]:[BONE]
      ex): H:F27
    M-DEC: [CHAIN]:[WT-RESN][RESI][MT-RESN]:[BONE]
      ex): H:F27A

    Args:
        residue (str): X-dec, M-dec の残基入力書式
                       * 1残基分の入力書式とする。

    Returns:
        int: 残基書式から取得された骨格情報（main/side）
    """
    items = residue.split(":")
    if len(items) == 3:
        return items[2].lower()
    else:
        return None


def get_res_df(df, chain_id, res_no):
    df_res = df[(df["chain_id"] == chain_id) & (df["residue_number"] == res_no)]
    return df_res


def create_result_file(output_dir, file_name):
    output_file_path = os.path.join(output_dir, file_name)
    check_file_path = os.path.join(output_dir, file_name)

    seq = 0
    while os.path.exists(check_file_path):
        seq += 1
        backup_file_name = "".join(["#", "{:0=2}".format(seq), "_", file_name])
        check_file_path = os.path.join(output_dir, backup_file_name)

    if 0 < seq:
        os.rename(output_file_path, check_file_path)

    if not os.path.exists(output_dir):
        Path(output_dir).mkdir()

    Path(output_file_path).touch()
    return output_file_path


def get_decoy_name(decoy):
    name = os.path.splitext(os.path.basename(decoy))[0]
    return name


def convert_atom_name(atom_name, atom_conv_table):
    atom_name_conv = [a[1] for a in atom_conv_table if a[0] == atom_name][0]
    return atom_name_conv


def remove_main_chain(df):
    main_chain = ["N", "CA", "C", "O"]
    df_side_chain = df
    for c in main_chain:
        drop_row = df_side_chain.index[df_side_chain["atom_name"] == c]
        df_side_chain = df_side_chain.drop(drop_row)
    return df_side_chain


def remove_hydrogen(ppdb, atom_conv_table):
    df_atom_exh = ppdb.df["ATOM"]
    df_hetatm_exh = ppdb.df["HETATM"]
    h_all = [a[0] for a in atom_conv_table if a[1] == "H"]
    for h in h_all:
        drop_row_atom = df_atom_exh.index[df_atom_exh["atom_name"] == h]
        df_atom_exh.drop(index=drop_row_atom, inplace=True)
        drop_row_hetatm = df_hetatm_exh.index[df_hetatm_exh["atom_name"] == h]
        df_hetatm_exh.drop(index=drop_row_hetatm, inplace=True)
    return (df_atom_exh, df_hetatm_exh)


def yaml_parse(mol_select_file, pdb, is_mutate=False):
    """相互作用対象分子指定ファイル構文解析関数

    Args:
        mol_select_file (str)): 相互作用対象分子指定ファイル
        pdb (str): PDBファイル
        is_mutate (bool, optional): True: m-dec使用時のフラグ. Defaults to False.
    """
    AMINO = {
        "PHE": "F",
        "LEU": "L",
        "ILE": "I",
        "MET": "M",
        "VAL": "V",
        "SER": "S",
        "PRO": "P",
        "THR": "T",
        "ALA": "A",
        "TYR": "Y",
        "HIS": "H",
        "GLN": "Q",
        "ASN": "N",
        "LYS": "K",
        "ASP": "D",
        "GLU": "E",
        "CYS": "C",
        "TRP": "W",
        "ARG": "R",
        "GLY": "G",
    }
    with open(mol_select_file, "r", encoding="utf8") as file:
        molculars = yaml.safe_load(file)
    pymol.pymol_argv = ["pymol", "-Qc"]
    pymol.finish_launching()
    pymol.cmd.load(pdb)

    res_list = []
    antibody_chain = None
    antigen_chain = None
    for key, val in molculars.items():
        if key.startswith("mutant"):
            chain = None
            if "chain" in val.keys():
                if isinstance(val["chain"], list):
                    assert len(val["chain"]) == 1, "'mutant's 'chain' must be single."
                    chain = val["chain"][0]
                else:
                    assert isinstance(val["chain"], str), "'mutant's 'chain' must be string."
                    chain = val["chain"]

            resi = None
            if "num" in val.keys():
                assert isinstance(val["num"], int), "'mutant's 'num' must be numerical."
                resi = val["num"]

            wt_resn = None
            mt_resn = None
            if "name" in val.keys():
                assert isinstance(val["name"], str), "'mutant's 'name' must be string."
                if is_mutate:
                    mt_resn = val["name"].upper()
                    resn = set()
                    pymol.cmd.iterate(
                        selection=f"resi {resi} & chain {chain}",
                        expression="set.add(resn)",
                        space={"set": resn},
                    )
                    wt_resn = list(resn)[0]
                else:
                    wt_resn = val["name"].upper()

            residue = None
            if is_mutate and resi is not None and mt_resn is not None and chain is not None:
                residue = f"{chain}:{AMINO[wt_resn]}{resi}{AMINO[mt_resn]}"

            elif not is_mutate and resi is not None and wt_resn is not None and chain is not None:
                residue = f"{chain}:{AMINO[wt_resn]}{resi}"

            else:
                raise ValueError("'num', 'name', and 'chain' must be specified for 'mutant'")

            bone = None
            if "type" in val.keys():
                bone = val["type"].lower()
            if bone is not None:
                residue = f"{residue}:{bone}"

            res_list.append(residue)

        elif key.startswith("antibody"):
            if "chain" not in val.keys():
                raise ValueError("'chain' must be specified for 'antibody'")

            antibody_chain = ",".join(val["chain"])

        elif key.startswith("antigen"):
            if "chain" not in val.keys():
                raise ValueError("'chain' must be specified for 'antigen'")

            antigen_chain = ",".join(val["chain"])

    pymol.cmd.reinitialize()
    return ",".join(res_list), antibody_chain, antigen_chain
