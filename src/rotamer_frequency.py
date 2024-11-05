"""rotamer frequency calculation program"""
#!/usr/bin/env python

import math
import os
import sys
import argparse
import traceback
from decimal import Decimal, ROUND_DOWN, ROUND_UP, ROUND_HALF_UP

import pandas
import numpy as np
from Bio.PDB.vectors import Vector, calc_dihedral
from biopandas.pdb import PandasPdb
import yaml

import descriptor_calc_common as dcom

OUTPUT_FILE_NAME = "descriptor_rotamer.csv"
AXIS = ["x_coord", "y_coord", "z_coord"]

np.seterr(all="raise")


class RotamerFrequency:
    """RotamerFrequency"""

    chi_atoms = {}
    rotamer_def = []

    def __init__(
        self,
        residues,
        fpath_model_list,
        fpath_xray,
        fpath_chi_atoms,
        fpath_rotamer_def,
        fpath_rotamer_lib,
        fpath_non_rotameric_def,
        output_dir=os.getcwd(),
    ):
        """Initial values

        Args:
            residues (str): a set of residue(s)
            fpath_model_list (str): a file that specifies input protein structure (model_list.txt)
            fpath_xray (str): path to input protein structure for the reference (e.g., wild type structure)
            fpath_chi_atoms (str): path to heavy atom names for dihedral calculation (../input/chi_atoms.txt)
            fpath_rotamer_def (str): path to rotamer definition (../input/rotamer_definitions.txt)
            fpath_rotamer_lib (str): path to rotamer library (../input/dunbrack2010-library)
            output_dir (str, optional): output directory path . Defaults to os.getcwd().
        """
        self.residues = residues
        self.fpath_model_list = fpath_model_list
        self.fpath_xray = fpath_xray
        self.output_dir = output_dir
        self.fpath_chi_atoms = fpath_chi_atoms
        self.fpath_rotamer_def = fpath_rotamer_def
        self.fpath_rotamer_lib = fpath_rotamer_lib
        self.fpath_non_rotameric_def = fpath_non_rotameric_def
        self.diff_log_rot_freqs = dict()
        self.non_rotameric_def = dict()

    def calculate(self, is_mutate=False):
        """Calculation of rotermer frequency
        Args:
            is_mutate (bool, optional): A flag if the input is mutated compared to the reference structure. Defaults to False.

        Returns:
            str: Rotamer descriptor file name
        """
        try:
            result_file = dcom.create_result_file(self.output_dir, OUTPUT_FILE_NAME)

            self.__load_chi_atoms()
            self.__load_rotamer_definition()
            self.__load_non_rotameric_definition()

            pxray = PandasPdb().read_pdb(self.fpath_xray)
            df_xray = dcom.get_coord_df(pxray)

            decoys = dcom.get_decoys(self.fpath_model_list)

            for decoy in decoys:
                print(f"[RotamerFrequency] {'Decoy':25} : {decoy}")
                pdecoy = PandasPdb().read_pdb(decoy)
                df_decoy = dcom.get_coord_df(pdecoy)

                len_result = None
                if len(self.residues.split(",")) == 1:
                    len_result = 2
                else:
                    len_result = len(self.residues.split(","))

                # Prepare an array that holds values for 
                # the maximum number of Chi (4) x number of target residues.
                chis_xray = [math.nan] * 4 * len_result
                chis_decoy = [math.nan] * 4 * len_result
                diffs = [math.nan] * 4 * len_result

                # Prepare an array for the number of target residues.
                xray_omegas = [math.nan] * len_result
                xray_phis = [math.nan] * len_result
                xray_psis = [math.nan] * len_result
                probs_xray = [math.nan] * len_result
                probs_max_xray = [math.nan] * len_result
                decoy_omegas = [math.nan] * len_result
                decoy_phis = [math.nan] * len_result
                decoy_psis = [math.nan] * len_result
                probs_decoy = [math.nan] * len_result
                probs_max_decoy = [math.nan] * len_result
                sum_d = math.nan if is_mutate else 0
                res_count = 0

                for residue in self.residues.split(","):
                    deg_list_xray = list()
                    deg_list_decoy = list()

                    chain_id = dcom.get_chain_id(residue)
                    res_no = dcom.get_res_no(residue)
                    print(f"[RotamerFrequency] {'Chain id':25} : {chain_id}")
                    print(f"[RotamerFrequency] {'Residue no':25} : {res_no}")

                    # Check the N and C terminal flag
                    n_terminal_flg, c_terminal_flg = self.__get_terminal_flg(
                        df_xray, chain_id, res_no
                    )

                    # DataFrame extraction before and after the target residue number
                    df_res_xray = dcom.get_res_df(df_xray, chain_id, res_no)
                    df_res_xray = pandas.concat(
                        [
                            df_res_xray,
                            dcom.get_res_df(df_xray, chain_id, res_no - 1),
                            dcom.get_res_df(df_xray, chain_id, res_no + 1),
                        ]
                    )
                    df_res_decoy = dcom.get_res_df(df_decoy, chain_id, res_no)
                    df_res_decoy = pandas.concat(
                        [
                            df_res_decoy,
                            dcom.get_res_df(df_decoy, chain_id, res_no - 1),
                            dcom.get_res_df(df_decoy, chain_id, res_no + 1),
                        ]
                    )

                    res_name_xray = df_res_xray.iloc[0]["residue_name"]
                    res_name_decoy = df_res_decoy.iloc[0]["residue_name"]

                    print(f"[RotamerFrequency] {'Residue name(Xray)':25} : {res_name_xray}")
                    if self.__is_ss_cystein(df=df_xray, chain_id=chain_id, resi=res_no):
                        raise ValueError("Not compatible residue(S-S Cystein)")

                    if is_mutate:
                        print(f"[RotamerFrequency] {'Residue name(Mod)':25} : {res_name_decoy}")
                        if self.__is_ss_cystein(df=df_decoy, chain_id=chain_id, resi=res_no):
                            raise ValueError("Not compatible residue(S-S Cystein)")

                    chi_count_xray = None
                    chi_count_decoy = None

                    if res_name_xray not in self.chi_atoms:
                        print("[RotamerFrequency] Rotamer not defined(Xray)")
                        chi_count_xray = 0
                    else:
                        chi_count_xray = int(len(self.chi_atoms[res_name_xray]))

                    if res_name_decoy not in self.chi_atoms:
                        print("[RotamerFrequency] Rotamer not defined(Mod)")
                        chi_count_decoy = 0
                    else:
                        chi_count_decoy = int(len(self.chi_atoms[res_name_decoy]))

                    print(f"[RotamerFrequency] {'Number of chi(Xray)':25} : {chi_count_xray}")
                    print(f"[RotamerFrequency] {'Number of chi(Mod)':25} : {chi_count_decoy}")

                    # Get omege
                    omega_xray, omega_decoy = None, None
                    if res_name_xray == "PRO":
                        omega_xray = self.__calc_omega(n_terminal_flg, df_res_xray, res_no)

                        print(f"[RotamerFrequency] {'Omega (Xray)':25} : {omega_xray:.1f} deg")
                        xray_omegas[res_count] = omega_xray

                    if res_name_decoy == "PRO":
                        omega_decoy = self.__calc_omega(n_terminal_flg, df_res_decoy, res_no)
                        omega_decoy = (
                            self.__round(omega_decoy, "0.1") if n_terminal_flg else omega_decoy
                        )
                        print(f"[RotamerFrequency] {'Omega (Mod)':25} : {omega_decoy:.1f} deg")
                        decoy_omegas[res_count] = omega_decoy

                    # Get phi
                    phi_xray = self.__calc_phi(n_terminal_flg, df_res_xray, res_no)
                    phi_decoy = self.__calc_phi(n_terminal_flg, df_res_decoy, res_no)

                    xray_phis[res_count] = phi_xray
                    decoy_phis[res_count] = phi_decoy
                    print(f"[RotamerFrequency] {'Phi(Xray)':25} : {phi_xray:.1f} deg")
                    print(f"[RotamerFrequency] {'Phi(Mod)':25} : {phi_decoy:.1f} deg")

                    # Get psi
                    psi_xray = self.__calc_psi(c_terminal_flg, df_res_xray, res_no)
                    psi_decoy = self.__calc_psi(c_terminal_flg, df_res_decoy, res_no)

                    xray_psis[res_count] = psi_xray
                    decoy_psis[res_count] = psi_decoy
                    print(f"[RotamerFrequency] {'Psi(Xray)':25} : {psi_xray:.1f} deg")
                    print(f"[RotamerFrequency] {'Psi(Mod)':25} : {psi_decoy:.1f} deg")

                    for idx in range(max(chi_count_xray, chi_count_decoy)):
                        rad_xray = None
                        if idx < chi_count_xray:
                            # Get chi
                            vectors = self.__create_vector(df_res_xray, res_name_xray, idx)
                            rad_xray = calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3])
                            deg_xray = math.degrees(rad_xray)
                            deg_list_xray.append(deg_xray)
                            chis_xray[idx + (4 * res_count)] = deg_xray
                        else:
                            rad_xray = math.nan

                        rad_decoy = None
                        if idx < chi_count_decoy:
                            # Get chi
                            vectors = self.__create_vector(df_res_decoy, res_name_decoy, idx)
                            rad_decoy = calc_dihedral(
                                vectors[0], vectors[1], vectors[2], vectors[3]
                            )
                            deg_decoy = math.degrees(rad_decoy)
                            deg_list_decoy.append(deg_decoy)
                            chis_decoy[idx + (4 * res_count)] = deg_decoy
                        else:
                            rad_decoy = math.nan

                        angle_diff = None
                        if is_mutate:
                            angle_diff = math.nan
                        else:
                            angle_diff = self.__calc_angle_diff(rad_xray, rad_decoy)
                            diffs[idx + (4 * res_count)] = angle_diff
                            sum_d += angle_diff

                        print(
                            f"[RotamerFrequency] {f'Chi_{idx + 1} angle_diff':25} : "
                            f"{angle_diff:.1f}"
                        )

                    if res_name_xray not in ["ALA", "GLY"]:
                        (
                            class_list,
                            class_num_list,
                            df_result_xray,
                            prob_max_xray,
                        ) = self.__get_library_record(
                            res_name=res_name_xray,
                            phi=phi_xray,
                            psi=psi_xray,
                            omega=omega_xray,
                            chis=deg_list_xray,
                        )

                        for idx in range(max(chi_count_xray, chi_count_xray)):
                            label = f"Chi_{idx + 1} class(Xray)"
                            if idx < len(class_list):
                                print(
                                    f"[RotamerFrequency] {label:25} : {class_num_list[idx]}"
                                    f"({class_list[idx]})({deg_list_xray[idx]:.1f} deg)"
                                )
                            else:
                                print(f"[RotamerFrequency] {label:25} : {math.nan}")

                    if res_name_decoy not in ["ALA", "GLY"]:
                        (
                            class_list,
                            class_num_list,
                            df_result_decoy,
                            prob_max_decoy,
                        ) = self.__get_library_record(
                            res_name=res_name_decoy,
                            phi=phi_decoy,
                            psi=psi_decoy,
                            omega=omega_decoy,
                            chis=deg_list_decoy,
                        )

                        for idx in range(max(chi_count_xray, chi_count_decoy)):
                            label = f"Chi_{idx + 1} class(Mod)"
                            if idx < len(class_list):
                                print(
                                    f"[RotamerFrequency] {label:25} : {class_num_list[idx]}"
                                    f"({class_list[idx]})({deg_list_decoy[idx]:.1f} deg)"
                                )
                            else:
                                print(f"[RotamerFrequency] {label:25} : {math.nan}")

                    # Get rotamer probabity
                    if res_name_xray in ["ALA", "GLY"]:
                        probabil_xray, prob_max_xray = 1.0, 1.0

                    else:
                        probabil_xray = self.__get_probability_from_library(
                            res_name_xray, df_result_xray
                        )

                    if probabil_xray == 0.0000001:
                        print(f"[RotamerFrequency] {'Probability(Xray)':25} : prob. < 0.000000")
                    else:
                        print(f"[RotamerFrequency] {'Probability(Xray)':25} : {probabil_xray:.6f}")

                    print(f"[RotamerFrequency] {'Probability.max(Xray)':25} : {prob_max_xray:.6f}")

                    probs_xray[res_count] = probabil_xray
                    probs_max_xray[res_count] = prob_max_xray

                    if res_name_decoy in ["ALA", "GLY"]:
                        probabil_decoy, prob_max_decoy = 1.0, 1.0

                    else:
                        probabil_decoy = self.__get_probability_from_library(
                            res_name_decoy, df_result_decoy
                        )

                    if probabil_decoy == 0.0000001:
                        print(f"[RotamerFrequency] {'Probability(Mod)':25} : prob. < 0.000000")
                    else:
                        print(f"[RotamerFrequency] {'Probability(Mod)':25} : {probabil_decoy:.6f}")
                    print(f"[RotamerFrequency] {'Probability.max(Mod)':25} : {prob_max_decoy:.6f}")

                    probs_decoy[res_count] = probabil_decoy
                    probs_max_decoy[res_count] = prob_max_decoy

                    res_count += 1

                print(f"[RotamerFrequency] {'Angle_diff summary':25} : {sum_d:.1f}")

                probability_xray = self.__calc_probability(
                    [p for p in probs_xray if not math.isnan(p)]
                )
                rot_energy_xray = self.__rotamer_energy(
                    probability_xray, [p for p in probs_max_xray if not math.isnan(p)]
                )

                probability_decoy = self.__calc_probability(
                    [p for p in probs_decoy if not math.isnan(p)]
                )
                rot_energy_decoy = self.__rotamer_energy(
                    probability_decoy, [p for p in probs_max_decoy if not math.isnan(p)]
                )

                print(f"[RotamerFrequency] {'Probability total(Xray)':25} : {probability_xray:.6f}")
                print(f"[RotamerFrequency] {'Probability total(Mod)':25} : {probability_decoy:.6f}")
                print(
                    f"[RotamerFrequency] {'Rotamer Energy(Xray)':25} : "
                    f"{rot_energy_xray:.1f} (kcal/mol)"
                )
                print(
                    f"[RotamerFrequency] {'Rotamer Energy(Mod)':25} : "
                    f"{rot_energy_decoy:.1f} (kcal/mol)"
                )
                print(
                    f"[RotamerFrequency] {'Rotamer Energy_diff(Mod)':25} : "
                    f"{rot_energy_decoy - rot_energy_xray:.1f} (kcal/mol)"
                )

                decoy_name = dcom.get_decoy_name(decoy)
                if any(
                    [
                        probability_decoy == 0,
                        probability_xray == 0,
                        math.isnan(probability_decoy),
                        math.isnan(probability_xray),
                    ]
                ):
                    diff_log_rot_freq = math.pow(10, -9)

                else:
                    diff_log_rot_freq = math.log(probability_decoy / probability_xray)

                self.diff_log_rot_freqs[decoy_name] = diff_log_rot_freq
                print(
                    f"[RotamerFrequency] {'Probability log(Mod/Xray)':25} : "
                    f"{diff_log_rot_freq:.1f}"
                )

                # Output results
                with open(result_file, "a", encoding="utf8") as f_result:
                    columns = [
                        dcom.get_decoy_name(self.fpath_xray),
                        str(self.__round(probability_xray, "0.000001")),
                        str(self.__round(rot_energy_xray, "0.1")),
                        ",".join(map(lambda x: str(self.__round(x, "0.000001")), probs_xray)),
                        ",".join(map(lambda x: str(self.__round(x, "0.000001")), probs_max_xray)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), xray_omegas)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), xray_phis)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), xray_psis)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), chis_xray)),
                        dcom.get_decoy_name(decoy),
                        str(self.__round(probability_decoy, "0.000001")),
                        str(self.__round(rot_energy_decoy, "0.1")),
                        ",".join(map(lambda x: str(self.__round(x, "0.000001")), probs_decoy)),
                        ",".join(map(lambda x: str(self.__round(x, "0.000001")), probs_max_decoy)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), decoy_omegas)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), decoy_phis)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), decoy_psis)),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), chis_decoy)),
                        str(self.__round(rot_energy_decoy - rot_energy_xray, "0.1")),
                        ",".join(map(lambda x: str(self.__round(x, "0.1")), diffs)),
                        str(self.__round(sum_d, "0.1")),
                    ]
                    f_result.write(",".join(columns) + "\n")
            return result_file
        except Exception:
            traceback.print_exc(file=sys.stderr)
            sys.exit(1)

    def __load_chi_atoms(self):
        """Get atom name definitions for chi"""
        with open(self.fpath_chi_atoms, "r", encoding="utf8") as file:
            raw_data = [line.strip().split(" ") for line in file.readlines()]
        for r_d in raw_data:
            key = r_d[0]
            data = [atoms.split("|") for atoms in r_d[1:]]
            self.chi_atoms[key] = data

    def __load_rotamer_definition(self):
        """Get rotamer definitions"""
        with open(self.fpath_rotamer_def, "r", encoding="utf8") as file:
            self.rotamer_def.extend([line.strip().split(" ") for line in file.readlines()])

    def __load_non_rotameric_definition(self):
        """Get non-rotamer definitions"""
        with open(self.fpath_non_rotameric_def, "r", encoding="utf8") as file:
            self.non_rotameric_def = yaml.safe_load(file)

    def __is_ss_cystein(self, df, chain_id, resi):
        """Check SS bond in cystein residues

        Cysteine residues with the both SG atoms located within 2.5 Å of each other are judged to be SS-bonded.

        Args:
            df (dataframe): Dataframe containing structural information
            chain_id (str): Chain ID of the target cysteine residue
            resi (int): Residue number of the target cysteine residue

        Returns:
            bool: True: SS-bonded
        """
        cond_target_res = (df["chain_id"] == chain_id) & (df["residue_number"] == resi)
        sg_coord = (
            df.loc[
                cond_target_res & (df["residue_name"] == "CYS") & (df["atom_name"] == "SG"), AXIS
            ]
            .to_numpy()
            .tolist()
        )
        if len(sg_coord) == 0:
            return False
        sg_coord = sg_coord[0]

        cond_ss_bond = (
            (
                np.sqrt(
                    np.sum(
                        df[["x_coord", "y_coord", "z_coord"]].subtract(sg_coord, axis=1) ** 2,
                        axis=1,
                    )
                )
                <= 2.5
            )
            & (df["atom_name"] == "SG")
            & (df["residue_name"] == "CYS")
            & ~cond_target_res
        )
        return not df.loc[cond_ss_bond].empty

    def __create_vector(self, df, res_name, chi_seq):
        """Vector calculation of chi

        Args:
            df (DataFrame): DataFrame of pdb
            res_name (str): residue name
            chi_seq (int): Number of chi

        Returns:
            list: vector list
        """
        vectors = []
        for atom_name in self.chi_atoms[res_name][chi_seq]:
            xyz = df.loc[df["atom_name"] == atom_name, AXIS].to_numpy().tolist()[0]
            vectors.append(Vector(xyz[0], xyz[1], xyz[2]))
        return vectors

    def __create_vector_proline(self, df, res_no):
        """Vector calculation of Proline (omega)

        Args:
            df (DataFrame): DataFrame of pdb
            res_no (int): Target residue number

        Returns:
            list: vector list
        """
        vectors = []
        xyz = (
            df.loc[(df["atom_name"] == "CA") & (df["residue_number"] == res_no - 1), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "C") & (df["residue_number"] == res_no - 1), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "N") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "CA") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        return vectors

    def __create_vector_phi(self, df, res_no):
        """Vector calculation of phi

        Args:
            df (DataFrame): DataFrame of pdb
            res_no (int): Target residue number

        Returns:
            list: vector list
        """
        vectors = []
        xyz = (
            df.loc[(df["atom_name"] == "C") & (df["residue_number"] == res_no - 1), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "N") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "CA") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "C") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))
        return vectors

    def __create_vector_psi(self, df, res_no):
        """Vector calculation of psi

        Args:
            df (DataFrame): DataFrame of pdb
            res_no (int): Target residue number

        Returns:
            list: vector list
        """
        vectors = []
        xyz = (
            df.loc[(df["atom_name"] == "N") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "CA") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "C") & (df["residue_number"] == res_no), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))

        xyz = (
            df.loc[(df["atom_name"] == "N") & (df["residue_number"] == res_no + 1), AXIS]
            .to_numpy()
            .tolist()
        )[0]
        vectors.append(Vector(xyz[0], xyz[1], xyz[2]))
        return vectors

    def __get_library_record(self, res_name, phi, psi, omega, chis):

        df = self.__get_dunbrack2010_library(res_name, omega)
        # Calculate the bin to which value corresponds to search the rotamer library
        phi_bin = math.floor(math.floor(phi) / 10) * 10
        psi_bin = math.floor(math.floor(psi) / 10) * 10
        df = df[(df["Phi"] == phi_bin) & (df["Psi"] == psi_bin)]
        if df.empty:
            raise ValueError(f"No record matching Phi and Psi: Phi={phi}, Psi={psi}")
        prob_max = df["Probabil"].max()

        non_rotameric_chi = None
        rotameric_defs = [rdef for rdef in self.rotamer_def if rdef[0] == res_name]
        class_list = [None] * len(chis)
        class_num_list = [None] * len(chis)
        for idx, deg in enumerate(chis, start=1):
            if (
                res_name in self.non_rotameric_def.keys()
                and self.non_rotameric_def[res_name] == idx
            ):
                non_rotameric_chi = idx
                class_list[idx - 1] = "-"
                class_num_list[idx - 1] = None
                continue

            # NOTE: Perform the angular transformations to unify the rotor class
            # g+ : [   0,  120)
            # t  : [-180, -120)
            # t  : [ 120,  180)
            # g- : [-120,    0)
            deg = -180 if deg == 180 else deg

            for rotameric in rotameric_defs:
                deg_min = float(rotameric[3])
                deg_max = float(rotameric[4])
                if deg_min <= deg < deg_max:
                    class_list[idx - 1] = rotameric[2]
                    class_num_list[idx - 1] = int(rotameric[1])
                    df = df.loc[df[f"r{idx}"] == int(rotameric[1])]

        if not df.empty and res_name in self.non_rotameric_def.keys():
            # Extract records from non-Rotameric conformations
            df_tmp = df.copy()
            if res_name in ["PHE", "TYR", "GLU", "ASP"]:
                df_tmp["subtraction"] = df[f"chi{non_rotameric_chi}Val"].apply(
                    lambda x: min(
                        abs(x - chis[non_rotameric_chi - 1]),
                        abs(x - chis[non_rotameric_chi - 1] + 180),
                        abs(x - chis[non_rotameric_chi - 1] - 180),
                    )
                )
            else:
                df_tmp["subtraction"] = df[f"chi{non_rotameric_chi}Val"].apply(
                    lambda x: abs(x - chis[non_rotameric_chi - 1])
                )
            df = df[df_tmp["subtraction"] == df_tmp["subtraction"].min()]

            class_num_list[non_rotameric_chi - 1] = df[f"r{non_rotameric_chi}"].iat[0]

        return class_list, class_num_list, df, prob_max

    def __calc_angle_diff(self, rad_xray, rad_mod):
        """Calculation of difference in angle

        Args:
            rad_xray (float): reference structure [rad]
            rad_mod (float): model structure [rad]

        Returns:
            float: difference in angle
        """
        diff = 2 * (1 - math.cos(rad_xray - rad_mod))
        return diff

    def __calc_probability(self, probs):
        """Calculation of the (joint) probability 

        Args:
            probs (list): Rotamer probabilities

        Returns:
            float: the (joint) probability 
        """
        probability = 1
        valid = 0
        for prob in probs:
            if np.isnan(prob):
                continue
            probability *= prob
            valid += 1
        if valid == 0:
            return np.nan
        else:
            return probability

    def __get_probability_from_library(self, res_name, df=None):
        """Get Rotamer probability from the rotamer library

        Args:
            res_name (str): residue name
            df (DataFrame): DataFrame
        Returns:
            float: Rotamer frequency
        """

        # if the rotamer frequency is not defined in the library, return 1e-7 as the probability.
        if df.empty:
            print(
                "[RotamerFrequency] [WARNING] Rotamer undefined in library, "
                "probability is set to 1e-7.",
                file=sys.stderr,
            )
            return 0.0000001

        prob = math.nan
        if any([len(df) == 1, len(df) == 2 and res_name in self.non_rotameric_def.keys()]):
            prob = sum(df["Probabil"].tolist()) / len(df)

        else:
            raise ValueError(f"Inappropriate library\n{df}")

        if prob == 0:
            print(
                "[RotamerFrequency] [WARNING] Probability is 0.000000, " "and set it to 1e-7.",
                file=sys.stderr,
            )
            return 0.0000001

        return prob

    def __get_dunbrack2010_library(self, res_name, omega=None):
        """Get Datafrom from the rotamer library

        Args:
            res_name (str): residue name
            omega (float), optional): omega value. Defaults to None.

        Returns:
            DataFrame:  DataFrame of the rotamer library
        """
        lib_header = {
            "residue_name": [0, 3],
            "Phi": [3, 9],
            "Psi": [9, 14],
            "Count": [14, 20],
            "r1": [20, 26],
            "r2": [26, 29],
            "r3": [29, 32],
            "r4": [32, 35],
            "Probabil": [35, 45],
            "chi1Val": [45, 53],
            "chi2Val": [53, 61],
            "chi3Val": [61, 69],
            "chi4Val": [69, 77],
            "chi1Sig": [77, 87],
            "chi2Sig": [87, 95],
            "chi3Sig": [95, 103],
            "chi4Sig": [103, 111],
        }

        if res_name == "CYS":
            library = os.path.join(self.fpath_rotamer_lib, "cyh.bbdep.rotamers.lib")

        elif res_name != "PRO":
            library = os.path.join(self.fpath_rotamer_lib, f"{res_name.lower()}.bbdep.rotamers.lib")

        elif res_name == "PRO":
            if -90 < omega < 90:
                library = os.path.join(self.fpath_rotamer_lib, "cpr.bbdep.rotamers.lib")

            else:
                library = os.path.join(self.fpath_rotamer_lib, "tpr.bbdep.rotamers.lib")
        else:
            raise ValueError(f"Invalid residue name: {res_name}")

        df = pandas.read_csv(library, comment="#", header=None, names=["0"])
        df_copy = df.copy()

        for key, value in lib_header.items():
            df_copy[key] = df["0"].str[value[0] : value[1]]
            if key == "residue_name":
                df_copy[key] = df_copy[key].astype(str)
            elif key in [
                "Probabil",
                "chi1Val",
                "chi2Val",
                "chi3Val",
                "chi4Val",
                "chi1Sig",
                "chi2Sig",
                "chi3Sig",
                "chi4Sig",
            ]:
                df_copy[key] = df_copy[key].astype(float)
            else:
                df_copy[key] = df_copy[key].astype(int)

        df_copy = df_copy.drop(columns="0")
        return df_copy

    def __rotamer_energy(self, rot_freq, probs_max_list):
        """Convert rotamer probability to energy

        Args:
            rot_freq (float): Rotamer frequency
            probs_max_list (list): The maximum rotamer probability

        Returns:
            float: rotamer energy (kcal/mol)
        """
        R = 0.001987
        T = 300
        rot_energy = 0

        if len(probs_max_list) == 1:
            assert probs_max_list[0] != 0, probs_max_list
            rot_energy += -1 * R * T * math.log(rot_freq / probs_max_list[0])
        else:
            assert probs_max_list[0] != 0 and probs_max_list[1] != 0, probs_max_list
            rot_energy += -1 * R * T * math.log(rot_freq / (probs_max_list[0] * probs_max_list[1]))
        return rot_energy

    def __get_terminal_flg(self, df, chain_id, res_no):
        """Get N and C terminal flag

        Args:
            df (DataFrame): DataFrame of pdb
            chain_id (str): CHAIN ID
            res_no (int): target residue number

        Returns:
            bool: N and C terminal flag
        """
        n_terminal_flg, c_terminal_flg = False, False
        df_res = df[(df["chain_id"] == chain_id) & (df["residue_number"] == res_no - 1)]
        if df_res.empty:
            n_terminal_flg = True
        df_res = df[(df["chain_id"] == chain_id) & (df["residue_number"] == res_no + 1)]
        if df_res.empty:
            c_terminal_flg = True

        return n_terminal_flg, c_terminal_flg

    def __calc_omega(self, n_terminal_flg, df_res, res_no):
        """Calculation of omega

        Args:
            n_terminal_flg (bool): N-terminal flag
            df_res (DataFrame): DataFrame of pdb
            res_no (int): target residue number

        Returns:
            float: omega value [deg]
        """
        if n_terminal_flg:
            return 180
        else:
            vectors = self.__create_vector_proline(df_res, res_no)
            return math.degrees(calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3]))

    def __calc_phi(self, n_terminal_flg, df_res, res_no):
        """Calculation of phi

        Args:
            n_terminal_flg (bool): N-terminal flag
            df_res (DataFrame): DataFrame of pdb
            res_no (int): target residue number

        Returns:
            float: phi value [deg]
        """
        if n_terminal_flg:
            return -60
        else:
            vectors = self.__create_vector_phi(df_res, res_no)
            return math.degrees(calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3]))

    def __calc_psi(self, c_terminal_flg, df_res, res_no):
        """Calculation of psi

        Args:
            c_terminal_flg (bool): C-terminal flag
            df_res (DataFrame): DataFrame of pdb
            res_no (int): target residue number

        Returns:
            float: psi value [deg]
        """
        if c_terminal_flg:
            return 60
        else:
            vectors = self.__create_vector_psi(df_res, res_no)
            return math.degrees(calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3]))

    def __round(self, value, digit="0.1"):
        """Round-off function

        Args:
            value (float): target value
            digit (str, optional): the number of digits to leave after the decimal point.
                                   Defaults to "0.1".

        Returns:
            float: rounded-off value
        """
        value = Decimal(str(value))
        return float(value.quantize(Decimal(digit), rounding=ROUND_HALF_UP))


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser(description=__doc__)

    argument_parser.add_argument(
        "residue_number",
        type=str,
        help="Specify the residue_number.",
    )
    argument_parser.add_argument(
        "model_list",
        type=str,
        help="Specify the model_list path.",
    )
    argument_parser.add_argument(
        "xray",
        type=str,
        help="Specify the xray.",
    )
    argument_parser.add_argument(
        "--chi_atom",
        "-c",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "../input/chi_atoms.txt"),
        help="Specify the chi_atom.",
    )
    argument_parser.add_argument(
        "--rotamer_lib",
        "-l",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "../input/dunbrack2010-library/"),
        help="Specify the rotamer_lib.",
    )
    argument_parser.add_argument(
        "--rotameric_def",
        "-r",
        type=str,
        default=os.path.join(
            os.path.dirname(__file__), "../input/rotamer_definitions-rotameric.txt"
        ),
        help="Specify the rotamer_dif.",
    )
    argument_parser.add_argument(
        "--non_rotameric_def",
        "-n",
        type=str,
        default=os.path.join(
            os.path.dirname(__file__), "../input/rotamer_definitions-nonrotameric.txt"
        ),
        help="Specify the rotamer_lib.",
    )
    argument_parser.add_argument(
        "--is_mutate",
        "-m",
        action="store_true",
        help="is_mutate flag.",
    )

    args = argument_parser.parse_args()

    r = RotamerFrequency(
        args.residue_number,
        args.model_list,
        args.xray,
        args.chi_atom,
        args.rotameric_def,
        args.rotamer_lib,
        args.non_rotameric_def,
    )

    r.calculate(args.is_mutate)
