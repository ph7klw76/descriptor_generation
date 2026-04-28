#!/usr/bin/env python3
"""
robust_descriptor_pipeline.py
══════════════════════════════════════════════════════════════════════
Fault-tolerant, restartable, large-scale molecular descriptor pipeline.

PROBLEM IT SOLVES:
    Computing 2870 descriptors for 134,000 QM9 molecules using PaDEL
    (Java), Mordred, and RDKit.  The naive approach crashes because:
      • PaDEL's Java JVM runs out of heap memory on large batches
      • Mordred infinite-loops on certain ring topologies
      • padelpy blocks the process with no timeout mechanism
      • A single crash after 8 hours loses ALL computed data

ARCHITECTURE:
    ┌─────────────────────────────────────────────────────────────┐
    │  SQLite checkpoint database (WAL mode — crash-safe)        │
    │  • One row per molecule                                     │
    │  • Three status flags: mordred_done, rdkit_done, padel_done│
    │  • Descriptor values stored as compressed JSON              │
    │  • Atomic commits — power loss mid-write cannot corrupt     │
    └────────────────────┬────────────────────────────────────────┘
                         │
    ┌────────────────────┴────────────────────────────────────────┐
    │  Stage 1: Mordred    (per-molecule, subprocess timeout)     │
    │  Stage 2: RDKit      (per-molecule, in-process with guard)  │
    │  Stage 3: PaDEL      (micro-batch SDF, subprocess timeout)  │
    │                                                             │
    │  Each stage:                                                │
    │   • Queries DB for pending molecules (stage_done = 0)       │
    │   • Processes them one-by-one or in micro-batches           │
    │   • Commits results to DB after EACH molecule/batch         │
    │   • On crash: restart → skips all completed molecules       │
    │   • On timeout: logs error, marks molecule, moves on        │
    └─────────────────────────────────────────────────────────────┘

INSTALL (Windows Anaconda Prompt):
    conda install -c conda-forge rdkit -y
    pip install mordred padelpy pandas openpyxl
    # Java ≥ 8: https://adoptium.net/ (check "Add to PATH")

USAGE IN SPYDER:
    Edit the paths at the bottom of this file, press F5.
    If it crashes or you stop it: just press F5 again — it resumes.

USAGE FROM COMMAND LINE:
    python robust_descriptor_pipeline.py --qm9 qm9.csv --ref reference.xlsx --out result.csv
"""

import os
import re
import sys
import gc
import json
import time
import zlib
import shutil
import sqlite3
import signal
import logging
import hashlib
import tempfile
import argparse
import traceback
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple, Any
import multiprocessing as mp
from multiprocessing import Process, Queue, Event
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FuturesTimeout
import queue as _queue_mod
import warnings
warnings.filterwarnings("ignore")

# ── Logging ──
LOG_FORMAT = "%(asctime)s │ %(levelname)-5s │ %(message)s"
DATE_FORMAT = "%H:%M:%S"


def setup_logging(log_file: str = "descriptor_pipeline.log"):
    """Configure dual logging: console (INFO) + file (DEBUG)."""
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers.clear()

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter(LOG_FORMAT, DATE_FORMAT))
    root.addHandler(console)

    fh = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s │ %(levelname)-7s │ %(funcName)s │ %(message)s"
    ))
    root.addHandler(fh)

    return logging.getLogger("pipeline")


# ═══════════════════════════════════════════════════════════════════
#  1. CHECKPOINT DATABASE
# ═══════════════════════════════════════════════════════════════════

class CheckpointDB:
    """
    SQLite-based crash-safe checkpoint store.

    Uses WAL (Write-Ahead Logging) journal mode so that:
      • Readers never block writers and vice versa
      • A crash mid-transaction cannot corrupt the database
      • Partial writes are automatically rolled back on recovery

    Descriptor values are stored as zlib-compressed JSON to keep the
    database under ~1 GB even for 134k molecules × 2870 descriptors.
    """

    SCHEMA = """
    CREATE TABLE IF NOT EXISTS molecules (
        name         TEXT PRIMARY KEY,
        smiles       TEXT NOT NULL,
        mordred_done INTEGER DEFAULT 0,
        rdkit_done   INTEGER DEFAULT 0,
        padel_done   INTEGER DEFAULT 0,
        mordred_data BLOB,
        rdkit_data   BLOB,
        padel_data   BLOB,
        error_log    TEXT DEFAULT '',
        created_at   REAL,
        updated_at   REAL
    );
    CREATE INDEX IF NOT EXISTS idx_mordred ON molecules(mordred_done);
    CREATE INDEX IF NOT EXISTS idx_rdkit   ON molecules(rdkit_done);
    CREATE INDEX IF NOT EXISTS idx_padel   ON molecules(padel_done);
    """

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path, timeout=30)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute("PRAGMA cache_size=-64000")  # 64 MB cache
        self.conn.executescript(self.SCHEMA)
        self.conn.commit()

    def register_molecules(self, smiles_dict: Dict[str, str]):
        """Register molecules (skip if already exist)."""
        now = time.time()
        data = [(name, smi, now) for name, smi in smiles_dict.items()]
        self.conn.executemany(
            "INSERT OR IGNORE INTO molecules (name, smiles, created_at) VALUES (?, ?, ?)",
            data
        )
        self.conn.commit()

    def get_pending(self, stage: str, limit: Optional[int] = None) -> List[Tuple[str, str]]:
        """Get molecules not yet completed for a stage."""
        col = f"{stage}_done"
        sql = f"SELECT name, smiles FROM molecules WHERE {col} = 0"
        if limit:
            sql += f" LIMIT {limit}"
        return self.conn.execute(sql).fetchall()

    def get_completed_count(self, stage: str) -> int:
        col = f"{stage}_done"
        row = self.conn.execute(f"SELECT COUNT(*) FROM molecules WHERE {col} = 1").fetchone()
        return row[0]

    def get_total_count(self) -> int:
        row = self.conn.execute("SELECT COUNT(*) FROM molecules").fetchone()
        return row[0]

    def get_failed(self, stage: str) -> List[Tuple[str, str]]:
        """Get molecules that failed (done = -1) for a stage."""
        col = f"{stage}_done"
        return self.conn.execute(
            f"SELECT name, error_log FROM molecules WHERE {col} = -1"
        ).fetchall()

    def save_result(self, name: str, stage: str, descriptors: Dict[str, float]):
        """Save computed descriptors for one molecule, one stage."""
        col_done = f"{stage}_done"
        col_data = f"{stage}_data"
        blob = zlib.compress(json.dumps(descriptors).encode(), level=1)
        self.conn.execute(
            f"UPDATE molecules SET {col_done}=1, {col_data}=?, updated_at=? WHERE name=?",
            (blob, time.time(), name)
        )
        self.conn.commit()

    def save_batch_results(self, results: List[Tuple[str, str, Dict[str, float]]]):
        """Save multiple (name, stage, descriptors) in one transaction."""
        for name, stage, descriptors in results:
            col_done = f"{stage}_done"
            col_data = f"{stage}_data"
            blob = zlib.compress(json.dumps(descriptors).encode(), level=1)
            self.conn.execute(
                f"UPDATE molecules SET {col_done}=1, {col_data}=?, updated_at=? WHERE name=?",
                (blob, time.time(), name)
            )
        self.conn.commit()

    def mark_failed(self, name: str, stage: str, error_msg: str):
        """Mark a molecule as failed for a stage."""
        col = f"{stage}_done"
        self.conn.execute(
            f"UPDATE molecules SET {col}=-1, error_log=error_log||?, updated_at=? WHERE name=?",
            (f"[{stage}] {error_msg}\n", time.time(), name)
        )
        self.conn.commit()

    def load_all_results(self, reference_columns: List[str]) -> pd.DataFrame:
        """
        Reconstruct the full DataFrame from the checkpoint DB.
        Merges mordred + rdkit + padel data for each molecule.

        NOTE: For large datasets (>50k molecules), use
        export_chunked_csv() instead to avoid out-of-memory errors.
        """
        rows = self.conn.execute(
            "SELECT name, mordred_data, rdkit_data, padel_data FROM molecules"
        ).fetchall()

        records = []
        for name, m_blob, r_blob, p_blob in rows:
            row = {"Name": name}
            for blob in [m_blob, r_blob, p_blob]:
                if blob:
                    try:
                        data = json.loads(zlib.decompress(blob))
                        row.update(data)
                    except Exception:
                        pass
            records.append(row)

        df = pd.DataFrame(records)

        # Align to reference columns
        for col in reference_columns:
            if col not in df.columns:
                df[col] = 0.0
        df = df[["Name"] + reference_columns].fillna(0.0)

        return df

    def export_chunked_csv(self, reference_columns: List[str],
                            output_csv: str, chunk_size: int = 5000) -> int:
        """
        Stream results from checkpoint DB to CSV in chunks.
        Never holds more than chunk_size rows in memory.

        This is the memory-safe alternative to load_all_results().

        Returns: total number of rows written
        """
        import csv

        header = ["Name"] + reference_columns
        total = self.get_total_count()
        written = 0

        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(header)

            cursor = self.conn.execute(
                "SELECT name, mordred_data, rdkit_data, padel_data FROM molecules"
            )

            while True:
                rows = cursor.fetchmany(chunk_size)
                if not rows:
                    break

                for name, m_blob, r_blob, p_blob in rows:
                    # Merge all descriptor blobs for this molecule
                    merged = {}
                    for blob in [m_blob, r_blob, p_blob]:
                        if blob:
                            try:
                                data = json.loads(zlib.decompress(blob))
                                merged.update(data)
                            except Exception:
                                pass

                    # Build row aligned to reference columns
                    csv_row = [name]
                    for col in reference_columns:
                        val = merged.get(col, 0.0)
                        try:
                            csv_row.append(float(val) if val and val == val else 0.0)
                        except (ValueError, TypeError):
                            csv_row.append(0.0)

                    writer.writerow(csv_row)
                    written += 1

                # Flush every chunk
                f.flush()

        return written

    def get_progress_summary(self) -> Dict:
        """Get completion statistics."""
        total = self.get_total_count()
        return {
            "total": total,
            "mordred_done": self.get_completed_count("mordred"),
            "mordred_failed": len(self.get_failed("mordred")),
            "rdkit_done": self.get_completed_count("rdkit"),
            "rdkit_failed": len(self.get_failed("rdkit")),
            "padel_done": self.get_completed_count("padel"),
            "padel_failed": len(self.get_failed("padel")),
        }

    def close(self):
        self.conn.close()


# ═══════════════════════════════════════════════════════════════════
#  2. ISOLATED MOLECULE COMPUTERS (run in subprocesses)
# ═══════════════════════════════════════════════════════════════════

def _compute_mordred_single(smiles: str) -> Optional[Dict[str, float]]:
    """
    Compute Mordred descriptors for ONE molecule.
    Runs in a subprocess so a hang/crash cannot kill the main process.
    """
    try:
        from rdkit import Chem
        from mordred import Calculator, descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        calc = Calculator(descriptors, ignore_3D=True)
        result = calc(mol)

        desc = {}
        for descriptor, value in zip(calc.descriptors, result):
            name = str(descriptor)
            if isinstance(value, (int, float, np.integer, np.floating)):
                v = float(value)
                if np.isfinite(v):
                    desc[name] = v
                else:
                    desc[name] = 0.0
            else:
                desc[name] = 0.0

        return desc
    except Exception:
        return None


def _compute_rdkit_single(smiles: str) -> Optional[Dict[str, float]]:
    """Compute RDKit descriptors for ONE molecule."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED
        from rdkit.ML.Descriptors import MoleculeDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        desc_names = [x[0] for x in Descriptors._descList]
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(desc_names)
        values = calculator.CalcDescriptors(mol)

        desc = {}
        for name, val in zip(desc_names, values):
            if isinstance(val, (int, float)) and np.isfinite(val):
                desc[name] = float(val)
            else:
                desc[name] = 0.0

        # Rename to match dataset conventions
        if "TPSA" in desc:
            desc["TPSA_y"] = desc.pop("TPSA")

        try:
            desc["qed"] = float(QED.qed(mol))
        except Exception:
            desc["qed"] = 0.0

        return desc
    except Exception:
        return None


# Default fingerprint groups to enable (padelpy ships them disabled).
# Used when no reference-driven restriction is supplied — produces full
# coverage matching modified_skewness_filtered_data_updated_smile.xlsx.
_PADEL_DEFAULT_ENABLE_GROUPS = (
    "Fingerprinter",
    "ExtendedFingerprinter",
    "GraphOnlyFingerprinter",
    "AtomPairs2DFingerprinter",
    "AtomPairs2DFingerprintCount",
)

# All PaDEL 2D + fingerprint group names (from descriptors.xml).
_PADEL_ALL_2D_GROUPS = (
    "AcidicGroupCount", "ALOGP", "APol", "AromaticAtomsCount",
    "AromaticBondsCount", "AtomCount", "Autocorrelation", "BaryszMatrix",
    "BasicGroupCount", "BCUT", "BondCount", "BPol",
    "BurdenModifiedEigenvalues", "CarbonTypes", "ChiChain", "ChiCluster",
    "ChiPathCluster", "ChiPath", "Constitutional", "Crippen",
    "DetourMatrix", "EccentricConnectivityIndex", "EStateAtomType",
    "ExtendedTopochemicalAtom", "FMF", "FragmentComplexity",
    "HBondAcceptorCount", "HBondDonorCount", "HybridizationRatio",
    "InformationContent", "KappaShapeIndices", "LargestChain",
    "LargestPiSystem", "LongestAliphaticChain", "MannholdLogP",
    "McGowanVolume", "MDE", "MLFER", "PathCount", "PetitjeanNumber",
    "RingCount", "RotatableBondsCount", "RuleOfFive", "Topological",
    "TopologicalCharge", "TopologicalDistanceMatrix", "TPSA", "VABC",
    "VAdjMa", "WalkCount", "Weight", "WeightedPath", "WienerNumbers",
    "XLogP", "ZagrebIndex",
)
_PADEL_ALL_FP_GROUPS = (
    "Fingerprinter", "ExtendedFingerprinter", "EStateFingerprinter",
    "GraphOnlyFingerprinter", "MACCSFingerprinter", "PubchemFingerprinter",
    "SubstructureFingerprinter", "SubstructureFingerprintCount",
    "KlekotaRothFingerprinter", "KlekotaRothFingerprintCount",
    "AtomPairs2DFingerprinter", "AtomPairs2DFingerprintCount",
)

# Map from PaDEL fingerprint family → output column-name regex.
# Used to map reference column names to required PaDEL groups without
# having to invoke PaDEL for each one.
_PADEL_FP_NAME_PATTERNS = {
    "Fingerprinter":                 re.compile(r"^FP\d+$"),
    "ExtendedFingerprinter":         re.compile(r"^ExtFP\d+$"),
    "EStateFingerprinter":           re.compile(r"^EStateFP\d+$"),
    "GraphOnlyFingerprinter":        re.compile(r"^GraphFP\d+$"),
    "MACCSFingerprinter":            re.compile(r"^MACCSFP\d+$"),
    "PubchemFingerprinter":          re.compile(r"^PubchemFP\d+$"),
    "SubstructureFingerprinter":     re.compile(r"^SubFP\d+$"),
    "SubstructureFingerprintCount":  re.compile(r"^SubFPC\d+$"),
    "KlekotaRothFingerprinter":      re.compile(r"^KRFP\d+$"),
    "KlekotaRothFingerprintCount":   re.compile(r"^KRFPC\d+$"),
    "AtomPairs2DFingerprinter":      re.compile(r"^AD2D\d+$"),
    "AtomPairs2DFingerprintCount":   re.compile(r"^APC2D\d+_"),
}


def _build_padel_descriptor_xml(enabled_groups: Optional[List[str]] = None) -> str:
    """Write a customised PaDEL descriptor-types XML.

    enabled_groups : list of PaDEL group names to enable (e.g. "Autocorrelation",
        "Fingerprinter"). Groups not in this list are forced false. If None,
        keep XML defaults and additionally enable the 5 fingerprint groups
        listed in _PADEL_DEFAULT_ENABLE_GROUPS.

    Returns the path to the customised XML file.
    """
    import padelpy
    src = os.path.join(
        os.path.dirname(padelpy.__file__),
        "PaDEL-Descriptor",
        "descriptors.xml",
    )
    with open(src, "r", encoding="utf-8") as f:
        xml = f.read()

    if enabled_groups is None:
        for group in _PADEL_DEFAULT_ENABLE_GROUPS:
            xml = re.sub(
                r'(<Descriptor name="' + re.escape(group) + r'" value=")false(")',
                r"\1true\2",
                xml,
            )
    else:
        wanted = set(enabled_groups)

        def _set(match):
            name = match.group(1)
            return (f'<Descriptor name="{name}" '
                    f'value="{"true" if name in wanted else "false"}"/>')

        xml = re.sub(
            r'<Descriptor name="([^"]+)" value="(?:true|false)"/>',
            _set, xml,
        )

    out_dir = tempfile.mkdtemp(prefix="padel_xml_")
    out_path = os.path.join(out_dir, "descriptors.xml")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(xml)
    return out_path


def _padel_probe_2d_group(group: str, probe_smi: str = "c1ccc(-c2ccccc2)cc1") -> set:
    """Run PaDEL with only one 2D group enabled, return the column names
    it produces. Used by discover_required_padel_groups()."""
    from rdkit import Chem
    from padelpy import padeldescriptor

    tmpdir = tempfile.mkdtemp()
    sdf = os.path.join(tmpdir, "p.sdf")
    out = os.path.join(tmpdir, "p.csv")
    xml = _build_padel_descriptor_xml(enabled_groups=[group])
    try:
        mol = Chem.MolFromSmiles(probe_smi)
        mol.SetProp("_Name", "p")
        w = Chem.SDWriter(sdf); w.write(mol); w.close()
        padeldescriptor(
            mol_dir=sdf, d_file=out, descriptortypes=xml,
            detectaromaticity=True, standardizenitro=True,
            standardizetautomers=True, threads=1, removesalt=True,
            maxruntime=60000, fingerprints=True, d_2d=True,
        )
        if os.path.isfile(out):
            df = pd.read_csv(out)
            return set(df.columns) - {"Name"}
        return set()
    except Exception:
        return set()
    finally:
        try: os.remove(xml)
        except Exception: pass
        try: os.remove(os.path.dirname(xml))
        except Exception: pass
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def discover_required_padel_groups(reference_columns: List[str],
                                    cache_path: Optional[str] = None,
                                    logger=None) -> List[str]:
    """Determine the minimal set of PaDEL descriptor groups whose
    outputs cover the supplied reference column names.

    Fingerprint groups are matched by regex on the column name (cheap).
    2D groups are discovered empirically by probing PaDEL once per
    group on a small molecule. Result is cached as JSON.
    """
    if cache_path and os.path.isfile(cache_path):
        try:
            with open(cache_path) as f:
                groups = json.load(f)
            if isinstance(groups, list) and groups:
                if logger:
                    logger.info(f"[PaDEL] Loaded restricted group set from cache "
                                f"({len(groups)} groups)")
                return groups
        except Exception:
            pass

    ref_set = set(reference_columns)
    needed = set()

    # Map FP families directly via name regex (also handle "_y" rename for nH/nC/nN).
    for group, pat in _PADEL_FP_NAME_PATTERNS.items():
        if any(pat.match(c) for c in ref_set):
            needed.add(group)

    # Probe each 2D group to see which produce names that appear in ref.
    if logger:
        logger.info(f"[PaDEL] Probing {len(_PADEL_ALL_2D_GROUPS)} 2D groups "
                    f"to discover required set (~1 min)…")
    # Add candidates: include _y-stripped versions (PaDEL produces nH not nH_y)
    ref_set_norm = set(ref_set)
    for k in ("nH_y", "nC_y", "nN_y"):
        if k in ref_set_norm:
            ref_set_norm.add(k[:-2])

    for g in _PADEL_ALL_2D_GROUPS:
        outs = _padel_probe_2d_group(g)
        if outs & ref_set_norm:
            needed.add(g)

    groups = sorted(needed)
    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(groups, f)
        except Exception:
            pass
    if logger:
        logger.info(f"[PaDEL] Restricted group set: {len(groups)} groups → {groups}")
    return groups


def _rename_padel_collisions(desc: Dict[str, float]) -> Dict[str, float]:
    """Rename PaDEL keys that collide with Mordred names to *_y to
    match the reference column convention (cf. the existing TPSA→TPSA_y
    rename in the RDKit stage)."""
    for key in ("nH", "nC", "nN"):
        if key in desc:
            desc[key + "_y"] = desc.pop(key)
    return desc


def _compute_padel_batch(smiles_names: List[Tuple[str, str]],
                          timeout_per_mol: int = 30,
                          enabled_groups: Optional[List[str]] = None) -> Dict[str, Dict[str, float]]:
    """
    Compute PaDEL descriptors for a MICRO-BATCH (10-50 molecules).
    Returns {name: {descriptor: value}}.

    If `enabled_groups` is provided, only those PaDEL descriptor groups
    are computed (massive speedup on large molecules — see
    discover_required_padel_groups).
    """
    try:
        from rdkit import Chem
        from padelpy import padeldescriptor

        tmpdir = tempfile.mkdtemp()
        sdf_path = os.path.join(tmpdir, "batch.sdf")
        csv_path = os.path.join(tmpdir, "output.csv")
        xml_path = _build_padel_descriptor_xml(enabled_groups=enabled_groups)

        writer = Chem.SDWriter(sdf_path)
        valid_names = []
        for smi, name in smiles_names:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mol.SetProp("_Name", name)
                writer.write(mol)
                valid_names.append(name)
        writer.close()

        if not valid_names:
            return {}

        max_runtime = timeout_per_mol * len(valid_names)

        padeldescriptor(
            mol_dir=sdf_path,
            d_file=csv_path,
            descriptortypes=xml_path,
            detectaromaticity=True,
            standardizenitro=True,
            standardizetautomers=True,
            threads=2,
            removesalt=True,
            maxruntime=max_runtime,
            fingerprints=True,
            d_2d=True,
        )

        results = {}
        if os.path.isfile(csv_path):
            df = pd.read_csv(csv_path)
            for _, row in df.iterrows():
                name = row.get("Name", "")
                desc = {}
                for col in df.columns:
                    if col != "Name":
                        val = row[col]
                        if pd.notna(val):
                            try:
                                desc[col] = float(val)
                            except (ValueError, TypeError):
                                desc[col] = 0.0
                results[name] = _rename_padel_collisions(desc)

        # Cleanup
        for f in [sdf_path, csv_path, xml_path]:
            if os.path.isfile(f):
                try:
                    os.remove(f)
                except Exception:
                    pass
        for d in [tmpdir, os.path.dirname(xml_path)]:
            try:
                os.rmdir(d)
            except Exception:
                pass

        return results
    except Exception:
        return {}


def _compute_padel_single_fallback(smiles: str,
                                    enabled_groups: Optional[List[str]] = None) -> Optional[Dict[str, float]]:
    """Fallback: compute PaDEL for ONE molecule using padeldescriptor with
    the same custom XML so fingerprint coverage matches the batch path
    (padelpy.from_smiles cannot accept a custom descriptortypes file)."""
    try:
        from rdkit import Chem
        from padelpy import padeldescriptor

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        tmpdir = tempfile.mkdtemp()
        sdf_path = os.path.join(tmpdir, "single.sdf")
        csv_path = os.path.join(tmpdir, "single.csv")
        xml_path = _build_padel_descriptor_xml(enabled_groups=enabled_groups)

        mol.SetProp("_Name", "mol")
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()

        padeldescriptor(
            mol_dir=sdf_path,
            d_file=csv_path,
            descriptortypes=xml_path,
            detectaromaticity=True,
            standardizenitro=True,
            standardizetautomers=True,
            threads=1,
            removesalt=True,
            maxruntime=120000,
            fingerprints=True,
            d_2d=True,
        )

        desc = None
        if os.path.isfile(csv_path):
            df = pd.read_csv(csv_path)
            if len(df) > 0:
                row = df.iloc[0]
                desc = {}
                for col in df.columns:
                    if col == "Name":
                        continue
                    val = row[col]
                    if pd.notna(val):
                        try:
                            desc[col] = float(val)
                        except (ValueError, TypeError):
                            desc[col] = 0.0

        for f in [sdf_path, csv_path, xml_path]:
            if os.path.isfile(f):
                try:
                    os.remove(f)
                except Exception:
                    pass
        for d in [tmpdir, os.path.dirname(xml_path)]:
            try:
                os.rmdir(d)
            except Exception:
                pass

        if desc is None:
            return None
        return _rename_padel_collisions(desc)
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════
#  3. SUBPROCESS EXECUTOR WITH HARD TIMEOUT
# ═══════════════════════════════════════════════════════════════════

def _timeout_worker(func, args, result_path):
    """Worker entry point: run func and pickle result to file.

    Writing to a file (instead of a multiprocessing.Queue) avoids the
    Windows pipe-buffer deadlock where a worker blocks on q.put() of a
    large result while the parent is still inside p.join() and not yet
    draining the queue.
    """
    import pickle
    try:
        result = func(*args)
        payload = ("ok", result)
    except Exception as e:
        payload = ("err", repr(e))
    try:
        with open(result_path, "wb") as f:
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception:
        pass


class TimeoutExecutor:
    """
    Runs a function in a separate process with a hard time limit.

    Uses multiprocessing.Process directly (NOT ProcessPoolExecutor) so
    that on timeout we can call .terminate() / .kill() to actually stop
    a hung worker. ProcessPoolExecutor's context-manager shutdown waits
    for running tasks to finish, which means a hung Java/PaDEL process
    would freeze the entire pipeline despite the future timeout firing.

    Result is transferred via a pickled tempfile (not a Queue) so that
    large descriptor dicts don't deadlock on Windows pipe buffers.
    """

    @staticmethod
    def run(func, args: tuple, timeout_sec: int,
            default=None) -> Any:
        """
        Execute func(*args) with hard timeout.

        Returns the function's return value, or `default` on timeout/error.
        On timeout, the worker process is forcefully terminated.
        """
        import pickle
        fd, result_path = tempfile.mkstemp(prefix="desc_result_", suffix=".pkl")
        os.close(fd)

        ctx = mp.get_context("spawn")
        p = ctx.Process(target=_timeout_worker, args=(func, args, result_path))
        try:
            p.start()
            p.join(timeout_sec)

            if p.is_alive():
                p.terminate()
                p.join(5)
                if p.is_alive():
                    try:
                        p.kill()
                    except Exception:
                        pass
                    p.join(2)
                return default

            try:
                with open(result_path, "rb") as f:
                    status, payload = pickle.load(f)
            except Exception:
                return default

            if status == "ok":
                return payload
            return default
        finally:
            try:
                if os.path.exists(result_path):
                    os.remove(result_path)
            except Exception:
                pass


# ═══════════════════════════════════════════════════════════════════
#  4. PROGRESS TRACKER WITH ETA
# ═══════════════════════════════════════════════════════════════════

class ProgressTracker:
    """Tracks throughput and estimates time remaining."""

    def __init__(self, total: int, stage_name: str):
        self.total = total
        self.stage = stage_name
        self.start_time = time.time()
        self.processed = 0
        self.last_report = 0

    def update(self, n: int = 1):
        self.processed += n

    def should_report(self, interval: int = 50) -> bool:
        if self.processed - self.last_report >= interval:
            self.last_report = self.processed
            return True
        return False

    def report(self, logger) -> str:
        elapsed = time.time() - self.start_time
        rate = self.processed / max(elapsed, 0.01)
        remaining = self.total - self.processed
        eta_sec = remaining / max(rate, 0.001)
        eta_str = str(timedelta(seconds=int(eta_sec)))

        pct = self.processed / max(self.total, 1) * 100
        bar_len = 30
        filled = int(bar_len * self.processed / max(self.total, 1))
        bar = "█" * filled + "░" * (bar_len - filled)

        msg = (f"[{self.stage}] {bar} {pct:5.1f}% │ "
               f"{self.processed:,}/{self.total:,} │ "
               f"{rate:.1f} mol/s │ ETA {eta_str}")
        logger.info(msg)
        return msg


# ═══════════════════════════════════════════════════════════════════
#  5. STAGE RUNNERS
# ═══════════════════════════════════════════════════════════════════

def run_mordred_stage(db: CheckpointDB, logger,
                      timeout_per_mol: int = 120,
                      report_every: int = 50):
    """
    Stage 1: Compute Mordred descriptors for all pending molecules.

    Each molecule is processed in a subprocess with a hard timeout.
    Results are committed to SQLite after each molecule.
    """
    pending = db.get_pending("mordred")
    if not pending:
        logger.info("[Mordred] All molecules already completed — skipping")
        return

    logger.info(f"[Mordred] {len(pending)} molecules pending")
    tracker = ProgressTracker(len(pending), "Mordred")

    for name, smiles in pending:
        result = TimeoutExecutor.run(
            _compute_mordred_single, (smiles,),
            timeout_sec=timeout_per_mol,
            default=None,
        )

        if result is not None and len(result) > 0:
            db.save_result(name, "mordred", result)
        else:
            db.mark_failed(name, "mordred", f"timeout({timeout_per_mol}s) or invalid")

        tracker.update()
        if tracker.should_report(report_every):
            tracker.report(logger)

        # Periodic garbage collection
        if tracker.processed % 200 == 0:
            gc.collect()

    tracker.report(logger)
    logger.info(f"[Mordred] Stage complete: "
                f"{db.get_completed_count('mordred')} done, "
                f"{len(db.get_failed('mordred'))} failed")


def run_rdkit_stage(db: CheckpointDB, logger,
                    timeout_per_mol: int = 60,
                    report_every: int = 100):
    """
    Stage 2: Compute RDKit descriptors for all pending molecules.

    RDKit is generally stable, so we run in-process with a try/except
    guard rather than a subprocess.  Falls back to subprocess on error.
    """
    pending = db.get_pending("rdkit")
    if not pending:
        logger.info("[RDKit] All molecules already completed — skipping")
        return

    logger.info(f"[RDKit] {len(pending)} molecules pending")
    tracker = ProgressTracker(len(pending), "RDKit")

    for name, smiles in pending:
        try:
            result = _compute_rdkit_single(smiles)
        except Exception as e:
            # Fallback to subprocess if in-process fails
            result = TimeoutExecutor.run(
                _compute_rdkit_single, (smiles,),
                timeout_sec=timeout_per_mol,
                default=None,
            )

        if result is not None and len(result) > 0:
            db.save_result(name, "rdkit", result)
        else:
            db.mark_failed(name, "rdkit", "computation error")

        tracker.update()
        if tracker.should_report(report_every):
            tracker.report(logger)

    tracker.report(logger)
    logger.info(f"[RDKit] Stage complete: "
                f"{db.get_completed_count('rdkit')} done, "
                f"{len(db.get_failed('rdkit'))} failed")


def run_padel_stage(db: CheckpointDB, logger,
                    batch_size: int = 25,
                    timeout_per_mol: int = 45,
                    report_every: int = 25,
                    enabled_groups: Optional[List[str]] = None):
    """
    Stage 3: Compute PaDEL descriptors in micro-batches.

    enabled_groups : if provided, only these PaDEL descriptor groups are
        computed. Use discover_required_padel_groups() to derive a minimal
        set from the reference column list — yields large speedups
        (~50x on 50-atom polyaromatics) and avoids Java timeouts in
        expensive descriptor classes.
    """
    pending = db.get_pending("padel")
    if not pending:
        logger.info("[PaDEL] All molecules already completed — skipping")
        return

    logger.info(f"[PaDEL] {len(pending)} molecules pending (batch_size={batch_size})")
    if enabled_groups:
        logger.info(f"[PaDEL] Restricted to {len(enabled_groups)} groups: {enabled_groups}")
    tracker = ProgressTracker(len(pending), "PaDEL")

    # Process in micro-batches
    for batch_start in range(0, len(pending), batch_size):
        batch = pending[batch_start:batch_start + batch_size]
        batch_pairs = [(smi, name) for name, smi in batch]
        batch_names = [name for name, _ in batch]

        batch_timeout = timeout_per_mol * len(batch) + 60  # extra buffer

        # Try batch processing
        results = TimeoutExecutor.run(
            _compute_padel_batch,
            (batch_pairs, timeout_per_mol, enabled_groups),
            timeout_sec=batch_timeout,
            default=None,
        )

        if results is not None:
            # Save successful batch results
            batch_results = []
            succeeded_names = set()
            for name in batch_names:
                if name in results and results[name]:
                    batch_results.append((name, "padel", results[name]))
                    succeeded_names.add(name)

            if batch_results:
                db.save_batch_results(batch_results)

            # Retry failed molecules individually
            failed_in_batch = [
                (name, smi) for name, smi in batch
                if name not in succeeded_names
            ]
        else:
            # Entire batch failed — retry all individually
            failed_in_batch = batch
            logger.warning(f"[PaDEL] Batch failed, retrying {len(batch)} individually")

        # Individual fallback for failed molecules
        for name, smiles in failed_in_batch:
            result = TimeoutExecutor.run(
                _compute_padel_single_fallback, (smiles, enabled_groups),
                timeout_sec=timeout_per_mol + 30,
                default=None,
            )
            if result is not None and len(result) > 0:
                db.save_result(name, "padel", result)
            else:
                db.mark_failed(name, "padel", "batch+individual failed")

        tracker.update(len(batch))
        if tracker.should_report(report_every):
            tracker.report(logger)

        # Force garbage collection between batches
        gc.collect()

    tracker.report(logger)
    logger.info(f"[PaDEL] Stage complete: "
                f"{db.get_completed_count('padel')} done, "
                f"{len(db.get_failed('padel'))} failed")


# ═══════════════════════════════════════════════════════════════════
#  6. EXPORT: CHECKPOINT DB → ALIGNED CSV
# ═══════════════════════════════════════════════════════════════════

def export_results(db: CheckpointDB, reference_columns: List[str],
                   output_csv: str, logger):
    """
    Export checkpoint DB to CSV using memory-safe chunked streaming.
    Writes 5000 rows at a time — never loads full dataset into RAM.
    """
    logger.info(f"Exporting results to {output_csv}")
    logger.info(f"  Using chunked export (memory-safe for large datasets)")

    total = db.get_total_count()
    logger.info(f"  Total molecules: {total:,}")

    n_written = db.export_chunked_csv(reference_columns, output_csv, chunk_size=5000)

    logger.info(f"  Written: {n_written:,} rows × {len(reference_columns)} columns")
    logger.info(f"  Saved: {output_csv}")

    return None


def export_partial(db: CheckpointDB, reference_columns: List[str],
                   output_csv: str, logger):
    """Export whatever has been computed so far (for use during long runs)."""
    logger.info(f"Exporting partial results...")
    return export_results(db, reference_columns, output_csv, logger)


# ═══════════════════════════════════════════════════════════════════
#  7. MASTER PIPELINE
# ═══════════════════════════════════════════════════════════════════

def run_pipeline(
    qm9_csv: str,
    reference_xlsx: str,
    output_csv: str,
    checkpoint_db: str = "descriptor_checkpoint.db",
    max_molecules: Optional[int] = None,
    mordred_timeout: int = 120,
    rdkit_timeout: int = 60,
    padel_timeout: int = 45,
    padel_batch_size: int = 25,
    skip_mordred: bool = False,
    skip_rdkit: bool = False,
    skip_padel: bool = False,
    padel_restrict_to_ref: bool = False,
    export_interval: int = 5000,
):
    """
    Master pipeline: orchestrates all three stages with full fault tolerance.

    Parameters
    ----------
    qm9_csv : str
        CSV with SMILES column.
    reference_xlsx : str
        MR-TADF descriptor file (for column alignment).
    output_csv : str
        Final output path.
    checkpoint_db : str
        SQLite checkpoint file (created automatically).
        THIS IS YOUR SAFETY NET — do not delete between runs.
    max_molecules : int
        Limit for testing (None = all).
    mordred_timeout : int
        Seconds per molecule for Mordred stage.
    rdkit_timeout : int
        Seconds per molecule for RDKit stage.
    padel_timeout : int
        Seconds per molecule for PaDEL stage.
    padel_batch_size : int
        Molecules per PaDEL micro-batch (lower = safer, slower).
    skip_mordred/rdkit/padel : bool
        Skip a stage entirely (useful for debugging).
    export_interval : int
        Export partial CSV every N molecules (0 = disabled).
    """
    logger = setup_logging()
    t_start = time.time()

    logger.info("═" * 65)
    logger.info("  FAULT-TOLERANT DESCRIPTOR PIPELINE")
    logger.info("═" * 65)
    logger.info(f"  QM9 file:      {qm9_csv}")
    logger.info(f"  Reference:     {reference_xlsx}")
    logger.info(f"  Output:        {output_csv}")
    logger.info(f"  Checkpoint DB: {checkpoint_db}")
    logger.info(f"  Max molecules: {max_molecules or 'ALL'}")
    logger.info("")

    # ── Load reference columns ──
    ref_df = pd.read_excel(reference_xlsx, nrows=0)
    ref_columns = [c for c in ref_df.columns if c not in ["Name", "smile"]]
    logger.info(f"Reference columns: {len(ref_columns)}")

    # ── Load QM9 SMILES ──
    try:
    	qm9_df = pd.read_csv(qm9_csv, encoding="utf-8")
    except UnicodeDecodeError:
    	try:
        		qm9_df = pd.read_csv(qm9_csv, encoding="latin1")
    	except UnicodeDecodeError:
        		qm9_df = pd.read_csv(qm9_csv, encoding="cp1252")

    smi_col = None
    for c in ["smiles", "SMILES", "Smiles", "canonical_smiles", "smi"]:
        if c in qm9_df.columns:
            smi_col = c
            break
    if smi_col is None:
        smi_col = qm9_df.columns[0]

    name_col = None
    for c in ["name", "Name", "mol_id", "id", "gdb_idx"]:
        if c in qm9_df.columns:
            name_col = c
            break

    smiles_list = qm9_df[smi_col].dropna().astype(str).tolist()
    if name_col:
        names = qm9_df[name_col].dropna().astype(str).tolist()
    else:
        names = [f"qm9_{i:06d}" for i in range(len(smiles_list))]

    min_len = min(len(smiles_list), len(names))
    smiles_list, names = smiles_list[:min_len], names[:min_len]

    if max_molecules:
        smiles_list = smiles_list[:max_molecules]
        names = names[:max_molecules]

    logger.info(f"QM9 molecules: {len(smiles_list)}")

    # ── Initialise checkpoint DB ──
    db = CheckpointDB(checkpoint_db)
    smiles_dict = dict(zip(names, smiles_list))
    db.register_molecules(smiles_dict)

    # Show resume status
    progress = db.get_progress_summary()
    logger.info(f"\nResume status:")
    logger.info(f"  Total registered: {progress['total']}")
    logger.info(f"  Mordred:  {progress['mordred_done']} done, {progress['mordred_failed']} failed")
    logger.info(f"  RDKit:    {progress['rdkit_done']} done, {progress['rdkit_failed']} failed")
    logger.info(f"  PaDEL:    {progress['padel_done']} done, {progress['padel_failed']} failed")

    if progress['mordred_done'] + progress['rdkit_done'] + progress['padel_done'] > 0:
        logger.info("  → Resuming from checkpoint (completed molecules will be skipped)")

    # ── Stage 1: Mordred ──
    if not skip_mordred:
        logger.info(f"\n{'━'*65}")
        logger.info("  STAGE 1/3: MORDRED DESCRIPTORS")
        logger.info(f"{'━'*65}")
        run_mordred_stage(db, logger, timeout_per_mol=mordred_timeout)
    else:
        logger.info("\n[Mordred] Skipped by user")

    # ── Stage 2: RDKit ──
    if not skip_rdkit:
        logger.info(f"\n{'━'*65}")
        logger.info("  STAGE 2/3: RDKIT DESCRIPTORS")
        logger.info(f"{'━'*65}")
        run_rdkit_stage(db, logger, timeout_per_mol=rdkit_timeout)
    else:
        logger.info("\n[RDKit] Skipped by user")

    # ── Stage 3: PaDEL ──
    if not skip_padel:
        logger.info(f"\n{'━'*65}")
        logger.info("  STAGE 3/3: PADEL DESCRIPTORS + FINGERPRINTS")
        logger.info(f"{'━'*65}")

        padel_groups = None
        if padel_restrict_to_ref:
            cache_path = checkpoint_db + ".padel_groups.json"
            padel_groups = discover_required_padel_groups(
                ref_columns, cache_path=cache_path, logger=logger,
            )

        run_padel_stage(
            db, logger,
            batch_size=padel_batch_size,
            timeout_per_mol=padel_timeout,
            enabled_groups=padel_groups,
        )
    else:
        logger.info("\n[PaDEL] Skipped by user")

    # ── Final export ──
    logger.info(f"\n{'━'*65}")
    logger.info("  FINAL EXPORT")
    logger.info(f"{'━'*65}")
    export_results(db, ref_columns, output_csv, logger)

    # ── Summary ──
    elapsed = time.time() - t_start
    progress = db.get_progress_summary()

    logger.info(f"\n{'═'*65}")
    logger.info(f"  PIPELINE COMPLETE")
    logger.info(f"{'═'*65}")
    logger.info(f"  Total time: {timedelta(seconds=int(elapsed))}")
    logger.info(f"  Molecules:  {progress['total']}")
    logger.info(f"  Mordred:    {progress['mordred_done']} OK, {progress['mordred_failed']} failed")
    logger.info(f"  RDKit:      {progress['rdkit_done']} OK, {progress['rdkit_failed']} failed")
    logger.info(f"  PaDEL:      {progress['padel_done']} OK, {progress['padel_failed']} failed")
    logger.info(f"  Output:     {output_csv}")
    logger.info(f"  Checkpoint: {checkpoint_db} (keep for debugging)")
    logger.info(f"{'═'*65}")

    db.close()
    logger.info(f"Done! Output saved to: {output_csv}")
    return output_csv


# ═══════════════════════════════════════════════════════════════════
#  RUN IN SPYDER: EDIT PATHS AND PRESS F5
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    # Check if running from command line with arguments
    if len(sys.argv) > 1 and sys.argv[1].startswith("--"):
        parser = argparse.ArgumentParser(
            description="Fault-tolerant molecular descriptor pipeline"
        )
        parser.add_argument("--qm9", required=True, help="QM9 SMILES CSV")
        parser.add_argument("--ref", required=True, help="MR-TADF reference xlsx")
        parser.add_argument("--out", default="qm9_descriptors.csv", help="Output CSV")
        parser.add_argument("--db", default="descriptor_checkpoint.db", help="Checkpoint DB")
        parser.add_argument("--max", type=int, default=None, help="Max molecules")
        parser.add_argument("--padel-batch", type=int, default=25, help="PaDEL batch size")
        parser.add_argument("--mordred-timeout", type=int, default=120, help="Seconds per molecule for Mordred")
        parser.add_argument("--rdkit-timeout", type=int, default=60, help="Seconds per molecule for RDKit")
        parser.add_argument("--padel-timeout", type=int, default=45, help="Seconds per molecule for PaDEL")
        parser.add_argument("--skip-mordred", action="store_true")
        parser.add_argument("--skip-rdkit", action="store_true")
        parser.add_argument("--skip-padel", action="store_true")
        parser.add_argument("--padel-restrict-to-ref", action="store_true",
                            help="Auto-discover and enable only the PaDEL "
                                 "descriptor groups whose outputs appear in "
                                 "the reference xlsx. Big speedup on large "
                                 "molecules; first run probes ~55 groups (~1 min).")
        args = parser.parse_args()

        run_pipeline(
            qm9_csv=args.qm9,
            reference_xlsx=args.ref,
            output_csv=args.out,
            checkpoint_db=args.db,
            max_molecules=args.max,
            mordred_timeout=args.mordred_timeout,
            rdkit_timeout=args.rdkit_timeout,
            padel_timeout=args.padel_timeout,
            padel_batch_size=args.padel_batch,
            skip_mordred=args.skip_mordred,
            skip_rdkit=args.skip_rdkit,
            skip_padel=args.skip_padel,
            padel_restrict_to_ref=args.padel_restrict_to_ref,
        )

    else:
        # ╔═══════════════════════════════════════════════════════╗
        # ║  SPYDER MODE: EDIT THESE PATHS, THEN PRESS F5        ║
        # ║                                                       ║
        # ║  If it crashes or you stop it: JUST PRESS F5 AGAIN.  ║
        # ║  It will resume from where it left off.               ║
        # ╚═══════════════════════════════════════════════════════╝

        result = run_pipeline(
            # ── File paths (EDIT THESE) ──
            qm9_csv=r"C:\Users\Woon\Documents\DICC\Inverteddesign\Book4.csv",
            reference_xlsx=r"C:\Users\Woon\Documents\DICC\Inverteddesign\modified_skewness_filtered_data_updated_smile.xlsx",
            output_csv=r"C:\Users\Woon\Documents\DICC\Inverteddesign\book4_descriptors.csv",

            # ── Checkpoint file (keeps progress across crashes) ──
            checkpoint_db=r"C:\Users\Woon\Documents\DICC\Inverteddesign\descriptor_checkpoint.db",

            # ── Settings ──
            max_molecules=100,      # START WITH 100 TO TEST, then set None for all
            mordred_timeout=120,    # seconds per molecule (increase for slow machines)
            rdkit_timeout=60,       # seconds per molecule
            padel_timeout=45,       # seconds per molecule
            padel_batch_size=25,    # molecules per PaDEL batch (lower = safer)

            # ── Skip stages (set True to skip) ──
            skip_mordred=False,
            skip_rdkit=False,
            skip_padel=False,
        )

        if result is not None:
            print(f"\nDone! Output: {result}")
