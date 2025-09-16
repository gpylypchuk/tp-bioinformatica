#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# msa_report.py

from pathlib import Path
import argparse
from Bio import AlignIO

def percent_identity(s1, s2):
    match = total = 0
    for a,b in zip(s1, s2):
        if a == '-' or b == '-':
            continue
        total += 1
        if a.upper() == b.upper():
            match += 1
    return (100.0 * match / total) if total else 0.0

def consensus_simple(columns):
    cons = []
    for col in columns:
        freq = {}
        for aa in col:
            if aa == '-': continue
            freq[aa.upper()] = freq.get(aa.upper(), 0) + 1
        cons.append(max(freq, key=freq.get) if freq else '-')
    return "".join(cons)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Resumen simple de un MSA (alignment.fasta)")
    ap.add_argument("-i","--input", default="msa_out/alignment.fasta")
    args = ap.parse_args()

    aln = AlignIO.read(args.input, "fasta")
    L = aln.get_alignment_length()
    print(f"Secuencias: {len(aln)} | Largo del alineamiento: {L}")

    # asumir que la primera es QUERY (por cómo armamos el FASTA)
    query = str(aln[0].seq)
    columns = [ [row.seq[i] for row in aln] for i in range(L) ]
    fully_conserved = sum(1 for col in columns
                          if len({aa.upper() for aa in col if aa!='-'}) == 1)
    print(f"Columnas 100% conservadas: {fully_conserved} ({100.0*fully_conserved/L:.1f}%)")

    print("\nIdentidad vs QUERY:")
    for rec in aln[1:]:
        pid = percent_identity(query, str(rec.seq))
        print(f"  {rec.id:40s}  {pid:6.2f}%")

    cons = consensus_simple(columns)
    with open(Path(args.input).with_suffix(".consensus.fa"), "w") as fh:
        fh.write(">CONSENSUS\n")
        fh.write(cons + "\n")
    print(f"\n[OK] Consenso → {Path(args.input).with_suffix('.consensus.fa')}")
