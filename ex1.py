#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq

# ---------- utilidades mínimas ----------

def translate_all_frames(nt_seq: Seq) -> List[Tuple[str, Seq]]:
    """Traduce las 6 frames (+1,+2,+3 y -1,-2,-3) sin cortar en '*'.
    Devuelve lista de (label, aa_seq)."""
    frames: List[Tuple[str, Seq]] = []
    # forward
    for off in range(3):
        aa = nt_seq[off:].translate(to_stop=False)
        frames.append((f"+{off+1}", aa))
    # reverse
    rc = nt_seq.reverse_complement()
    for off in range(3):
        aa = rc[off:].translate(to_stop=False)
        frames.append((f"-{off+1}", aa))
    return frames

def all_orfs_from_aa(aa_seq: Seq) -> List[Tuple[int, int, Seq]]:
    """Devuelve TODAS las ORFs canónicas M...*
    Cada item: (aa_start_idx, aa_end_idx_exclusive, aa_subseq) en coords de AA (0-based)."""
    text = str(aa_seq)
    out: List[Tuple[int, int, Seq]] = []
    i, n = 0, len(text)
    while i < n:
        if text[i] == 'M':
            j = i
            while j < n and text[j] != '*':
                j += 1
            out.append((i, j, Seq(text[i:j])))
            i += 1  # permitir ORFs solapadas
        else:
            i += 1
    return out

def pick_best_orf(orfs: List[Tuple[int, int, Seq]]) -> Tuple[int, int, Seq] | None:
    """Elige la mejor ORF como la MÁS LARGA."""
    if not orfs:
        return None
    return max(orfs, key=lambda o: len(o[2]))

def wrap60(s: str) -> str:
    return "\n".join(s[i:i+60] for i in range(0, len(s), 60))

# ---------- programa principal ----------

def main():
    ap = argparse.ArgumentParser(
        description="Ejercicio 1: producir FASTA apto para BLAST con ORFs M...*"
    )
    ap.add_argument("-i", "--input", required=True, help="Archivo GenBank (.gb/.gbk) con 1+ mRNAs")
    ap.add_argument("-o", "--output", required=True, help="Archivo FASTA de salida (.faa)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--all-orfs", action="store_true", help="Emitir TODAS las ORFs de las 6 frames")
    mode.add_argument("--best", action="store_true", help="Emitir UNA ORF por registro (la más larga)")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)

    records = list(SeqIO.parse(str(in_path), "genbank"))
    if not records:
        raise SystemExit("No se encontraron registros GenBank en el input.")

    out_lines: List[str] = []

    for rec in records:
        rec_id = rec.id or getattr(rec, "name", "") or "record"
        nt_seq: Seq = rec.seq.upper()

        # Buscar ORFs en las 6 frames
        frames = translate_all_frames(nt_seq)
        all_orfs: List[Tuple[str, int, int, int, Seq]] = []  # (frame_label, idx, a, b, s)

        for label, aa in frames:
            for idx, (a, b, s) in enumerate(all_orfs_from_aa(aa), start=1):
                all_orfs.append((label, idx, a, b, s))

        if not all_orfs:
            # Nada encontrado para este registro; continuar con el siguiente
            continue

        if args.all_orfs:
            # Emitir todas las ORFs (sin '*')
            for label, idx, a, b, s in all_orfs:
                hdr = (f">{rec_id}|frame={label}|orf={idx}"
                       f"|aa_start={a+1}|aa_end={b}|aa_len={len(s)}")
                out_lines.append(hdr)
                out_lines.append(wrap60(str(s)))
        else:
            # Elegir la mejor ORF (más larga) y emitir una sola
            best = pick_best_orf([(a, b, s) for (_l, _i, a, b, s) in all_orfs])
            if not best:
                continue
            a, b, s = best
            # Encontrar su frame para la cabecera
            best_label = next(
                (lbl for (lbl, _i, aa, bb, ss) in all_orfs if aa == a and bb == b and str(ss) == str(s)),
                "NA"
            )
            hdr = f">{rec_id}|best_frame={best_label}|len={len(s)}"
            out_lines.append(hdr)
            out_lines.append(wrap60(str(s)))

    if not out_lines:
        raise SystemExit("No se generó ninguna ORF. Verificá el GenBank.")

    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    nseq = sum(1 for l in out_lines if l.startswith(">"))
    print(f"[OK] Escribí {out_path} con {nseq} secuencia(s) aptas para BLAST.")

if __name__ == "__main__":
    main()
