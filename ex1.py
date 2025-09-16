#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq

def translate_all_frames(nt_seq: Seq, to_stop: bool = False) -> List[Tuple[str, Seq]]:
    """
    Devuelve traducciones para las 6 frames:
      +1,+2,+3 y -1,-2,-3 (reverse complement).
    """
    frames = []
    # forward frames
    for offset in range(3):
        aa = nt_seq[offset:].translate(to_stop=to_stop)
        frames.append((f"+{offset+1}", aa))
    # reverse frames
    rc = nt_seq.reverse_complement()
    for offset in range(3):
        aa = rc[offset:].translate(to_stop=to_stop)
        frames.append((f"-{offset+1}", aa))
    return frames

def longest_orf_from_aa(aa_seq: Seq) -> Optional[Seq]:
    """
    Dado un AA con '*' marcando stops, devuelve el ORF más largo
    desde 'M' hasta antes de '*' (si existe). Si no hay 'M',
    toma el tramo más largo sin '*'.
    """
    text = str(aa_seq)
    best = ""
    i = 0
    n = len(text)
    while i < n:
        if text[i] == 'M':
            j = i
            while j < n and text[j] != '*':
                j += 1
            cand = text[i:j]
            if len(cand) > len(best):
                best = cand
            i = j + 1
        else:
            i += 1
    if not best:
        # fallback: fragmento más largo sin '*'
        for frag in text.split('*'):
            if len(frag) > len(best):
                best = frag
    return Seq(best) if best else None

def extract_cds_translation(record) -> List[Tuple[str, Seq]]:
    """
    Si el GenBank trae anotada la CDS y/o /translation, la agregamos.
    """
    out = []
    if hasattr(record, "features"):
        for feat in record.features:
            if feat.type == "CDS":
                # 1) si hay /translation, la usamos
                tr = feat.qualifiers.get("translation", [])
                if tr:
                    out.append(("cds", Seq(tr[0].replace(" ", "").replace("\n", ""))))
                    continue
                # 2) si no hay /translation, recortamos por localización y traducimos
                try:
                    cds_nt = feat.extract(record.seq)
                    out.append(("cds", cds_nt.translate(to_stop=True)))
                except Exception:
                    pass
    return out

def write_fasta(lines: List[str], out_path: Path):
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def main():
    p = argparse.ArgumentParser(
        description="Ejercicio 1: traducir mRNA(s) GenBank a AA (6 frames) y exportar FASTA."
    )
    p.add_argument("-i", "--input", required=True, help="Archivo GenBank (.gb/.gbk) con uno o más mRNAs")
    p.add_argument("-o", "--output", help="Archivo FASTA de salida (por defecto: <input>.fas)")
    p.add_argument("--longest", action="store_true",
                   help="Exportar solo el ORF más largo por registro (buscado entre las 6 frames)")
    p.add_argument("--include-cds", action="store_true",
                   help="Si hay CDS anotada en GenBank, incluir su traducción como extra")
    p.add_argument("--to-stop", action="store_true",
                   help="Traducir cortando en primer stop por frame (BioPython translate(to_stop=True))")
    args = p.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_suffix(".fas")

    records = list(SeqIO.parse(str(in_path), "genbank"))
    if not records:
        raise SystemExit("No se encontraron registros GenBank en el input.")

    fasta_lines = []

    for rec in records:
        rec_id = rec.id if rec.id else (rec.name or "record")
        seq_nt: Seq = rec.seq.upper()

        frames = translate_all_frames(seq_nt, to_stop=args.to_stop)

        if args.longest:
            # buscar ORF más largo entre las 6 frames
            best_label = None
            best_seq = Seq("")
            for label, aa in frames:
                orf = longest_orf_from_aa(aa)
                if orf and len(orf) > len(best_seq):
                    best_seq = orf
                    best_label = label
            if best_label is None or len(best_seq) == 0:
                # si no encontramos nada decente, como fallback tomamos la +1 completa
                best_label = "+1"
                best_seq = frames[0][1]
            header = f">{rec_id}|frame={best_label}|len={len(best_seq)}"
            fasta_lines.append(header)
            fasta_lines.append(str(best_seq))
        else:
            # exportar todas las frames
            for label, aa in frames:
                header = f">{rec_id}|frame={label}|len={len(aa)}"
                fasta_lines.append(header)
                fasta_lines.append(str(aa))

        if args.include_cds:
            cds_list = extract_cds_translation(rec)
            for label, aa in cds_list:
                header = f">{rec_id}|{label}=annot|len={len(aa)}"
                fasta_lines.append(header)
                fasta_lines.append(str(aa))

    write_fasta(fasta_lines, out_path)
    print(f"[OK] Escribí {out_path} con {sum(1 for l in fasta_lines if l.startswith('>'))} secuencia(s).")

if __name__ == "__main__":
    main()
