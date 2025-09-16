#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq

# ---------- utilidades ----------

def translate_all_frames(nt_seq: Seq) -> List[Tuple[str, Seq]]:
    """Traduce las 6 frames (+1,+2,+3 y -1,-2,-3) sin cortar en '*'.
    Devuelve lista de (label, aa_seq)."""
    frames = []
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
    """Devuelve TODAS las ORFs canónicas M...* (sin incluir '*').
    Cada item: (aa_start_idx, aa_end_idx_exclusive, aa_subseq)
    - aa_start_idx y aa_end_idx_exclusive están en coords de aminoácido (0-based)."""
    text = str(aa_seq)
    out = []
    i, n = 0, len(text)
    while i < n:
        if text[i] == 'M':
            j = i
            while j < n and text[j] != '*':
                j += 1
            # segmento M... (hasta antes de '*')
            out.append((i, j, Seq(text[i:j])))
            i += 1  # permitir ORFs solapadas
        else:
            i += 1
    return out

def pick_best_orf(orfs: List[Tuple[int,int,Seq]],
                  length_target: Optional[int] = None,
                  length_minmax: Optional[Tuple[int,int]] = None) -> Optional[Tuple[int,int,Seq]]:
    """Elige la 'mejor' ORF:
       - si length_minmax=(min,max): filtra por rango
       - si length_target: elige la más cercana a target
       - si nada: elige la MÁS LARGA
    """
    if not orfs:
        return None
    candidates = orfs
    if length_minmax:
        mn, mx = length_minmax
        candidates = [o for o in orfs if mn <= len(o[2]) <= mx] or orfs
    if length_target:
        # elegir por proximidad al objetivo
        return sorted(candidates, key=lambda o: abs(len(o[2]) - length_target))[0]
    # por defecto: más larga
    return max(candidates, key=lambda o: len(o[2]))

def extract_cds_aa(record) -> List[Seq]:
    """Si el GenBank trae características CDS, devuelve sus traducciones (usa /translation si está)."""
    aas = []
    if hasattr(record, "features"):
        for feat in record.features:
            if feat.type == "CDS":
                tr = feat.qualifiers.get("translation", [])
                if tr:
                    aas.append(Seq(tr[0].replace(" ", "").replace("\n", "")))
                else:
                    try:
                        cds_nt = feat.extract(record.seq)
                        aas.append(cds_nt.translate(to_stop=True))
                    except Exception:
                        pass
    return aas

def wrap60(s: str) -> str:
    return "\n".join(s[i:i+60] for i in range(0, len(s), 60))

# ---------- programa principal ----------

def main():
    ap = argparse.ArgumentParser(
        description="Ejercicio 1 (revisado): producir FASTA apto para BLAST con ORFs M...*"
    )
    ap.add_argument("-i","--input", required=True, help="Archivo GenBank (.gb/.gbk) con 1+ mRNAs")
    ap.add_argument("-o","--output", help="Archivo FASTA salida (default: <input>.faa)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--all-orfs", action="store_true", help="Emitir TODAS las ORFs M...* de las 6 frames")
    mode.add_argument("--best", action="store_true", help="Emitir UNA ORF por registro (mejor)")
    ap.add_argument("--prefer-cds", action="store_true",
                    help="Si hay CDS anotada, usar su traducción (ignora frames); si no, buscar ORFs en 6 frames")
    ap.add_argument("--length-target", type=int, default=None,
                    help="Longitud objetivo para elegir la mejor ORF (p.ej., 147 para HBB)")
    ap.add_argument("--length-range", type=str, default=None,
                    help="Rango de longitudes 'min,max' para filtrar ORFs (p.ej., 140,160)")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_suffix(".faa")

    length_minmax = None
    if args.length_range:
        try:
            mn, mx = args.length_range.split(",")
            length_minmax = (int(mn), int(mx))
        except Exception:
            raise SystemExit("Formato de --length-range inválido. Ej: 140,160")

    records = list(SeqIO.parse(str(in_path), "genbank"))
    if not records:
        raise SystemExit("No se encontraron registros GenBank en el input.")

    out_lines = []

    for rec in records:
        rec_id = rec.id or rec.name or "record"
        nt_seq: Seq = rec.seq.upper()

        # 1) ¿Usar CDS anotada?
        cds_aas: List[Seq] = extract_cds_aa(rec) if args.prefer_cds else []

        if args.prefer_cds and cds_aas:
            # Emitimos CDS (todas si all-orfs, o la 'mejor' si best)
            if args.all_orfs:
                for k, aa in enumerate(cds_aas, start=1):
                    hdr = f">{rec_id}|source=CDS|idx={k}|len={len(aa)}"
                    out_lines.append(hdr); out_lines.append(wrap60(str(aa)))
            else:
                # elegir una CDS (si hay más de una)
                best_cds = pick_best_orf([(0,len(a),a) for a in cds_aas],
                                         length_target=args.length_target,
                                         length_minmax=length_minmax)
                if not best_cds:
                    continue
                aa = best_cds[2]
                hdr = f">{rec_id}|source=CDS|len={len(aa)}"
                out_lines.append(hdr); out_lines.append(wrap60(str(aa)))
            continue  # saltar búsqueda de frames si usamos CDS

        # 2) Si no usamos/tenemos CDS: buscar ORFs en las 6 frames
        frames = translate_all_frames(nt_seq)
        all_orfs = []
        for label, aa in frames:
            for idx, (a, b, s) in enumerate(all_orfs_from_aa(aa), start=1):
                # Guardamos ORF de esta frame
                all_orfs.append((label, idx, a, b, s))

        if not all_orfs:
            # Fallback: nada encontrado; no emitimos
            continue

        if args.all_orfs:
            # Emitir todas las ORFs (sin '*')
            for label, idx, a, b, s in all_orfs:
                hdr = (f">{rec_id}|frame={label}|orf={idx}"
                       f"|aa_start={a+1}|aa_end={b}|aa_len={len(s)}")
                out_lines.append(hdr); out_lines.append(wrap60(str(s)))
        else:
            # Elegir la mejor ORF
            pick = pick_best_orf([(a,b,s) for (_l,_i,a,b,s) in all_orfs],
                                 length_target=args.length_target,
                                 length_minmax=length_minmax)
            if not pick:
                continue
            # Encontrar su frame para la cabecera
            a,b,s = pick
            label = None
            for _label,_idx,_a,_b,_s in all_orfs:
                if _a==a and _b==b and str(_s)==str(s):
                    label = _label; break
            hdr = f">{rec_id}|best_frame={label}|len={len(s)}"
            out_lines.append(hdr); out_lines.append(wrap60(str(s)))

    if not out_lines:
        raise SystemExit("No se generó ninguna ORF. Verificá el GenBank o los filtros de longitud.")

    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    nseq = sum(1 for l in out_lines if l.startswith(">"))
    print(f"[OK] Escribí {out_path} con {nseq} secuencia(s) aptas para BLAST.")

if __name__ == "__main__":
    main()
