#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Ejercicio 2.a - BLAST remoto (NCBI) sobre Swiss-Prot (proteínas)
# Requiere: biopython
# Uso:
#   python ex2_remote.py -i ex1_orfmax.fasta -o blast_remote.xml
#   (para salida tabular: -f 6 -o blast_remote.tsv)

import argparse, time
from pathlib import Path
from Bio.Blast import NCBIWWW, NCBIXML

def run_remote_blast(fasta_path: Path, out_path: Path, outfmt: int = 5, top_hits: int = 10):
    seqs = fasta_path.read_text().strip().split(">")
    seqs = [(">"+s).strip() for s in seqs if s.strip()]
    # Formatos: 5 = XML (estándar BLAST), 6 = tabular (TSV)
    program = "blastp"      # buscamos proteínas (ORFs AA del Ej. 1)
    database = "swissprot"  # Swiss-Prot en NCBI
    with out_path.open("w", encoding="utf-8") as fout:
        for idx, fa in enumerate(seqs, start=1):
            print(f"[{idx}/{len(seqs)}] Enviando a NCBI...", flush=True)
            fmt = "XML" if outfmt == 5 else "Tabular"
            result_handle = NCBIWWW.qblast(program=program,
                                           database=database,
                                           sequence=fa,
                                           format_type=fmt,
                                           hitlist_size=top_hits,
                                           expect=10.0)
            fout.write(result_handle.read())
            fout.write("\n")
            result_handle.close()
            time.sleep(2)

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="BLAST remoto (NCBI) sobre Swiss-Prot")
    p.add_argument("-i", "--input", required=True, help="FASTA AA (Ej.1)")
    p.add_argument("-o", "--output", default="blast_remote.xml", help="Salida (XML o TSV)")
    p.add_argument("-f", "--outfmt", type=int, default=5, choices=[5,6], help="5=XML, 6=Tabular")
    p.add_argument("-n", "--num_hits", type=int, default=10, help="Cantidad de top hits (default 10)")
    args = p.parse_args()
    run_remote_blast(Path(args.input), Path(args.output), args.outfmt, args.num_hits)
    print(f"[OK] Escribí {args.output}")
