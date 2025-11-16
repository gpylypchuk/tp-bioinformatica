#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejercicio 3 – MSA a partir del XML del Ej. 2 (sin re-BLASTear)
- Lee BLAST XML (outfmt 5) generado en el Ej. 2
- Extrae top-10 hits (no redundantes) por consulta
- Descarga FASTA de los hits (NCBI Protein)
- Crea FASTA: QUERY + top10
- Corre MSA con MUSCLE/Clustal Ω

Uso típico:
  export NCBI_EMAIL="tu_email@dominio.com"
  python ex3_from_xml.py --xml blast_remote.xml --query query.faa --outdir msa_out --msa-tool muscle
"""

import os, sys, time, argparse, subprocess
from pathlib import Path
from typing import List, Tuple, Optional, Set
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

def die(msg: str, code: int = 1):
    print(f"[ERROR] {msg}", file=sys.stderr); sys.exit(code)

def which(cmd: str) -> Optional[str]:
    from shutil import which as _which
    return _which(cmd)

def read_fasta_one(path: Path):
    recs = list(SeqIO.parse(str(path), "fasta"))
    if not recs: die(f"Sin secuencias en {path}")
    if len(recs) > 1: print(f"[WARN] {path} tiene {len(recs)} secuencias; uso la primera ({recs[0].id}).")
    return recs[0]

def write_fasta(path: Path, records: List[SeqIO.SeqRecord]):
    SeqIO.write(records, str(path), "fasta")

def parse_top_hits(xml_path: Path, topn: int = 10) -> List[Tuple[str, str, str, float]]:
    """
    Devuelve hits como lista de tuplas:
      (query_id, accession, title, min_evalue)
    Soporta XML con 1+ resultados concatenados (usa NCBIXML.parse).
    """
    hits = []
    with open(xml_path) as fh:
        for blast_record in NCBIXML.parse(fh):
            qid = getattr(blast_record, "query", None) or getattr(blast_record, "query_id", "QUERY")
            seen: Set[str] = set()
            for aln in blast_record.alignments:
                acc = aln.accession or aln.hit_id.split("|")[-1]
                if acc in seen: 
                    continue
                min_eval = min(hsp.expect for hsp in aln.hsps) if aln.hsps else 1e9
                hits.append((qid, acc, aln.title, min_eval))
                seen.add(acc)
                if len(seen) >= topn:
                    break
    return hits

def fetch_proteins_ncbi(accessions: List[str], email: str) -> List[SeqIO.SeqRecord]:
    """
    Descarga FASTA de proteínas por accession desde NCBI Protein.
    (P68871, etc.). Requiere Entrez.email.
    """
    Entrez.email = email
    out: List[SeqIO.SeqRecord] = []
    for i in range(0, len(accessions), 50):
        chunk = accessions[i:i+50]
        with Entrez.efetch(db="protein", id=",".join(chunk), rettype="fasta", retmode="text") as h:
            out.extend(list(SeqIO.parse(h, "fasta")))
        time.sleep(0.34)  # suave con NCBI
    return out

def run_msa(tool: str, in_fa: Path, out_aln: Path):
    tool = tool.lower()
    if tool == "none":
        print("[MSA] Saltado.")
        return
    if tool == "muscle":
        exe = which("muscle") or which("muscle5") or which("muscle.exe")
        if not exe: die("No se encontró MUSCLE en el PATH. Instalalo o usa --msa-tool none.")
        cmd = [exe, "-align", str(in_fa), "-output", str(out_aln)]
    elif tool in ("clustalo","clustalomega"):
        exe = which("clustalo") or which("clustalo.exe")
        if not exe: die("No se encontró Clustal Ω en el PATH. Instalalo o usa --msa-tool none.")
        cmd = [exe, "-i", str(in_fa), "-o", str(out_aln), "--force", "--threads=2"]
    else:
        die(f"MSA tool desconocido: {tool}")
    print("[MSA CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"[OK] Alineamiento → {out_aln}")

def main():
    ap = argparse.ArgumentParser(description="Ej.3: MSA usando SOLO el XML del Ej.2 (sin re-BLAST)")
    ap.add_argument("--xml", required=True, help="Reporte BLAST XML (outfmt 5) del Ej. 2 (con top-10)")
    ap.add_argument("--query", help="FASTA AA de la consulta para incluir en el MSA (recomendado)")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL",""), help="Email para Entrez (o export NCBI_EMAIL)")
    ap.add_argument("--outdir", default="msa_out", help="Directorio de salida")
    ap.add_argument("--msa-tool", default="muscle", choices=["muscle","clustalo","none"], help="Herramienta de MSA")
    ap.add_argument("--topn", type=int, default=10, help="Cantidad de hits a tomar (por si el XML trae >10)")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    xml_path = Path(args.xml)
    if not xml_path.exists():
        die(f"No existe el XML: {xml_path}")

    # 1) Parsear hits del XML
    hits = parse_top_hits(xml_path, topn=args.topn)
    if not hits:
        die("No hay hits en el XML (¿el Ej.2 devolvió resultados?)")

    # Si el XML tuviera múltiples queries, agrupamos por query_id y procesamos cada una
    from collections import defaultdict
    by_query = defaultdict(list)
    for qid, acc, title, ev in hits:
        by_query[qid].append((acc, title, ev))

    if not args.email:
        die("Seteá --email o export NCBI_EMAIL para Entrez (requisito de NCBI).")

    # 2) Para cada query en el XML, armar FASTA y MSA
    for qid, tuples in by_query.items():
        accessions = [t[0] for t in tuples]
        print(f"[INFO] {qid}: {len(accessions)} hit(s) del XML")

        # Descargar secuencias de los hits
        hit_records = fetch_proteins_ncbi(accessions, email=args.email)

        # Armar FASTA: QUERY (si la pasaron) + hits
        records = []
        if args.query:
            qrec = read_fasta_one(Path(args.query))
            qrec.id = f"QUERY|{qrec.id}"
            qrec.description = "QUERY"
            records.append(qrec)

        # Ordenar los hits en el orden del XML
        rank = {acc:i for i,acc in enumerate(accessions)}
        hit_records.sort(key=lambda r: rank.get(r.id.split("|")[0], 10**9))
        records.extend(hit_records)

        # Guardar FASTA y correr MSA
        fa = outdir / f"{qid.replace(' ','_')}_query_plus_top10.fasta"
        write_fasta(fa, records)
        print(f"[OK] FASTA → {fa}")

        aln = outdir / f"{qid.replace(' ','_')}_alignment.fasta"
        run_msa(args.msa_tool, fa, aln)

if __name__ == "__main__":
    main()
