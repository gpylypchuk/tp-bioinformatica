#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BLAST Unificado - Remoto y Local
# Uso:
#   python blast.py --mode remote -i secuencias.fasta -o resultado.xml
#   python blast.py --mode local -i secuencias.fasta -o resultado.xml

import argparse
import subprocess
import tempfile
import time
from pathlib import Path
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

def run_remote_blast(fasta_path: Path, out_path: Path, outfmt: int = 5, top_hits: int = 10):
    """
    Ejecuta BLAST remoto contra NCBI Swiss-Prot
    """
    print("🌐 Ejecutando BLAST REMOTO contra NCBI...")
    
    seqs = fasta_path.read_text().strip().split(">")
    seqs = [(">"+s).strip() for s in seqs if s.strip()]
    
    program = "blastp"
    database = "swissprot"
    
    with out_path.open("w", encoding="utf-8") as fout:
        for idx, fa in enumerate(seqs, start=1):
            print(f"[{idx}/{len(seqs)}] Enviando a NCBI...", flush=True)
            fmt = "XML" if outfmt == 5 else "Tabular"
            
            result_handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=fa,
                format_type=fmt,
                hitlist_size=top_hits,
                expect=10.0
            )
            
            fout.write(result_handle.read())
            fout.write("\n")
            result_handle.close()
            time.sleep(2)  # Respeto a los servidores de NCBI

def run_local_blast(fasta_path: Path, out_path: Path, outfmt: int = 5, top_hits: int = 10, db_path: str = "swissprot"):
    """
    Ejecuta BLAST local contra base de datos local - robusto.
    Crea un archivo temporal por query y lo concatena al out_path final.
    """
    print("💻 Ejecutando BLAST LOCAL (robusto)...")

    # Normalize db_path: allow "swissprot", "swissprot.fasta", "/path/to/swissprot.fasta"
    db_p = Path(db_path).resolve()
    if db_p.suffix in {".fasta", ".fa", ".fas", ".gz"}:
        db_stem = db_p.stem
        db_dir = db_p.parent
    else:
        db_stem = db_p.name
        db_dir = db_p.parent if db_p.parent != Path("") else Path.cwd()

    # prefer absolute base path for -db (no extension)
    db_full = str((db_dir / db_stem).resolve())

    # Verify DB files exist (phr/pin/psq)
    expected = [(db_dir / f"{db_stem}.{ext}").exists() for ext in ("phr", "pin", "psq")]
    if not all(expected):
        raise FileNotFoundError(f"BLAST DB files not found for base '{db_full}'. "
                                f"Checked: {(db_dir / f'{db_stem}.phr')}, {(db_dir / f'{db_stem}.pin')}, {(db_dir / f'{db_stem}.psq')}")

    print(f"Usando DB base: {db_full}")
    print(f"DB directory: {db_dir}")

    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    if not sequences:
        raise ValueError(f"No sequences found in {fasta_path}")
    print(f"Procesando {len(sequences)} secuencias...")

    # Ensure output parent dir exists
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # We'll write results incrementally into out_path (overwrite at start)
    with out_path.open("w", encoding="utf-8") as fout:
        fout.write("")  # truncate

    for idx, record in enumerate(sequences, start=1):
        print(f"[{idx}/{len(sequences)}] {record.id} ...", flush=True)

        # write each query to a temp file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qtmp:
            SeqIO.write(record, qtmp.name, "fasta")
            qtmp_path = Path(qtmp.name)

        # result goes to a temp file to avoid stdout/cwd weirdness
        with tempfile.NamedTemporaryFile(mode="w+b", suffix=".out", delete=False) as rtmp:
            rtmp_path = Path(rtmp.name)

        cmd = [
            "blastp",
            "-query", str(qtmp_path.resolve()),
            "-db", db_full,            # base name (no extension)
            "-out", str(rtmp_path.resolve()),
            "-outfmt", str(outfmt),
            "-max_target_seqs", str(top_hits),
            "-evalue", "10.0"
        ]

        print("Ejecutando:", " ".join(cmd))
        # Run without changing cwd (we use full path db)
        proc = subprocess.run(cmd, capture_output=True, text=True)

        # Debug info
        print(f"Return code: {proc.returncode}")
        if proc.stderr:
            print("blast stderr:", proc.stderr.strip())

        if proc.returncode != 0:
            print(f"Error en blastp para {record.id} (returncode {proc.returncode}). Se omitirá esta secuencia.")
            # append stderr details to output file for debugging
            with out_path.open("a", encoding="utf-8") as fout:
                fout.write(f"\n# ERROR for {record.id}\n")
                fout.write(proc.stderr or "No stderr\n")
        else:
            # read the temporary result file and append to final out
            out_text = rtmp_path.read_text(encoding="utf-8")
            if not out_text.strip():
                print(f"⚠️  Resultado vacío para {record.id} (posible sin hits por evalue/short seq).")
                # still append a comment so user knows a query was run
                with out_path.open("a", encoding="utf-8") as fout:
                    fout.write(f"\n# No hits for {record.id}\n")
            else:
                with out_path.open("a", encoding="utf-8") as fout:
                    fout.write(out_text)
                    # Add newline separator in case outfmt is tabular
                    fout.write("\n")

        # cleanup temps
        try:
            qtmp_path.unlink(missing_ok=True)
            rtmp_path.unlink(missing_ok=True)
        except Exception:
            pass

    print("Local BLAST finalizado.")


def check_blast_installation():
    """Verifica si BLAST+ está instalado y disponible"""
    try:
        result = subprocess.run(["blastp", "-version"], 
                               capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False

def check_database(db_path: str):
    """Verifica si la base de datos BLAST existe"""
    db_files = [f"{db_path}.{ext}" for ext in ["phr", "pin", "psq"]]
    return all(Path(f).exists() for f in db_files)

def print_mode_info(mode: str):
    """Imprime información sobre el modo seleccionado"""
    if mode == "remote":
        print("🌐 MODO REMOTO:")
        print("   ✓ Usa servidores de NCBI")
        print("   ✓ No requiere instalación local")
        print("   ⚠️ Requiere conexión a internet")
        print("   ⚠️ Tiene límites de velocidad")
        print("   ⚠️ Datos se envían a NCBI")
    else:
        print("💻 MODO LOCAL:")
        print("   ✓ Ejecuta en tu computadora")
        print("   ✓ Muy rápido")
        print("   ✓ Sin límites de velocidad")
        print("   ✓ Datos permanecen locales")
        print("   ⚠️ Requiere BLAST+ instalado")
        print("   ⚠️ Requiere base de datos local")

def main():
    parser = argparse.ArgumentParser(
        description="BLAST Unificado - Remoto y Local",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python blast.py --mode remote -i secuencias.fasta -o resultado.xml
  python blast.py --mode local -i secuencias.fasta -o resultado.xml
  python blast.py --mode local -f 6 -o resultado.tsv -i secuencias.fasta
        """
    )
    
    # Argumentos principales
    parser.add_argument("--mode", required=True, choices=["remote", "local"], 
                       help="Modo de ejecución: 'remote' (NCBI) o 'local' (computadora)")
    parser.add_argument("-i", "--input", required=True, 
                       help="Archivo FASTA con secuencias")
    parser.add_argument("-o", "--output", default="blast_results.xml", 
                       help="Archivo de salida")
    
    # Argumentos de formato
    parser.add_argument("-f", "--outfmt", type=int, default=5, choices=[5,6], 
                       help="Formato de salida: 5=XML, 6=Tabular")
    parser.add_argument("-n", "--num_hits", type=int, default=10, 
                       help="Número máximo de hits (default: 10)")
    
    # Argumentos específicos para modo local
    parser.add_argument("-d", "--database", default="swissprot", 
                       help="Base de datos BLAST local (solo para modo local)")
    
    # Argumentos de información
    parser.add_argument("--info", action="store_true", 
                       help="Mostrar información detallada del modo seleccionado")
    
    args = parser.parse_args()
    
    # Mostrar información del modo
    print_mode_info(args.mode)
    print()
    
    if args.info:
        return
    
    # Verificaciones previas según el modo
    if args.mode == "local":
        print("🔍 Verificando instalación local...")
        
        if not check_blast_installation():
            print("\n❌ BLAST+ no está instalado o no está en el PATH")
            print("📋 Instrucciones de instalación:")
            print("   Windows: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download")
            print("   Linux: sudo apt-get install ncbi-blast+")
            print("   macOS: brew install blast")
            return
        
        if not check_database(args.database):
            print(f"\n❌ Base de datos '{args.database}' no encontrada")
            print("📋 Para configurar la base de datos:")
            print("   1. Descargar Swiss-Prot: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz")
            print("   2. Crear base de datos: makeblastdb -in swissprot.fasta -dbtype prot -out swissprot")
            print("   3. O ejecutar: python setup_local_blast.py")
            return
        
        print("✅ Todo listo para BLAST local")
    else:
        print("✅ Todo listo para BLAST remoto")
    
    print()
    
    # Ejecutar BLAST según el modo
    try:
        if args.mode == "remote":
            run_remote_blast(Path(args.input), Path(args.output), args.outfmt, args.num_hits)
        else:
            run_local_blast(Path(args.input), Path(args.output), args.outfmt, args.num_hits, args.database)
        
        print(f"\n🎉 ¡Completado! Resultados guardados en: {args.output}")
        
        # Mostrar estadísticas del archivo de salida
        if Path(args.output).exists():
            size_mb = Path(args.output).stat().st_size / (1024 * 1024)
            print(f"📊 Tamaño del archivo: {size_mb:.2f} MB")
            
    except KeyboardInterrupt:
        print("\n⏹️ Operación cancelada por el usuario")
    except Exception as e:
        print(f"\n❌ Error durante la ejecución: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())