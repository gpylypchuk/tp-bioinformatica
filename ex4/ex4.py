import argparse
import subprocess
from pathlib import Path
import sys
import tempfile
import shutil


def run(cmd):
    """Ejecuta un comando externo y corta si hay error."""
    print("[RUN]", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"Error ejecutando: {' '.join(cmd)}")


def parse_fasta(path: Path):
    """Devuelve una lista de (header, sequence) a partir de un FASTA."""
    records = []
    header = None
    seq_lines = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue  # salteo líneas vacías
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


def main():
    parser = argparse.ArgumentParser(
        description=("Ejercicio 4 - EMBOSS + PROSITE sobre ORFs")
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Archivo FASTA de aminoácidos con TODOS los ORFs (ex1_all_orfs.faa)",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Archivo de salida con los dominios encontrados (se sobreescribe)",
    )
    parser.add_argument(
        "--prosite-dir",
        default=str(Path.home() / "prosite"),
        help="Directorio con prosite.dat/prosite.doc/prosite.lines (por defecto: ~/prosite)",
    )
    args = parser.parse_args()

    input_faa = Path(args.input).resolve()
    prosite_dir = Path(args.prosite_dir).resolve()
    output_txt = Path(args.output).resolve()

    if not input_faa.exists():
        sys.exit(f"No encuentro el archivo de entrada: {input_faa}")

    # Chequeo PROSITE
    if not (prosite_dir / "prosite.lines").exists():
        print("No encuentro prosite.lines, corro prosextract...")
        run(["prosextract", "-prositedir", str(prosite_dir), "-auto"])

    # Parseo el FASTA a nivel Python
    records = parse_fasta(input_faa)
    print(f"Leí {len(records)} ORFs del archivo {input_faa.name}")

    # Creo/limpio el archivo de salida
    if output_txt.exists():
        output_txt.unlink()

    tmpdir = Path(tempfile.mkdtemp(prefix="ej4_prosite_"))
    print(f"Usando directorio temporal: {tmpdir}")

    try:
        for idx, (header, seq) in enumerate(records, start=1):
            orf_id = f"orf_{idx:02d}"
            print(f"\n=== Analizando {orf_id}: {header} (len={len(seq)}) ===")

            # 1) Escribo un FASTA temporal para este ORF
            tmp_faa = tmpdir / f"{orf_id}.faa"
            with tmp_faa.open("w") as f:
                f.write(f">{header}\n")
                # Parto la secuencia en líneas de 60 aa (opcional)
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")

            # 2) Archivo de salida temporal de EMBOSS
            tmp_out = tmpdir / f"{orf_id}_prosite.txt"

            # 3) Corro patmatmotifs SOLO para este ORF
            run([
                "patmatmotifs",
                "-sequence", str(tmp_faa),
                "-outfile", str(tmp_out),
                "-full", "Y",
                "-prune", "N",
            ])

            # 4) Añado el resultado al archivo general
            with output_txt.open("a") as fout, tmp_out.open() as fin:
                fout.write(f"\n################ ORF {idx}: {header} ################\n")
                shutil.copyfileobj(fin, fout)

    finally:
        # Borro los temporales
        shutil.rmtree(tmpdir)

    print("\nListo!")
    print(f"  Entrada : {input_faa}")
    print(f"  Salida  : {output_txt}")


if __name__ == "__main__":
    main()
