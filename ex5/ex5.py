#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def gc_content(seq):
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq)


def tm_wallace(seq):
    A = seq.count("A") + seq.count("a")
    T = seq.count("T") + seq.count("t")
    G = seq.count("G") + seq.count("g")
    C = seq.count("C") + seq.count("c")
    return 2*(A+T) + 4*(G+C)


def valid_primer(seq, cfg):
    L = len(seq)
    if not (cfg["primer_min_len"] <= L <= cfg["primer_max_len"]):
        return False

    gc = gc_content(seq)
    if not (cfg["gc_min"] <= gc <= cfg["gc_max"]):
        return False

    if cfg["avoid_terminal_gc"]:
        if seq[-1] in "GCgc":
            return False

    tm = tm_wallace(seq)
    if tm > cfg["max_tm"]:
        return False

    return True


def design_primers(transcript, cfg):
    primers = []

    for i in range(len(transcript)):
        for L in range(cfg["primer_min_len"], cfg["primer_max_len"] + 1):
            if i + L > len(transcript):
                continue

            candidate = transcript[i:i+L]

            if valid_primer(candidate, cfg):
                primers.append(candidate)
                if len(primers) == cfg["num_primers"]:
                    return primers

    return primers


def main():
    parser = argparse.ArgumentParser(description="Diseño automático de primers.")
    parser.add_argument("-i", "--input", required=True, help="Transcripto en formato FASTA")
    parser.add_argument("-c", "--config", required=True, help="Archivo JSON con parámetros")
    parser.add_argument("-o", "--output", required=True, help="Archivo de salida con primers diseñados")
    args = parser.parse_args()

    # Leer FASTA (primera secuencia)
    seq = ""
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq += line.strip()

    # Configuración
    with open(args.config) as f:
        cfg = json.load(f)

    primers = design_primers(seq, cfg)

    if not primers:
        print("⚠ No se encontraron primers que cumplan los criterios.")
        return

    with open(args.output, "w") as out:
        for i, p in enumerate(primers, start=1):
            out.write(f">primer_{i}\n{p}\n")

    print(f"✔ Diseño completado. Se generaron {len(primers)} primers.")
    print(f"→ Archivo generado: {args.output}")


if __name__ == "__main__":
    main()
