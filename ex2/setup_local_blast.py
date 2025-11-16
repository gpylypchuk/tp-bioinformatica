#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script de configuración para BLAST local
# Este script ayuda a descargar e instalar los componentes necesarios

import os
import sys
import subprocess
import urllib.request
import gzip
import shutil
from pathlib import Path

# Configurar encoding para Windows
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    
    def get_short_path(path):
        """Obtiene la ruta corta (8.3) de Windows para evitar problemas con espacios y paréntesis"""
      
        try:
                import ctypes
                from ctypes import wintypes
                _GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW
                _GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
                _GetShortPathNameW.restype = wintypes.DWORD
                
                path_str = str(path)
                buf = ctypes.create_unicode_buffer(260)
                result = _GetShortPathNameW(path_str, buf, 260)
                if result:
                    return buf.value
        except:
                pass
            # Si todo falla, devolver la ruta original
        return str(path)
else:
    def get_short_path(path):
        """En sistemas no-Windows, devolver la ruta normal"""
        return str(path)

def download_file(url: str, filename: str):
    """Descarga un archivo desde URL"""
    print(f"Descargando {filename}...")
    try:
        urllib.request.urlretrieve(url, filename)
        print(f"✓ Descargado: {filename}")
        return True
    except Exception as e:
        print(f"✗ Error descargando {filename}: {e}")
        return False

def extract_gz(gz_file: str, output_file: str):
    """Extrae archivo .gz"""
    print(f"Extrayendo {gz_file}...")
    try:
        with gzip.open(gz_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"✓ Extraído: {output_file}")
        return True
    except Exception as e:
        print(f"✗ Error extrayendo {gz_file}: {e}")
        return False

def run_makeblastdb(fasta_file: str, db_name: str):
    """Ejecuta makeblastdb y genera también índices extendidos (.pdb, .pjs, .ptf, .pto, .pot)"""
    print(f"Creando base de datos BLAST: {db_name}")
    try:
        current_dir = Path.cwd().resolve()
        fasta_path = Path(fasta_file).resolve()

        if not fasta_path.exists():
            print(f"✗ Error: Archivo FASTA no encontrado: {fasta_path}")
            return False

        db_path = current_dir / db_name

        if sys.platform == 'win32':
            fasta_str = get_short_path(fasta_path)
            db_str = get_short_path(db_path)
            print(f"Usando rutas cortas de Windows para evitar problemas con espacios/paréntesis")
        else:
            fasta_str = str(fasta_path)
            db_str = str(db_path)

        # 1️⃣ Crear base de datos principal
        cmd = [
            "makeblastdb",
            "-in", fasta_str,
            "-dbtype", "prot",
            "-parse_seqids",  # permite crear índices adicionales
            "-out", db_str
        ]

        print(f"Ejecutando: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(current_dir), encoding='utf-8', errors='replace')

        if result.returncode != 0:
            print(f"✗ Error creando base de datos:")
            print(result.stderr or result.stdout)
            return False

        print(f"✓ Base de datos principal creada")

        # 2️⃣ Crear índices extendidos (.pdb, .pjs, .pot, .ptf, .pto)
        alias_cmd = [
            "blastdb_aliastool",
            "-db", db_str,
            "-dbtype", "prot",
            "-out", db_str,
            "-title", f"{db_name} protein index"
        ]

        print(f"Creando índices extendidos...")
        alias_result = subprocess.run(alias_cmd, capture_output=True, text=True, cwd=str(current_dir), encoding='utf-8', errors='replace')

        if alias_result.returncode == 0:
            print(f"✓ Índices extendidos (.pdb, .pjs, .pot, .ptf, .pto) creados correctamente")
        else:
            print(f"⚠ No se pudieron crear los índices extendidos")
            print(alias_result.stderr or alias_result.stdout)

        return True

    except Exception as e:
        print(f"✗ Error inesperado: {e}")
        return False


def check_blast_installation():
    """Verifica si BLAST+ está instalado"""
    try:
        result = subprocess.run(["blastp", "-version"], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ BLAST+ está instalado")
            return True
    except FileNotFoundError:
        pass
    
    print("✗ BLAST+ no está instalado")
    return False

def check_database_exists(db_name: str):
    """Verifica si la base de datos BLAST ya existe"""
    db_files = [f"{db_name}.{ext}" for ext in ["phr", "pin", "psq"]]
    all_exist = all(Path(f).exists() for f in db_files)
    
    if all_exist:
        print(f"✓ Base de datos '{db_name}' ya existe")
        for db_file in db_files:
            if Path(db_file).exists():
                size_mb = Path(db_file).stat().st_size / (1024 * 1024)
                print(f"  - {Path(db_file).name} ({size_mb:.1f} MB)")
        return True
    return False

def main():
    print("=== Configuración de BLAST Local ===\n")
    
    # Verificar BLAST+
    if not check_blast_installation():
        print("\n❌ BLAST+ no está instalado. Instalación requerida:")
        print("\n📋 Instrucciones de instalación:")
        print("   Windows:")
        print("   1. Descargar desde: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download")
        print("   2. Extraer y agregar al PATH del sistema")
        print("   3. Reiniciar terminal/IDE")
        print("\n   Linux:")
        print("   sudo apt-get update && sudo apt-get install ncbi-blast+")
        print("\n   macOS:")
        print("   brew install blast")
        print("\n   Después de instalar, ejecuta este script nuevamente.")
        return
    
    # Verificar si la base de datos ya existe
    db_name = "swissprot"
    print(f"\n🔍 Verificando si la base de datos '{db_name}' ya existe...")
    
    if check_database_exists(db_name):
        print(f"\n✅ La base de datos ya está configurada y lista para usar!")
        print(f"\n🚀 Puedes ejecutar BLAST local ahora:")
        print(f"   python blast.py --mode local -i tu_archivo.fasta -o resultado.xml")
        print(f"\n💡 Si quieres recrear la base de datos, borra los archivos:")
        print(f"   - {db_name}.phr")
        print(f"   - {db_name}.pin")
        print(f"   - {db_name}.psq")
        return
    
    # Si no existe, proceder con la descarga y creación
    swissprot_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz"
    swissprot_gz = "swissprot.gz"
    swissprot_fasta = "swissprot.fasta"
    
    print(f"\n📥 Descargando base de datos Swiss-Prot...")
    if not download_file(swissprot_url, swissprot_gz):
        print("❌ No se pudo descargar Swiss-Prot")
        return
    
    # Extraer archivo
    if not extract_gz(swissprot_gz, swissprot_fasta):
        print("❌ No se pudo extraer Swiss-Prot")
        return
    
    # Crear base de datos BLAST
    if not run_makeblastdb(swissprot_fasta, db_name):
        print("❌ No se pudo crear la base de datos")
        return
    
    # Limpiar archivos temporales
    print("\n🧹 Limpiando archivos temporales...")
    for temp_file in [swissprot_gz, swissprot_fasta]:
        if Path(temp_file).exists():
            Path(temp_file).unlink()
            print(f"✓ Eliminado: {temp_file}")
    
    print(f"\n✅ Configuración completada!")
    print(f"   Base de datos: {db_name}")
    print(f"   Archivos creados: {db_name}.phr, {db_name}.pin, {db_name}.psq")
    print(f"\n🚀 Ahora puedes ejecutar:")
    print(f"   python blast.py --mode local -i tu_archivo.fasta -o resultado.xml")

if __name__ == "__main__":
    main()
