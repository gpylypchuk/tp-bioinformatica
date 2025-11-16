# 🧬 TP Bioinformática

Este proyecto contiene los programas y scripts desarrollados para el trabajo práctico de **Bioinformática**.  
Incluye ejercicios de análisis de secuencias genómicas/proteicas, BLAST, MSA, búsqueda de motivos y diseño de primers.

---

## ⚙️ Configuración del entorno

### Crear entorno virtual
```bash
python3 -m venv venv
```

### Activar entorno

- **Mac/Linux**
  ```bash
  source venv/bin/activate
  ```

- **Windows**
  ```bash
  venv\Scripts\activate.bat
  ```

### Instalar dependencias
```bash
pip install -r requirements.txt
```

### Entorno BioEMBOSS (solo para Ej. 4)
Si usás EMBOSS dentro de Conda:
```bash
conda create -n bioemboss emboss
conda activate bioemboss
```

---

# ▶️ Ejecución de los programas

---

## 📌 Ejercicio 1 – Traducción en los 6 frames

Convierte un archivo GenBank del gen **HBB** a aminoácidos en los seis marcos de lectura.

```bash
python3 ex1.py -i sequence.gb -o query_HBB_147aa.fasta
```

**Salida:** archivo FASTA con 16 ORFs (frames +1, +2, +3, -1, -2, -3).

---

## 📌 Ejercicio 2 – BLAST remoto y local

### BLAST remoto (NCBI)
```bash
python3 blast.py --mode remote -i seq.fasta -o blast_remote.xml
```

### BLAST local
```bash
python3 blast.py --mode local -i seq.fasta -o blast_local.xml -d ruta/a/db
```

**Salida:** archivos `.xml` con formato outfmt=5.

---

## 📌 Ejercicio 3 – Top 10 BLAST + Alineamiento Múltiple (MSA)

```bash
export NCBI_EMAIL="tu_email@dominio.com"

python3 blast_top10_msa.py \
  --query query.faa \
  --remote \
  --outdir out \
  --msa-tool muscle
```

O utilizando un archivo XML ya existente:

```bash
python3 blast_top10_msa.py \
  --xml blast.xml \
  --query query.faa \
  --outdir out \
  --msa-tool clustalo
```

---

# 📌 Ejercicio 4 – Búsqueda de motivos PROSITE usando EMBOSS

### 1️⃣ Preparar PROSITE
```bash
prosextract -prositedir ~/prosite
```

### 2️⃣ Ejecutar el script

```bash
python3 ex4.py \
  -i ex1_all_orfs.faa \
  -o ex4_prosite_results.txt \
  --prosite-dir ~/prosite
```

**Salida:** resultados por ORF con `HitCount` y coincidencias PROSITE.

---

# 📌 Ejercicio 5 – Diseño automático de primers

### 1️⃣ Archivo `config.json`
```json
{
  "primer_min_len": 18,
  "primer_max_len": 24,
  "gc_min": 0.50,
  "gc_max": 0.60,
  "max_tm": 67,
  "avoid_terminal_gc": true,
  "num_primers": 5
}
```

### 2️⃣ Convertir GenBank → FASTA

```bash
seqret -sequence sequence.gb -sformat1 genbank -osformat2 fasta -outseq transcripto.fa
```

### 3️⃣ Ejecutar el script

```bash
python3 ex5.py \
  -i transcripto.fa \
  -c config.json \
  -o primers_generados.txt
```

---

# 👥 Contribuidores

| [![agmiguel-k](https://github.com/agmiguel-k.png?size=60)](https://github.com/agmiguel-k) | [![gpylypchuk](https://github.com/gpylypchuk.png?size=60)](https://github.com/gpylypchuk) | [![thiagosayegh](https://github.com/thiagosayegh.png?size=60)](https://github.com/thiagosayegh) | [![s-gruss](https://github.com/s-gruss.png?size=60)](https://github.com/s-gruss) |
|:---:|:---:|:---:|:---:|
| [agmiguel-k](https://github.com/agmiguel-k) | [gpylypchuk](https://github.com/gpylypchuk) | [thiagosayegh](https://github.com/thiagosayegh) | [s-gruss](https://github.com/s-gruss) |
