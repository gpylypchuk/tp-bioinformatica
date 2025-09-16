# tp-bioinformatica

`python3 -m venv venv`

## Mac/linux

`source venv/bin/activate`

## Windows

`venv\Scripts\activate.bat`

## Install dependencies

`pip install -r requirements.txt`

# Run the programs

## Ejercicio 1

`python3 ex1.py -i HBB_NM000518_5.gbk -o HBB_AA_6frames.fasta`

## Ejercicio 2

`python3 ex2.py -i seq_1.fasta -o blast_remote.xml`

## Ejercicio 3

1.  Tené una consulta válida (AA). Ej: la ORF correcta de HBB (~147 aa)
    Supongamos que está en query.faa

2.  BLAST remoto + top10 + MSA con MUSCLE

`export NCBI_EMAIL="tu_email@dominio.com" # requerido por Entrez`

`python blast_top10_msa.py --query query.faa --remote --outdir out --msa-tool muscle`

3. Si ya tenés un BLAST XML (outfmt 5), reusalo:

`python blast_top10_msa.py --xml blast.xml --query query.faa --outdir out --msa-tool clustalo`
