# frippa
Fungal Ribosimally synthesized and Post-translationally modified Peptide (RiPP) assayer

## Dependencies
fRiPPa requires all following tools to be available on the system $PATH:

### lfasta from the FASTA2 package

```
mkdir fasta2
cd fasta2
wget http://faculty.virginia.edu/wrpearson/fasta/fasta2/fasta2.shar.Z
gunzip fasta2.shar.Z
sh fasta2.shar
make lfasta
```

### RADAR

Download latest release from GitHub (https://github.com/AndreasHeger/radar/releases)

```
tar -xvzf radar-x.x.x.tgz
cd radar-x.x.x.tgz
./configure
make; make install
```

### SignalP 5.0

Obtain executable from: https://services.healthtech.dtu.dk/software.php

### hmmsearch from HMMER

Obtain from: http://hmmer.org/
