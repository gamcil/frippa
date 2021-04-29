# frippa
Fungal Ribosimally synthesized and Post-translationally modified Peptide (RiPP) assayer

## Installation
The quickest way to start using fRiPPa is by using Docker. Note that due to licensing,
you have to grab your own copy of SignalP 5.0b first ([here](https://services.healthtech.dtu.dk/software.php)).

First, download the [fRiPPa Dockerfile](https://raw.githubusercontent.com/gamcil/frippa/master/Dockerfile)
and place it in the same directory as your downloaded copy of SignalP. You should have
the following files:

```
Dockerfile
signalp-5.0b.Linux.tar.gz
```

Then, you can build a Docker image by:

```
docker build -t frippa -f Dockerfile .
```

Once you have built the image, the easiest way to use it is by using the [fRiPPA wrapper
script](https://raw.githubusercontent.com/gamcil/frippa/master/frippa-docker). For
example:

```sh
wget https://raw.githubusercontent.com/gamcil/frippa/master/frippa-docker
chmod +x frippa-docker
frippa-docker genome.gbk
```

This will download the wrapper script, ``frippa-docker``, make it executable, and run
it on ``genome.gbk``. You can put this script somewhere on your system ``$PATH`` to
access it from anywhere.

The wrapper script will mount the present working directory (``$PWD``; the directory
you're in when you launched ``frippa-docker``) as the working directory, allowing you to
pass files between the host system and the docker container. Note: files must be
somewhere **inside** the ``$PWD``; using e.g. ``../path/outside/pwd`` will not work.

Under the hood, this is the equivalent of:

```sh
WORKDIR="${PWD}"
MOUNT="type=bind,source=${WORKDIR},target=${WORKDIR}"

docker run --rm "${IT[@]}" \
  --workdir "${WORKDIR}" \
  --mount "${MOUNT}" \
  frippa \
  frippa genome.gbk
```

## Usage
```
usage: frippa [-h] [-o OUTPUT] [-t THREADS] [-e EVALUE] [-n NEIGHBOURS]
              [-v OVERLAP] [-ra] [-cp CUTSITE_PCT] [-mr MIN_REPEATS]
              [-min_rl MIN_REPEAT_LENGTH] [-max_rl MAX_REPEAT_LENGTH]
              [-rs REPEAT_SIMILARITY] [-zc Z_SCORE_CUTOFF]
              genbank

positional arguments:
  genbank               Input GenBank file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path to write results to
  -t THREADS, --threads THREADS
                        How many threads to use when running hmmsearch and
                        SignalP
  -e EVALUE, --evalue EVALUE
                        E-value cutoff to use when filtering hmmsearch results

Cluster prediction settings:
  -n NEIGHBOURS, --neighbours NEIGHBOURS
                        Number of base pair either side of a DUF3328 hit to
                        grab proteins
  -v OVERLAP, --overlap OVERLAP
                        Maximum percentage of overlapping proteins before
                        clusters are merged
  -ra, --report_all     Report all clusters with DUF3328 hits, even if no
                        repeat predictions

Repeat filtering:
  -cp CUTSITE_PCT, --cutsite_pct CUTSITE_PCT
                        Threshold for total cut sites in individual repeats
  -mr MIN_REPEATS, --min_repeats MIN_REPEATS
                        Minimum number of detected repeats
  -min_rl MIN_REPEAT_LENGTH, --min_repeat_length MIN_REPEAT_LENGTH
                        Minimum length of a repeat sequence.
  -max_rl MAX_REPEAT_LENGTH, --max_repeat_length MAX_REPEAT_LENGTH
                        Maximum length of a repeat sequence. Long sequences
                        tend to be red herrings
  -rs REPEAT_SIMILARITY, --repeat_similarity REPEAT_SIMILARITY
                        Average (hamming) similarity threshold for repeat
                        sequences
  -zc Z_SCORE_CUTOFF, --z_score_cutoff Z_SCORE_CUTOFF
                        Minimum z-score for a repeat sequence
```

Note that you will have to provid

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
