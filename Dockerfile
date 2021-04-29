FROM ubuntu:latest
RUN apt-get update && apt-get install -y \
  git hmmer python3 python3-pip wget build-essential --no-install-recommends && \
  rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Set up SignalP
# Copy from Dockerfile folder, need user to download because of license
COPY signalp-5.0b.Linux.tar.gz /app
RUN tar xzvf signalp-5.0b.Linux.tar.gz && rm signalp-5.0b.Linux.tar.gz

# Set up FASTA
RUN mkdir fasta2 && \
  cd fasta2 && \
  wget https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta2/fasta2.shar.Z --no-check-certificate && \
  gunzip fasta2.shar.Z && \
  sh fasta2.shar && \
  make lfasta && \
  rm $(find . -type f ! -name "lfasta") && \ 
  cd ..

# Set up RADAR
RUN wget https://sourceforge.net/projects/repeatradar/files/radar-1.1.5.tar.gz/download \
    -O radar-1.1.5.tar.gz --no-check-certificate && \
  tar xzvf radar-1.1.5.tar.gz && \
  cd radar-1.1.5 && \
  ./configure && \
  make && \
  mv src/radar . && \
  rm -rf $(find . ! -name "radar" ! -path .) && \
  cd .. && rm radar-1.1.5.tar.gz

# Get latest fRiPPa master
RUN python3 -m pip install git+https://github.com/gamcil/frippa.git && \
 apt-get remove -y build-essential

# Add everything to PATH
ENV PATH="$PATH:/app/signalp-5.0b/bin/:/app/fasta2/:/app/radar-1.1.5/"

WORKDIR /data
