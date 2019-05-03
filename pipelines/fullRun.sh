set -e
data_DIR=$1
KMER_SIZE=$2

data_DIR=$(realpath ${data_DIR})

ref=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.fa(sta)?")
vcf=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.vcf")
reads=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.(fastq|bam|fq)")

gram_DIR=${data_DIR}/gram_k${KMER_SIZE}

existing_build=${gram_DIR}/kmers_stats

gramtools_loc="/home/brice/Desktop/git_repos/gramtools"
source "${gramtools_loc}"/venv/bin/activate

run_DIR=${data_DIR}/run_k${KMER_SIZE}

# For reference
quasimap_DIR=${run_DIR}/quasimap_outputs
infer_DIR=${run_DIR}/infer_outputs
discover_DIR=${run_DIR}/discover_outputs


## Build
gramtools build --gram-dir ${gram_DIR} --reference ${ref} --vcf ${vcf} \
	--kmer-size ${KMER_SIZE} --all-kmers


## Quasimap

gramtools quasimap --gram-directory ${gram_DIR} --run-dir ${run_DIR} --reads ${reads} 


## Infer
gramtools infer --run-dir ${run_DIR} 


## Discover
gramtools discover --run-dir ${run_DIR} 


