data_DIR=$1
samples_DIR=$2
num_samples=$3

usage(){
	echo "##############################"
	echo "usage: $0 data_DIR samples_DIR num_samples"
	echo "data_DIR must contain a fasta, and a vcf describing variation against it"
	echo "samples_DIR is where the samples will be placed, under a directory called \"generated_samples\""
	echo "##############################"
	exit 0
}

if [ -z "${data_DIR}" -o -z "${num_samples}" -o -z "${samples_DIR}" ]; then
	usage
fi

ref=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.fa(sta)?")
vcf=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.vcf(.gz)?")


if [ -z "${ref}" -o -z "${vcf}" ] ; then
	echo "Missing either fasta ref or vcf file in ${data_DIR} provided"
	usage
fi


script_DIR=$(dirname $0)

samples_DIR="${samples_DIR}/generated_samples"
mkdir -p ${samples_DIR}

sample_num=0
while [ $sample_num -lt ${num_samples} ]; do
	sample_num=$[$sample_num + 1]
	echo -e "Building sample number ${sample_num}"
	python3 ${script_DIR}/generate_samples/one_sample.py ${samples_DIR} ${sample_num} ${data_DIR}
done

