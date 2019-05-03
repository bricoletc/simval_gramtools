data_DIR=$1


mkdir -p ${slice_DIR}
contig_names=("Pf3D7_01_v3" "Pf3D7_02_v3")
start="105000"
end="130000"

ref=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.fa(sta)?")
vcf=$(find ${data_DIR} -maxdepth 1 -regextype posix-egrep -regex ".*.vcf(.gz)?")

sliced_ref=${slice_DIR}/slice.fa; touch "${sliced_ref}"; truncate -s 0 "${sliced_ref}"
sliced_vcf=${slice_DIR}/slice.vcf; touch "${sliced_vcf}"; truncate -s 0 "${sliced_vcf}"
rebased_sliced_vcf=${slice_DIR}/rebased_slice.vcf

# Need indexes of the fasta and vcf files for below.
num=0
for contig_name in ${contig_names[@]}; do
	echo $contig_name
	num=$((num+1))
	samtools faidx ${ref} ${contig_name}:${start}-${end} >> ${sliced_ref}
	if [ $num -gt 1 ]; then
		bcftools filter ${vcf} ${contig_name}:${start}-${end} | grep -v "^#" >> ${sliced_vcf}
	else
		bcftools filter ${vcf} ${contig_name}:${start}-${end} >> ${sliced_vcf}
	fi
done

# Rebase the vcf coords relative to the slice start offset
python3 ${script_DIR}/rebase_vcf.py ${sliced_vcf} ${start} ${rebased_sliced_vcf} # Second arg is offset to remove from vcf records.

# Make a space in fasta headers 
sed -i '/^>/ s/:/ :/' ${sliced_ref}
