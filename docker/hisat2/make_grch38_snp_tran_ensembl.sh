#!/bin/sh

#
# Downloads sequence for the GRCh38 release N version of H. sapiens (human) from
# Ensembl.
#
# Note that Ensembl's GRCh38 build has three categories of compressed fasta
# files:
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
# INPUT:
# ENSEMBL_RELEASE: release version , ex 90
# DBSNP_RELEASE : release version of dbsnp, ex 150

ENSEMBL_RELEASE=$1
ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
ENSEMBL_GRCh38_GTF_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens
GTF_FILE=Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf

DBSNP_RELEASE=$2
SNP_FILE=snp${DBSNP_RELEASE}Common.txt
UCSC_COMMON_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/${SNP_FILE}
HISAT2_DIR='/opt/tools/hisat2-2.1.0'
get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

HISAT2_BUILD_EXE=${HISAT2_DIR}/hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

HISAT2_SNP_SCRIPT=${HISAT2_DIR}/hisat2_extract_snps_haplotypes_UCSC.py
if [ ! -x "$HISAT2_SNP_SCRIPT" ] ; then
	if ! which hisat2_extract_snps_UCSC.py ; then
		echo "Could not find hisat2_extract_snps_haplotypes_UCSC.py in current directory or in PATH"
		exit 1
	else
		HISAT2_SNP_SCRIPT=`which hisat2_extract_snps_haplotypes_UCSC.py`
	fi
fi

HISAT2_SS_SCRIPT=${HISAT2_DIR}/hisat2_extract_splice_sites.py
if [ ! -x "$HISAT2_SS_SCRIPT" ] ; then
	if ! which hisat2_extract_splice_sites.py ; then
		echo "Could not find hisat2_extract_splice_sites.py in current directory or in PATH"
		exit 1
	else
		HISAT2_SS_SCRIPT=`which hisat2_extract_splice_sites.py`
	fi
fi

HISAT2_EXON_SCRIPT=${HISAT2_DIR}/hisat2_extract_exons.py
if [ ! -x "$HISAT2_EXON_SCRIPT" ] ; then
	if ! which hisat2_extract_exons.py ; then
		echo "Could not find hisat2_extract_exons.py in current directory or in PATH"
		exit 1
	else
		HISAT2_EXON_SCRIPT=`which hisat2_extract_exons.py`
	fi
fi

rm -f genome.fa
F=Homo_sapiens.GRCh38.dna.primary_assembly.fa
if [ ! -f $F ] ; then
	get ${ENSEMBL_GRCh38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi


if [ ! -f $SNP_FILE ] ; then
       get ${UCSC_COMMON_SNP}.gz || (echo "Error getting ${UCSC_COMMON_SNP}" && exit 1)
       gunzip ${SNP_FILE}.gz || (echo "Error unzipping ${SNP_FILE}" && exit 1)
       awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' ${SNP_FILE} > ${SNP_FILE}.tmp
       mv ${SNP_FILE}.tmp ${SNP_FILE}
       ${HISAT2_SNP_SCRIPT} genome.fa ${SNP_FILE} genome
fi

if [ ! -f $GTF_FILE ] ; then
       get ${ENSEMBL_GRCh38_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
       gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
       ${HISAT2_SS_SCRIPT} ${GTF_FILE} > genome.ss
       ${HISAT2_EXON_SCRIPT} ${GTF_FILE} > genome.exon
fi

CMD="${HISAT2_BUILD_EXE} -p 8 genome.fa --snp genome.snp --haplotype genome.haplotype --ss genome.ss --exon genome.exon genome_snp_tran"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
