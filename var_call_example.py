
from Bio import SeqIO
from shutil import copyfile
import os, argparse, sys, re, subprocess

def main():

	parser = argparse.ArgumentParser(description="The Program will call variants given Illumina Data and a reference")

	parser.add_argument('-1', action="store", help="The path to a set of reads", dest='r1')
	parser.add_argument('-2', action="store", help="The path to a second set of reads if one exists.", dest='r2')
	parser.add_argument('-r', action="store", help="The path to your reference genome", dest='ref')
	parser.add_argument('-o', action="store", help="The path to where you would like an output file generated", dest="out")

	args = parser.parse_args()

	r1 = args.r1
	r2 = args.r2
	ref = args.ref
	out = args.out

	print "\n"
	if ref == None:
		print "No reference file found, try -h for help!\n"
		exit(1)

	IsRef = os.path.isfile(ref)

	if IsRef in [False]:
		print "The reference file was not found\n"
		exit(1)


        if r1 and r2 == None:

                print "No reads were found!\n"
                exit(1)

	IsR1=False
	IsR2=False

	if r1 is not None:
		IsR1 = os.path.isfile(r1.encode('utf-8'))

	if r2 is not None:
		IsR2 = os.path.isfile(r2.encode('utf-8'))

	if IsR1 == False and IsR2 == False:
		print "The path to your read(s) does not exist!"
		exit(1)

	if out in [None]:
		print "No output path specified\n"
		exit(1)

	OutExist = os.path.isdir(out)
	if OutExist == True:
		print "Could not create the output directory, is already exists!\n"
		exit(1)

	try:
		os.mkdir(out)

	except:
		print "The output path could not be created\n"
		exit(1)

	print "Reference Found\n"
	print "Read(s) Found\n"
	print "Output Folder Sucessfully Created\n"

	ReferenceCreator(r1, r2, out, ref, IsR1, IsR2)

def ReferenceCreator(r1, r2, out, ref, IsR1, IsR2):

	FileFormat_fasta = ref.find('fasta')
	FileFormat_fa = ref.find ('fa')

	if FileFormat_fasta == -1 and FileFormat_fa ==-1:
		print "Attempting to convert reference to a fasta format\n"
		FastaExtension = ref.find(".fastq")
		FaExtension = ref.find(".fa")
		#FE = ref[FileExtension+1:]
		###print FE
		if FastaExtension ==-1:
			with open(ref, "rU") as input_handle, open (out+"/reference.fastq", 'w') as output_handle:
				sequences = SeqIO.parse(input_handle, "genbank")

				SeqIO.write(sequences, output_handle, "fastq")
				output_handle.close()

	if FileFormat_fasta !=-1 or FileFormat_fa !=-1:
		print "A copy of your file has been brought to your output directory for manipulation\n"
		copyfile(ref,out+"/reference.fastq")

	reference = out+"/reference.fasta"

	print "Creating a the reference sequence index\n"
	subprocess.check_call("samtools faidx "+reference, shell=True)

	print "Creating a reference sequence dictionary\n"
	subprocess.check_call("picard-tools CreateSequenceDictionary REFERENCE="+reference+" OUTPUT="+out+"/reference.dict", shell=True)

	BWA(r1, r2, out, ref, IsR1, IsR2,reference)

def BWA(r1, r2, out, ref, IsR1, IsR2, reference):

	print "Giving BWA the Index\n"

	bwa_index= out+"/bwa_index"
	subprocess.check_call("bwa index -p "+ bwa_index + " -a is  "+ reference, shell=True)

	print "Mapping Reads to Reference with BWA mem\n"

	if IsR1 == True and IsR2 == True:

		subprocess.check_call("bwa mem -t 4 " +bwa_index+" "+r1+" "+r2+" > "+out+"/aln_reads.sam", shell=True)

	if IsR1 == True and IsR2 == False:

		subprocess.check_call("bwa mem -t 4 " +bwa_index+" "+r1+" > "+out+"/aln_reads.sam", shell=True)

	if IsR1 == False and IsR2 == True:

		subprocess.check_call("bwa mem -t 4 " +bwa_index+" "+r2+" > "+out+"/aln_reads.sam", shell=True)

	print "Converting sam output to bam output\n"

	subprocess.check_call("samtools view -S "+out+"/aln_reads.sam -b -o "+out+"/aln_reads.bam -@ 8", shell=True)

	print "Sorting bam output\n"

	subprocess.check_call("samtools sort -@ 8 "+out+"/aln_reads.bam -o "+out+"/aln_sorted_reads.bam", shell=True)

	print "Indexing the bam file\n"

	subprocess.check_call("samtools index "+out+"/aln_sorted_reads.bam", shell=True)

	print "Creating a pileup using samtools mpileup, which then pipes stdout to bcftools to call variants\n"

	subprocess.check_call("samtools mpileup -ugf "+reference+" "+out+"/aln_sorted_reads.bam | bcftools call -vmO z -o "+out+"/raw_snps.vcf.gz", shell=True)

	subprocess.check_call("samtools mpileup -f "+reference+" -s "+out+"/aln_sorted_reads.bam > "+out+"/mpileup.tsv", shell=True)

	print "Filtering variant calls.  They will appear in filtered_snps.vcf\n"

	subprocess.check_call("bcftools filter -O z -o "+out+"/filtered_snps.vcf.gz -s LOWQUAL -i'%QUAL>20' " +out+"/raw_snps.vcf.gz", shell=True)

	print "Opening your output with gunzip\n"

	subprocess.check_call("gunzip "+out+"/raw_snps.vcf.gz", shell=True)
	subprocess.check_call("gunzip "+out+"/filtered_snps.vcf", shell=True)

	print "Program Complete!  Thanks for using!\n"

if __name__ == "__main__":
    main()
