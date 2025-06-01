#USO: prende un fastq e seleziona le reads che hanno contenuto di GC tra   
import os 
from Bio import SeqIO

path = os.getcwd()
fastq_file = input("Inserire il nome del file fastq da analizzare: ")
lower_GC = float(input("inserisci il minimo valore % di contenuto di GC delle reads che vuoi copiare nel file di output: "))
higher_GC = float(input("inserisci il massimo valore % di contenuto di GC delle reads che vuoi copiare nel file di output: "))
file_format = input("In che formato desideri l'output? fasta/fastq")
print("le seguenti sequenze hanno la percentuale di GC selezionata:")
for record in SeqIO.parse(fastq_file, "fastq"):
#SeqInputOutput.parse legge i dati del fastq (con "fastq" specifichi in che formato è il file con le sequenze)
	#print(record)
	a = record.seq
	split_a = list(record.seq)
	countA = 0
	countT = 0
	countG = 0
	countC = 0
	for nucleotide in split_a:
		if nucleotide == "A":
			countA = countA +1
		elif nucleotide == "T":
			countT = countT + 1
		elif nucleotide == "G":
			countG = countG + 1 
		elif nucleotide == "C":
			countC = countC + 1
	GC_percentage = ((countC + countG) / (countA + countT + countC + countG))*100	
	if GC_percentage  >= lower_GC and GC_percentage <= higher_GC:
		print(record.seq)
		print(GC_percentage)
		if file_format == 'fastq':
			with open(str(lower_GC) + "_" + str(higher_GC) +"%_GC_" + fastq_file +".fastq", "a") as output_handle:
				SeqIO.write(record, output_handle, "fastq")
		elif file_format == 'fasta':
			with open(str(lower_GC) + "_" + str(higher_GC) +"%_GC_" + fastq_file +".fasta", "a") as output_handle:
				SeqIO.write(record, output_handle, "fasta")
	#print(countG)
	#print(countC)			
	#print(split_a)
