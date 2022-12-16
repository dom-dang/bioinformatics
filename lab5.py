#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import sys

#error checker
def lineChecker (line):
    validChar = ["A", "T", "C", "G", "a", "t", "c", "g"," ", "/t"]
    flag = 0
    for letter in line: 
        for char in validChar: 
            if letter == char: 
                flag +=1
    if (flag == getLength(line)):
        return (True)
    else:
        return (False)

def getLength(word):
    num = 0
    for thing in word:
         num +=1
    return (num)

def getCount(word, char):
    num = 0
    for letter in word: 
        if (letter == char):
            num +=1
    return (num)

def getUpper (word):
    result = ""
    for char in word:
        if (ord(char) >= 97):
            result += chr(ord(char) -32)
        else:
            result += char
    return (result)



def calcGC (header, sequence):
    lenSeq = getLength (sequence)
    sequence = getUpper(sequence)
    if (lineChecker(sequence)):
        gCount = getCount (sequence, "G")
        cCount = getCount(sequence, "C")
        GCpercent = ((float(gCount)+float(cCount))/lenSeq) * 100
        return (header+"\t"+str(GCpercent) + "%\n")
    else: 
        return (header + "\t" + "ERROR\n") 

# use calcGC above for this function
def gInput(inFile, outFile):
    lineCount = 1
    with open(inFile) as input, open(outFile, "w") as output: 
        output.write("ID\tGC%\n")
        for line in input:
            line = line.strip()
            line = re.sub("\t", "", line)
            if (lineCount % 2 == 1):
                header = line[1:]
            else: 
                seq = re.sub(" ", "", line)
            if (lineCount % 2 == 0):
                output.write(calcGC(header, seq))
                header = ""
                seq = ""     
            lineCount +=1



def revComp (header, sequence):
    complement = ""
    sequence = getUpper(sequence)
    if (lineChecker(sequence)):     
        for letter in sequence: 
            if (letter == "A"):
                complement += "T"
            elif (letter == "T"):
                complement += "A"
            elif (letter == "C"):
                complement += "G"
            elif (letter == "G"):
                complement += "C"
        reverse = complement[::-1]
        return (header + "\n" + reverse + "\n")
    else: 
        return (header + "\n" + "ERROR\n")

def rInput(inFile, outFile):
    lineCount = 1
    with open(inFile) as input, open(outFile, "w") as output:
        for line in input:
            line = line.strip()
            line = re.sub("\t", "", line)
            if (lineCount % 2 == 1):
                header = line
            else:
                seq = re.sub(" ", "", line)
                output.write(revComp(header, seq))
                header = ""
                seq = ""
            lineCount +=1
    


def transcribe (header, sequence):
    lenSeq = getLength (sequence)
    sequence = getUpper(sequence)
    if (lineChecker(sequence)):
        sequence = re.sub("T", "U", sequence)
        return (header + "\n" + sequence + "\n")
    else:
        return (header + "\n" + "ERROR\n")

def sInput(inFile, outFile):
    lineCount = 1
    with open(inFile) as input, open(outFile, "w") as output:
        for line in input:
            line = line.strip()
            line = re.sub("\t", "", line)
            if (lineCount % 2 == 1):
                header = line
            else:
                seq = re.sub(" ", "", line)
                output.write(transcribe(header, seq))
                header = ""
                seq = ""
            lineCount +=1



def translate(header, sequence, codons_dict):
    sequence = getUpper(sequence)
    if (lineChecker(sequence)):
        protein = ""
        marker = 0
        rna = re.sub("T", "U", sequence)
        while marker < (getLength(rna) - 3):
            third = rna[marker:marker+3]
            seq = codons_dict[third]
            if (seq == "M"):
                protein +=seq
                marker +=3
            elif (seq!= "*" and protein[0] == "M"):
                protein += seq
                marker +=3
            elif (seq == "*"):
                marker = getLength(rna)
            else:
                marker +=1
        return(header + "\n" + protein + "\n")
    else:
        return (header + "\n" + "ERROR\n")

def lInput(inFile, outFile, codFile):
    lineCount = 1
    codons = {}
    with open(codFile) as input:
        for line in input:
            line = line.strip()
            x = line.split()
            codons[x[0]] = x[1]
    with open(inFile) as input, open(outFile, "w") as output:
        for line in input:
            line = line.strip()
            line = re.sub("\t", "", line)
            if (lineCount % 2 == 1):
                header = line
            else:
                seq = re.sub(" ", "", line)
                output.write(translate(header, seq, codons))
                header = ""
                seq = ""
            lineCount +=1



def countNucs (header, seq):
    lenSeq = getLength (seq)
    sequence = getUpper(seq)
    if (lineChecker(sequence)):
        aCount = str(getCount(sequence, "A"))
        tCount = str(getCount(sequence, "T"))
        cCount = str(getCount(sequence, "C"))
        gCount = str(getCount(sequence, "G")) 
        aPerc = aCount + "("+ str((float(aCount)/lenSeq) * 100) +"%)" + "\t"
        tPerc = tCount + "("+str((float(tCount)/lenSeq) * 100) +"%)" 
        cPerc = cCount + "(" + str((float(cCount)/lenSeq) * 100) +"%)" + "\t"
        gPerc = gCount + "("+str((float(gCount)/lenSeq) * 100) +"%)" + "\t" 
        return (header + "\t"+str(lenSeq)+ "\t"+ aPerc + cPerc + gPerc + tPerc + "\n")
    else:
        return (header + "\t" + "ERROR\n")    

def nInput(inFile, outFile):
    lineCount = 1 
    with open(inFile) as input, open(outFile, "w") as output:
       output.write("ID\tLength\tA(%A)\tC(%C)\tG(%G)\tT(%T)\n") 
       for line in input:
           line = line.strip()
           line = re.sub("\t", "", line)
           if (lineCount % 2 == 1):
               header = line[1:]
           else:
                seq = re.sub(" ", "", line)
                output.write(countNucs(header, seq))
                header = ""
                seq = ""
           lineCount +=1    



#start of main function line 215
 
if (getLength(sys.argv) <5):
    print("Error: USAGE: python myprogram.py: option –i input.fasta –o output.txt\n\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement\n\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides\n")
    sys.exit()

option = sys.argv[1]
inputLetter = sys.argv[2]
outputLetter = sys.argv[4]

if (option== "-g" or option == "-r" or option == "-s" or option == "-n"):
    if (getLength(sys.argv) !=6 or inputLetter != "-i" or outputLetter != "-o"):
        print("Error: USAGE: python myprogram.py: option –i input.fasta –o output.txt\n\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement\n\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides\n")
        sys.exit()
    else:
        inputFile = sys.argv[3]
        outputFile = sys.argv[5]
        if (option == "-g"):
            gInput(inputFile, outputFile)
        elif (option =="-r"):
            rInput(inputFile, outputFile)
        elif (option == "-s"):
            sInput(inputFile, outputFile)
        elif (option == "-n"):
            nInput(inputFile, outputFile)

elif (option =="-l"):
    if (getLength(sys.argv) !=8 or inputLetter != "-i" or outputLetter != "-o"):
        print("Error: USAGE: python myprogram.py: option –i input.fasta –o output.txt\n\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement\n\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides\n")
        sys.exit()
    else:
        inputFile = sys.argv[3]
        outputFile = sys.argv[5]
        codonFile = sys.argv[7]
        lInput(inputFile, outputFile, codonFile)
else: 
     print("Error: USAGE: python myprogram.py: option –i input.fasta –o output.txt\n\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement\n\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides")


