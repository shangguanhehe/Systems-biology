from Bio import SeqIO
from collections import OrderedDict
records = [r for r in SeqIO.parse("F:/data/wheat.fa", "fasta")]
fre_file = open('wheat_codon.txt', 'w')
fre_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('name', "AAA", "TAA", "CAA", "GAA", "AAC", "TAC", "CAC", "GAC", "AAG", "TAG", "CAG", "GAG", "AAT", "TAT", "CAT", "GAT", "ACA", "TCA", "CCA", "GCA", "ACC", "TCC", "CCC", "GCC", "ACG", "TCG", "CCG", "GCG", "ACT", "TCT", "CCT", "GCT", "AGA", "TGA", "CGA", "GGA", "AGC", "TGC", "CGC", "GGC", "AGG", "TGG", "CGG", "GGG", "AGT", "TGT", "CGT", "GGT", "ATA", "TTA", "CTA", "GTA", "ATC", "TTC", "CTC", "GTC", "ATG", "TTG", "CTG", "GTG", "ATT", "TTT", "CTT", "GTT", "total", "AAAfre", "TAAfre", "CAAfre", "GAAfre", "AACfre", "TACfre", "CACfre", "GACfre", "AAGfre", "TAGfre", "CAGfre", "GAGfre", "AATfre", "TATfre", "CATfre", "GATfre", "ACAfre", "TCAfre", "CCAfre", "GCAfre", "ACCfre", "TCCfre", "CCCfre", "GCCfre", "ACGfre", "TCGfre", "CCGfre", "GCGfre", "ACTfre", "TCTfre", "CCTfre", "GCTfre", "AGAfre", "TGAfre", "CGAfre", "GGAfre", "AGCfre", "TGCfre", "CGCfre", "GGCfre", "AGGfre", "TGGfre", "CGGfre", "GGGfre", "AGTfre", "TGTfre", "CGTfre", "GGTfre", "ATAfre", "TTAfre", "CTAfre", "GTAfre", "ATCfre", "TTCfre", "CTCfre", "GTCfre", "ATGfre", "TTGfre", "CTGfre", "GTGfre", "ATTfre", "TTTfre", "CTTfre", "GTTfre"))
for i in records:
    CodonsDict = OrderedDict([("AAA", 0), ("TAA", 0), ("CAA", 0), ("GAA", 0), ("AAC", 0), ("TAC", 0), ("CAC", 0), ("GAC", 0), ("AAG", 0), ("TAG", 0), ("CAG", 0), ("GAG", 0), ("AAT", 0), ("TAT", 0), ("CAT", 0), ("GAT", 0), ("ACA", 0), ("TCA", 0), ("CCA", 0), ("GCA", 0), ("ACC", 0), ("TCC", 0), ("CCC", 0), ("GCC", 0), ("ACG", 0), ("TCG", 0), ("CCG", 0), ("GCG", 0), ("ACT", 0), ("TCT", 0), ("CCT", 0), ("GCT", 0), ("AGA", 0), ("TGA", 0), ("CGA", 0), ("GGA", 0), ("AGC", 0), ("TGC", 0), ("CGC", 0), ("GGC", 0), ("AGG", 0), ("TGG", 0), ("CGG", 0), ("GGG", 0), ("AGT", 0), ("TGT", 0), ("CGT", 0), ("GGT", 0), ("ATA", 0), ("TTA", 0), ("CTA", 0), ("GTA", 0), ("ATC", 0), ("TTC", 0), ("CTC", 0), ("GTC", 0), ("ATG", 0), ("TTG", 0), ("CTG", 0), ("GTG", 0), ("ATT", 0), ("TTT", 0), ("CTT", 0), ("GTT", 0)])
    if i.seq.startswith('ATG') or i.seq.startswith('ATG'):
        for j in range(0, len(str(i.seq)), 3):
            codon = str(i.seq)[j:j+3]
            if codon in CodonsDict.keys():
                CodonsDict[codon] += 1
        total = sum([CodonsDict[key] for key in CodonsDict.keys()])
        fre_file.writelines('%s\t' % i.id)
        for key in CodonsDict.keys():
            fre_file.writelines('%i\t' % (CodonsDict[key]))
        fre_file.writelines('%s\t' % total)
        for key in CodonsDict.keys():
            fre_file.writelines('%f\t' % (CodonsDict[key] / float(total)))
        fre_file.writelines('\n')
fre_file.close()
