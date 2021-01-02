from Bio import SeqIO
in_file = open("F:/data/wheat.fa", "r")
out_file = open("wheat.txt", "w")
out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('name', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'total', 'Afre', 'Cfre', 'Dfre', 'Efre', 'Ffre', 'Gfre', 'Hfre', 'Ifre', 'Kfre', 'Lfre', 'Mfre', 'Nfre', 'Pfre', 'Qfre', 'Rfre', 'Sfre', 'Tfre', 'Vfre', 'Wfre', 'Yfre'))
for record in SeqIO.parse(in_file, "fasta"):
    if not record.seq.startswith('M'):
        continue
    total = float(len(record))
    amino = record.seq
    out_file.write('%s\t' % record.id)
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        n = amino.count(aa)
        out_file.write('%i\t' % n)
    out_file.write('%i\t' % total)
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        n = amino.count(aa)
        if total != 0.00:
            frq = n/total
            out_file.write('%f\t' % frq)
    out_file.writelines('\n')
in_file.close()
out_file.close()
