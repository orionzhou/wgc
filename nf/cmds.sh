s=maize
s=rice261
s=rice12
s=wheatD
s=common
s=Laburnicola

#excel.py tsv comparisons.xlsx --sheet $s 01.$s.tsv

for s in common maize rice261 rice12 wheatD Laburnicola
do
    excel.py tsv comparisons.xlsx --sheet $s 01.$s.tsv
done
