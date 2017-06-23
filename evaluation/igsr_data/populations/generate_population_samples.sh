#1
#African-Caribbean
rm acb.txt
grep ACB *.tsv | awk '{ print $1 }' > acb.txt

#2
#African-American SW
rm asw.txt
grep ASW *.tsv | awk '{ print $1 }' > asw.txt

#3
#Bengali
rm beb.txt
grep BEB *.tsv | awk '{ print $1 }' > beb.txt

#4
#Dai Chinese
rm cdx.txt
grep CDX *.tsv | awk '{ print $1 }' > cdx.txt

#5
#CEPH
rm ceu.txt
grep CEU *.tsv | awk '{ print $1 }' > ceu.txt

#6
#Han Chinese
rm chb.txt
grep CHB *.tsv | awk '{ print $1 }' > chb.txt

#7
#Denver Chinese
rm chd.txt
grep CHD *.tsv | awk '{ print $1 }' > chd.txt

#8
#Southern Han Chinese
rm chs.txt
grep CHS *.tsv | awk '{ print $1 }' > chs.txt

#9
#Colombian
rm clm.txt
grep CLM *.tsv | awk '{ print $1 }' > clm.txt

#10
#Esan
rm esn.txt
grep ESN *.tsv | awk '{ print $1 }' > esn.txt

#11
#Finnish
rm fin.txt
grep FIN *.tsv | awk '{ print $1 }' > fin.txt

#12
#British
rm gbr.txt
grep GBR *.tsv | awk '{ print $1 }' > gbr.txt

#13
#Gujarati
rm gih.txt
grep GIH *.tsv | awk '{ print $1 }' > gih.txt

#14
#Gambian Mandinka
rm gwd.txt
grep GWD *.tsv | awk '{ print $1 }' > gwd.txt

#14
#Gambian Fula
rm gwf.txt
grep GWF *.tsv | awk '{ print $1 }' > gwf.txt

#15
#Gambian Jola
rm gwj.txt
grep GWJ *.tsv | awk '{ print $1 }' > gwj.txt

#15
#Gambian Wolof
rm gww.txt
grep GWW *.tsv | awk '{ print $1 }' > gww.txt

#16
#Spanish
rm ibs.txt
grep IBS *.tsv | awk '{ print $1 }' > ibs.txt

#17
#Indian
rm itu.txt
grep ITU *.tsv | awk '{ print $1 }' > itu.txt

#18
#Japanese
rm jpt.txt
grep JPT *.tsv | awk '{ print $1 }' > jpt.txt

#19
#Kinh Vietnamese
rm khv.txt
grep KHV *.tsv | awk '{ print $1 }' > khv.txt

#20
#Luhya
rm lwk.txt
grep LWK *.tsv | awk '{ print $1 }' > lwk.txt

#21
#Mende
rm msl.txt
grep MSL *.tsv | awk '{ print $1 }' > msl.txt

#22
#Mexican-American
rm mxl.txt
grep MXL *.tsv | awk '{ print $1 }' > mxl.txt

#23
#Peruvian
rm pel.txt
grep PEL *.tsv | awk '{ print $1 }' > pel.txt

#24
#Punjabi
rm pjl.txt
grep PJL *.tsv | awk '{ print $1 }' > pjl.txt

#25
#Puerto Rican
rm pur.txt
grep PUR *.tsv | awk '{ print $1 }' > pur.txt

#26
#Sri Lankan
rm stu.txt
grep STU *.tsv | awk '{ print $1 }' > stu.txt

#27
#Tuscan
rm tsi.txt
grep TSI *.tsv | awk '{ print $1 }' > tsi.txt

#28
#Yoruba
rm yri.txt
grep YRI *.tsv | awk '{ print $1 }' > yri.txt


wc -l *.txt
wc -l *.tsv

