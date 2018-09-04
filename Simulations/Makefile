ART=/opt/art_bin_MountRainier/art_illumina
ARTOPTS = --seqSys HS25 --len 75 --noALN --samout
REF=~/Develop/Phix_mutagenesis/ExternalData/Illumina_WholeGenomeFasta/genome.fa

# simulate w/no errrors
COV=200000
OUTNAME =  Illumina_$(COV)_HS25_75

sim: gen mp ce
gen:
	/usr/bin/time $(ART) --fcov $(COV) --in $(REF) $(ARTOPTS) --out $(OUTNAME)
	/usr/bin/time samtools sort --reference $(REF) -O CRAM -m 5G  $(OUTNAME).sam  -o $(OUTNAME).cram
	samtools index  $(OUTNAME).cram

mp:
	samtools mpileup  $(OUTNAME).cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 35 -a | gzip --fast - > $(OUTNAME)_20_35.mpileup.gz &
	samtools mpileup  $(OUTNAME).cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 36 -a | gzip --fast - > $(OUTNAME)_20_36.mpileup.gz &
	samtools mpileup  $(OUTNAME).cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 37 -a | gzip --fast - > $(OUTNAME)_20_37.mpileup.gz &
	samtools mpileup  $(OUTNAME).cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 38 -a | gzip --fast - > $(OUTNAME)_20_38.mpileup.gz
	samtools mpileup  $(OUTNAME).cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 39 -a | gzip --fast - > $(OUTNAME)_20_39.mpileup.gz

ce:
	gunzip -c $(OUTNAME)_20_35.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $(OUTNAME)_20_35.tab &
	gunzip -c $(OUTNAME)_20_36.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $(OUTNAME)_20_36.tab &
	gunzip -c $(OUTNAME)_20_37.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $(OUTNAME)_20_37.tab &
	gunzip -c $(OUTNAME)_20_38.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $(OUTNAME)_20_38.tab
	gunzip -c $(OUTNAME)_20_39.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $(OUTNAME)_20_39.tab

syn:
	make _syn S=1pct COV=2000
	make _syn S=0.1pct COV=200
_syn:
	/usr/bin/time $(ART) --fcov $${COV} --in $${S}.fasta $(ARTOPTS) --out $${S}
	bowtie2 --threads 4 --very-fast --time --rg-id sim  -x $(REF) -U $${S}.fq \
		| samtools sort --reference $(REF) -O CRAM --threads 4 -m 5G - -o $${S}.cram
	rm -f $${S}.sam
	samtools index $${S}.cram
	samtools mpileup  $${S}.cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 35 -a | gzip --fast - > $${S}_20_35.mpileup.gz &
	samtools mpileup  $${S}.cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 36 -a | gzip --fast - > $${S}_20_36.mpileup.gz &
	samtools mpileup  $${S}.cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 37 -a | gzip --fast - > $${S}_20_37.mpileup.gz &
	samtools mpileup  $${S}.cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 38 -a | gzip --fast - > $${S}_20_38.mpileup.gz
	samtools mpileup  $${S}.cram -C 50 -d 1000000 -f $(REF) --min-MQ 20  --min-BQ 39 -a | gzip --fast - > $${S}_20_39.mpileup.gz
	gunzip -c $${S}_20_35.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $${S}_20_35.tab &
	gunzip -c $${S}_20_36.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $${S}_20_36.tab &
	gunzip -c $${S}_20_37.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $${S}_20_37.tab &
	gunzip -c $${S}_20_38.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $${S}_20_38.tab
	gunzip -c $${S}_20_39.mpileup.gz | python2.7 ~/Develop/python/Sequence/pileupCountAllAlleles.py | cut -f 1-8 > $${S}_20_39.tab







#download:
#	wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
#	tar zxf artbinmountrainier2016.06.05linux64.tgz
#	rm -f artbinmountrainier2016.06.05linux64.tgz
#	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz
#	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/NCBI/1993-04-28/PhiX_NCBI_1993-04-28.tar.gz
#	tar zxf PhiX_Illumina_RTA.tar.gz
#	tar zxf PhiX_NCBI_1993-04-28.tar.gz
#	rm -f PhiX_Illumina_RTA.tar.gz PhiX_NCBI_1993-04-28.tar.gz
#	mv PhiX/NCBI/1993-04-28/Sequence/WholeGenomeFasta NCBI_WholeGenomeFasta
#	mv PhiX/Illumina/RTA/Sequence/WholeGenomeFasta Illumina_WholeGenomeFasta
#	rm -fr PhiX


