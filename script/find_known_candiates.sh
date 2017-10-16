 fout=PAH_CHD.CHD.known.vcf.gz
 pref="${fout%.vcf.gz}"
 tabix ../../vcf/PAH_07102017.rawvariants.vcf.gz -R  CHD_candidate_genes.bed -h|bgzip -c > $fout
 bash /home/local/ARCS/nz2274/Pipeline/NA_script//mappability.annotation.sh 152 /home/local/ARCS/nz2274/Resources/mappability/hg19_152bp_mappability.bed.gz $fout
 mv $pref.mappability.vcf.gz $fout
 bash ~/Pipeline/NA_script/gnomad_coverage.annotation.sh $fout
 mv $pref.GnomADcov.vcf.gz $fout
 perl ~/Pipeline/NA_script/vcf2bed_full.pl  $fout
 file="$fout.2.txt"
 perl $ANNOVAR/table_annovar.pl $file $ANNHDB --buildver hg19 --remove -protocol revel,genomicSuperDups -operation f,r -otherinfo  -nastring .
 k=$(head -n 2 $fout.2.txt.hg19_multianno.txt|cut -f 6-7|sed 's/\t/\\t/g')
 awk -v h=$k 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){next}if(NR==2){$7="";$6=h;}print}' $fout.2.txt.hg19_multianno.txt >$fout.temp
 mv $fout.temp $fout.2.txt.hg19_multianno.txt
 Rscript  ~/PAH/Result/script/known_filter.R  ~/PAH/PAH-CHD/Data/$fout.2.txt.hg19_multianno.txt  ~/PAH/PAH-CHD/Data/PAH-CHD-list.csv  ~/PAH/PAH-CHD/Data/CHD.253kngenes.list  ~/PAH/PAH-CHD/Data/PAH-CHD.253KNG.csv