awk '{print $1,":",$4, $2}' dbsnp37.txt > dbsnp37_chpos.txt
sed -i -e 's/ : /:/g' dbsnp37_chpos.txt
sort -k2 dbsnp37_chpos.txt > dbsnp37_chrpos_sorted.txt

sort sorted_rsid_list.txt > rsid_list.txt

join -1 2 -2 1 dbsnp37_chpos.txt rsid_list.txt > rsid_37build.txt
join -1 2 -2 1 dbsnp37_chrpos_sorted.txt rsid_list.txt > rsid_37build.txt

