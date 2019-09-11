 cd /data2/shane/Documents/SLC_id

for i in /data2/shane/Documents/SLC_id/phylogeny/phylip/*.phy
do
  b=$(echo $(basename $i) | cut -d '_' -f 1,2) 
  raxfile=$i
  raxdir=/data2/shane/Documents/SLC_id/SLC_phylogeny/
  #rm ./SLC_phylogeny/RAxML*
  /data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T 12 -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
  #/data2/shane/Applications/standard-RAxML-master/raxmlHPC-AVX -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAAUTO -s $raxfile -n $i'.tre' -w $raxdir ## LOCAL
done