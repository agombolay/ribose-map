#Output directory
output=$directory/Ribose-Map/Results/$reference/$sample/Coordinates

sed 's/chr1/chrI/'      $output/temp4.txt > $output/temp5.txt
sed 's/chr2/chrII/'     $output/temp4.txt > $output/temp5.txt
sed 's/chr3/chrIII/'    $output/temp4.txt > $output/temp5.txt
sed 's/chr4/chrIV/'     $output/temp4.txt > $output/temp5.txt
sed 's/chr5/chrV/'      $output/temp4.txt > $output/temp5.txt
sed 's/chr6/chrVI/'     $output/temp4.txt > $output/temp5.txt
sed 's/chr7/chrVII/'    $output/temp4.txt > $output/temp5.txt
sed 's/chr8/chrVIII/'   $output/temp4.txt > $output/temp5.txt
sed 's/chr9/chrIX/'     $output/temp4.txt > $output/temp5.txt
sed 's/chr10/chrX/'     $output/temp4.txt > $output/temp5.txt
sed 's/chr11/chrXI/'    $output/temp4.txt > $output/temp5.txt
sed 's/chr12/chrXII/'   $output/temp4.txt > $output/temp5.txt
sed 's/chr13/chrXIII/'  $output/temp4.txt > $output/temp5.txt
sed 's/chr14/chrXIV/'   $output/temp4.txt > $output/temp5.txt
sed 's/chr15/chrXV/'    $output/temp4.txt > $output/temp5.txt
sed 's/chr16/chrXVI/'   $output/temp4.txt > $output/temp5.txt
sed 's/chr17/chrXVII/'  $output/temp4.txt > $output/temp5.txt
sed 's/chr18/chrXVIII/' $output/temp4.txt > $output/temp5.txt
sed 's/chr19/chrXIX/'   $output/temp4.txt > $output/temp5.txt
sed 's/chr20/chrXX/'    $output/temp4.txt > $output/temp5.txt
sed 's/chr21/chrXXI/'   $output/temp4.txt > $output/temp5.txt
sed 's/chr22/chrXXII/'  $output/temp4.txt > $output/temp5.txt
