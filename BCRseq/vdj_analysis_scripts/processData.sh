# Kenneth B. Hoehn
# 2/1/2021
# Run BCR data through immcantation pipeline
# First, open immcantation Docker image to run data
# sudo docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.0.0 bash
# toggle which sets to run by commenting out dirs lines
 
dirs=("Pam1_CITE_multi5P06" "Pam1_CITE_multi5P07" "Pam1_CITE_multi5P08")

for wdir in "${dirs[@]}"
do
	fdir="../raw/$wdir/vdj_b"
	echo $wdir
	echo $fdir
	changeo-10x -s $fdir/filtered_contig.fasta -a $fdir/filtered_contig_annotations.csv -o . \
    	-g human -t ig -p 3 -o $wdir
done

