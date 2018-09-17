#!/usr/bin/env bash

for i in /path/to/deseq/*.rnk; do

echo ${i::-4}

java -cp ~/gsea-3.0.jar -Xmx1024m xtools.gsea.GseaPreranked \
	-rnk ${i} \
	-gmx /path/to/h.all.v6.2.symbols.gmt \
	-rpt_label ${i::-4} \
	-zip_report false \
	-out /path/to/gsea_out -gui false

done

echo Done all
