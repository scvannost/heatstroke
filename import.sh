for i in {1..48}
do
if (($i < 10)); then
        scp 'scv@luria.mit.edu:~/data/mapped/0'$i'/0'$i'.isoforms.results' ~/heat_stroke/mapped/$i.isoforms.results
        scp 'scv@luria.mit.edu:~/data/mapped/0'$i'/0'$i'.genes.results' ~/heat_stroke/mapped/$i.genes.results
else
        scp 'scv@luria.mit.edu:~/data/mapped/'$i'/'$i'.isoforms.results' ~/heat_stroke/mapped/$i.isoforms.results
        scp 'scv@luria.mit.edu:~/data/mapped/'$i'/'$i'.genes.results' ~/heat_stroke/mapped/$i.genes.results
fi
done
