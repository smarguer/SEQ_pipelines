 bits=`ls x* | wc -l`
#echo $bits
        let "bits=$bits-1"
#echo $bits
        echo "mapping reads"

        for i in `seq 0 $bits`
        do
                if [ $i -lt 10 ]
                then
                        count="0$i"
                else
                        count="$i"
                fi
                mv x$count x$count.fasta
#oowriter x$count.fasta &
#oowriter &
<<<<<<< HEAD
                exonerate --query x$count.fasta --target /jurg/group/SAM_RNA-SEQ_PIPELINE/GENOMES/ALLchromosomes.090511.ERCC.fsa --model ungapped --showalignment FALSE --bestn 1 > x$count.ALLchr &
=======
                exonerate --query x$count.fasta --target ~/POMBE_SEQ/analysis/GENOMES/ALLchromosomes.090511.ERCC.fsa --model ungapped --showalignment FALSE --bestn 1 > x$count.ALLchr &
>>>>>>> 2b5f1ec794c8b158eadf5de448c55d5efaa080e1
#echo "$i"
PID[$i]=$!
done
wait $PID
wait

