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
                exonerate --query x$count.fasta --target /jurg/group/SAM_RNA-SEQ_PIPELINE/GENOMES/ALL_TRS_160511.fsa --model ungapped --showalignment FALSE --bestn 1 > x$count.ALLtrs &
#echo "$i"
PID[$i]=$!
done
wait $PID
wait

