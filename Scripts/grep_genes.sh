#!/bin/bash

PROJECT="/home/daniele/Desktop/IV_course/II_semester/TF_analysis";
HOMER_RESULTS="${PROJECT}/homerResults";
MOTIFS=`ls ${HOMER_RESULTS} | grep -E '^motif[0-9]+.motif$'`;

> "$PROJECT/genes.txt"

for	motif in $MOTIFS;
    do
        motif="$HOMER_RESULTS/$motif";
        grepped_motif=`head -n 1 $motif | awk '{print $2}' | sed 's/:/ /g' |
                       awk '{print $2}' | sed 's/\// /g' | awk '{print $1}'`;

        if [[ $grepped_motif =~ [()] ]]; then
            fixed_motifs=`echo $grepped_motif | sed 's/(/ /g' |
                          awk '{print $1}'`;
            echo $fixed_motifs >> "$PROJECT/genes.txt";
        else
            echo $grepped_motif >> "$PROJECT/genes.txt"
        fi
    done