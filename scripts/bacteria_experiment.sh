#!/bin/bash
#shopt -s extglob


TIMES="100 1000 10000 100000 1000000"
SEEDS="1 2 3 4 5 6 7 8 9 10"

for t in $TIMES; do
    for s in $SEEDS; do
        if [ $t -gt 1001 ]
        then
            seqgenOptions=""
        else
            seqgenOptions="--seqgenSim"
        fi         
        echo experiment $t $s $seqgenOptions
        python ../scripts/compareSimulators.py --nLeaves $t $seqgenOptions --phastSim --noHierarchy --seed $s --reference ./sequence.fasta | tail -1 >> indel_data.txt
        #rm -v !("sequence.fasta"|"indel_data.txt")
    done
done

# this line below will reformat the lines created above into something that can they be copied and pasted back into compareSimulators.py (as that script also does plotting)
awk -F "," '{print $1 $2 $4 $6}' indel_data.txt | sed 's/[][]//g' | paste -s -d '         \n' - | awk -F ' ' '{print "[["$1 ", "$5 ", "$9 ", "$13 ", "$17 ", "$21 ", "$25 ", "$29 ", "$33 ", "$37 "], ["$2 ", "$6 ", "$10 ", "$14 ", "$18 ", "$22 ", "$26 ", "$30 ", "$34 ", "$38 "], ["$3 ", "$7 ", "$11 ", "$15 ", "$19 ", "$23 ", "$27 ", "$31 ", "$35 ", "$39 "], ["$4 ", "$8 ", "$12 ", "$16 ", "$20 ", "$24 ", "$28 ", "$32 ", "$36 ", "$40 "]]"}' > indel_data_2.txt
