#!/bin/bash

TIMES="500000 1000000" #"10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000"
SEEDS="1 2 3 4 5 6 7 8 9 10"


for t in $TIMES; do
    for s in $SEEDS; do
        echo "Working on experiment "$t", "$s
        if [ $t -gt 5000 ]
        then
            indelibleOptions=""
        else
            indelibleOptions="--indelibleSim --indelibleSim2"
        fi 
        echo $indelibleOptions
        python compareSimulators.py --nLeaves $t --createFasta $indelibleOptions --phastSim --useIndels --seed $s | tail -1 >> indel_data.txt
    done
done

# this line below will reformat the lines created above into something that can they be copied and pasted back into compareSimulators.py (as that script also does plotting)
awk -F "," '{print $1 $2 $3 $7 $8}' indel_data.txt | sed 's/[][]//g' | paste -s -d '         \n' - | awk -F ' ' '{print "[["$1 ", "$6 ", "$11 ", "$16 ", "$21 ", "$26 ", "$31 ", "$36 ", "$41 ", "$46 "], ["$2 ", "$7 ", "$12 ", "$17 ", "$22 ", "$27 ", "$32 ", "$37 ", "$42 ", "$47 "], ["$3 ", "$8 ", "$13 ", "$18 ", "$23 ", "$28 ", "$33 ", "$38 ", "$43 ", "$48 "], ["$4 ", "$9 ", "$14 ", "$19 ", "$24 ", "$29 ", "$34 ", "$39 ", "$44 ", "$49 "], ["$5 ", "$10 ", "$15 ", "$20 ", "$25 ", "$30 ", "$35 ", "$40 ", "$45 ", "$50" ]]"}' > indel_data_2.txt