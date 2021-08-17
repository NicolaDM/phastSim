#!/bin/bash

TIMES="10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000"
SEEDS="1 2 3 4 5 6 7 8 9 10"


for t in $TIMES; do
    for s in $SEEDS; do
        echo "Working on experiment "$t", "$s
        if [ $t -gt 51 ]
        then
            pyvolveOptions=""
        else
            pyvolveOptions="--pyvolveSim"
        fi 
        echo $pyvolveOptions
        
        if [ $t -gt 10001 ]
        then
            indelibleOptions=""
        else
            indelibleOptions="--indelibleSim --indelibleSim2"
        fi 
        echo $indelibleOptions
        
        if [ $t -gt 20001 ]
        then
            seqgenOptions=""
        else
            seqgenOptions="--seqgenSim"
        fi         
        echo $seqgenOptions
        
        if [ $t -gt 50001 ]
        then
            phastSimOptions="--phastSim"
        else
            phastSimOptions="--phastSim --createFasta"
        fi 
        echo $phastSimOptions
        python compareSimulators.py --nLeaves $t $pyvolveOptions $seqgenOptions $phastSimOptions $indelibleOptions --seed $s | tail -1 >> indel_data.txt
    done
done

# this line below will reformat the lines created above into something that can they be copied and pasted back into compareSimulators.py (as that script also does plotting)
awk -F "," '{print $1 $2 $3 $7 $8}' indel_data.txt | sed 's/[][]//g' | paste -s -d '         \n' - | awk -F ' ' '{print "[["$1 ", "$9 ", "$17 ", "$25 ", "$33 ", "$41 ", "$49 ", "$57 ", "$65 ", "$73 "], ["$2 ", "$10 ", "$18 ", "$26 ", "$34 ", "$42 ", "$50 ", "$58 ", "$66 ", "$74 "], ["$3 ", "$11 ", "$19 ", "$27 ", "$35 ", "$43 ", "$51 ", "$59 ", "$67 ", "$75 "], ["$4 ", "$12 ", "$20 ", "$28 ", "$36 ", "$44 ", "$52 ", "$60 ", "$68 ", "$76 "], ["$5 ", "$13 ", "$21 ", "$29 ", "$37 ", "$45 ", "$53 ", "$61 ", "$69 ", "$77" ], ["$6 ", "$14 ", "$22 ", "$30 ", "$38 ", "$46 ", "$54 ", "$62 ", "$70 ", "$78 "], ["$7 ", "$15 ", "$23 ", "$31 ", "$39 ", "$47 ", "$55 ", "$63 ", "$71 ", "$79 "],["$8 ", "$16 ", "$24 ", "$32 ", "$40 ", "$48 ", "$56 ", "$64 ", "$72 ", "$80 "]]"}' > indel_data_2.txt