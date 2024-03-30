#########################################################################
# Author: Xingcheng Lin
# Created Time: Sun Dec  5 11:18:36 2021
# File Name: cmd.copyprotein_list.sh
# Description: proteinList.txt includes all the strong binder labels and the template PDB
#########################################################################
#!/bin/bash

# We copy the updated phi1 parameters into the protocol
cp phi1_list.txt optimization/for_bindingE/2bu1/
cp phi1_list.txt optimization/for_bindingE/loocv/
