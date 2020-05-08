# Binder for BVCN Comparative Genomics lesson

Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!

[![Binder](https://mybinder.org/badge_logo.svg)](https://gesis.mybinder.org/binder/v2/gh/Arkadiy-Garber/bvcn-binder-paml/master?urlpath=lab)

Part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)


## Walkthrough

Enter the MagicCave

    cd MagicCave/

print the MagicLamp help menu

    MagicLamp.py help

print WspGenie help menu

    MagicLamp.py WspGenie -h

run WspGenie on test dataset

    MagicLamp.py WspGenie -bin_dir test_dataset/ -bin_ext fna -out wspgenie_out


go into the wspgenie output directory and check out the output file

    cd wspgenie_out/
    less -S wspgenie-summary.csv

check out the gene predictions

    cd ORF_calls/
    cd ../../

mv ORF calls to the main directory

    mv wspgenie_out/ORF_calls/ ./

print LithoGenie help menu

    MagicLamp.py LithoGenie -h

run LithoGenie on ORF calls

    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out

check out the output

    cd lithogenie_out/
    less -S lithogenie-summary.csv
    less lithogenie.ALL.heatmap.csv
    cd ../

re-run LithoGenie to create a .heatmap.csv for an element-of-interest

    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out --skip -cat sulfur
    # answer 'y' to the question
    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out --skip -cat iron

check out the new results

    cd lithogenie_out/
    less lithogenie.sulfur.heatmap.csv
    less lithogenie.iron.heatmap.csv

print the HmmGenie help menu

    MagicLamp.py HmmGenie -h

run HmmGenie with a set of HMMs for gas vesicle formation

    MagicLamp.py HmmGenie -hmm_dir MagicCave/hmms/gas/ -hmm_ext hmm -bin_dir test_dataset/ -bin_ext fna -out gas_out

check out the results and re-run HmmGenie with more stringent parameters

    MagicLamp.py HmmGenie -hmm_dir MagicCave/hmms/gas/ -hmm_ext hmm -bin_dir test_dataset/ -bin_ext fna -out gas_out -clu 5

check out the results

    cd gas_out/
    less -S genie-summary.csv
