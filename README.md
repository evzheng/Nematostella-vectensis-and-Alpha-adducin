# Nematostella vectensis and its gene family, Alpha-adducin

This repository is the methods of my Final Project of *Nematostella vectensis* and its gene family, Alpha-adducin. Feel free to reference the Final Paper pdf attached to the github repository prior to viewing the contents of the information below. Many contents were taken from Professor Joshua Rest.

**Contents**

1. Homolog alignment
2. Make phylogenetic tree for alpha-adducin
3. Reconcile species and gene tree
4. Graph the predicted Pfam domains onto the alpha-adducin tree phylogeny

## I. Homolog alignment

     ncbi-acc-download -F fasta -m protein XP_001634797.2

This downloads the file of the alpha-adducin protein.

The alpha-adducin sequence was obtained from NCBI by using blastp on the protein sequence.

     blastp -db allprotein.fas -query XP_001634797.2.fa -outfmt 0 -max_hsps 1 -out alpha_adducin.blastp.typical.out

To look at the output, type in:

     less alpha_adducin.blastp.typical.out
     
     blastp -db allprotein.fas -query XP_001634797.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out alpha_adducin.blastp.detail.out
     
This blast reformats our protein sequence.

To look at the output file, type in:

     less alpha_adducin.blastp.detail.out
     
To look at how many total human hits in the file, type in:

     grep -c Hsapiens alpha_adducin.blastp.detail.out
     
We get 3.

Then, awk was used to filter the initial blastp results with all of the homologs, even those with e values above 1e-14. By setting the e-value to 1e-14, the accurate hits of the protein homologs can be obtained.

     awk '{if ($6<0.00000000000001)print $1 }' alpha_adducin.blastp.detail.out > alpha_adducin.blastp.detail.filtered.out
     
     wc -l alpha_adducin.blastp.detail.filtered.out == 16
     
The above determines how many proteins are left in the final alignment after filtering for a reasonable e-value, with an e-value of 1e-14.

Word count line command was done for the original file (.detail.out file) and got 25 hits. Then, the word count line command was ran for the new filtered out file (.detail.filtered.out) and got 16 hits. The difference of the two ( 9) would be the total hits removed by this filter. This counts how many were filtered out.
     
Seqkit extracts the alpha-adducin proteins from the filtered out homologs.

     seqkit grep --pattern-file alpha_adducin.blastp.detail.filtered.out allprotein.fas > alpha_adducin.blastp.detail.filtered.fas

This command takes out alpha_adducin from the allprotein.fas file and then puts it into a new file called alpha_adducin.blastp.detail.filtered.fas.

Muscle was able to perform a global sequence alignment.

     muscle -in alpha_adducin.blastp.detail.filtered.fas -out alpha_adducin.blastp.detail.filtered.aligned.fas
     
This command does sequence alignment for everything within the file. It adds the gaps and outputs the new fasta file.

     t_coffee -other_pg seq_reformat -in alpha_adducin.blastp.detail.filtered.aligned.fas -output sim
     
This performs statistics on the new fasta file.

     alv -k alpha_adducin.blastp.detail.filtered.aligned.fas | less -r
     
     alv -kli --majority alpha_adducin.blastp.detail.filtered.aligned.fas | less -RS
     
Both these files view the alignments in alv.

T-Coffee performs statistics on the alpha-adducin alignment, removing positions in the alignment that were more than 50% gapped.

     t_coffee -other_pg seq_reformat -in alpha_adducin.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out alpha_adducin.blastp.detail.filtered.aligned.r50.fas

     alv -kli --majority alpha_adducin.blastp.detail.filtered.aligned.r50.fas | less -RS
     
This remove highly gapped positions in t_coffee.

This value and then compared to the previous most-recent alv one (alv -kli --majority alpha_adducin.blastp.detail.filtered.aligned.fas | less -RS) will give how many columns have been removed from the alignment.

The majority sequence is (926-592) 334. Looking at the alignment width can determine how many gaps are removed.

## II. Make phylogenetic tree for alpha-adducin

     cp ../lab5-L05-Evelyn/alpha_adducin.blastp.detail.filtered.aligned.fas . sed "s/ //g" alpha_adducin.blastp.detail.filtered.aligned.fas > alpha_adducin.blastp.detail.filtered.aligned.fas

After the gene family homologs for alpha-adducin were aligned, the phylogenetic tree was made. Iqtree obtained the unrooted alpha-adducin tree. 

     iqtree -s alpha_adducin.blastp.detail.filtered.aligned_.fas -nt 2
     
By using the unrooted tree, the midpoint rooted tree of alpha-adducin was generated using gotree.

     gotree reroot midpoint -i alpha_adducin.blastp.detail.filtered.aligned_.fas.treefile -o alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile
     
     nw_display -s alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile -w 1000 -b 'opacity:0' > alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile.svg
     
Lastly, nw_display was used to display the tree. 

## III. Reconcile species and gene tree

     cp ~/labs/lab6-L05-Evelyn/alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile
     
This copies over the treefile from gene's midpoint root.

First, notung reconciles the alpha-adducin tree with the species tree that includes Bilaterian, Cnidaria, Chordata, Echinodermata, and Mollusca.

     java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile --reconcile --speciestag prefix --savepng --events
     
This reconciles the alpha-adducin gene tree to species tree (speciesTreeBilateriaCnidaria.tre).

     python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile.reconciled --include.species
     
This creates a gene-species tree XML from alpha-adducin reconciled tree.

Notung root reroots the alpha-adducin gene tree to minimize the amount of duplications and deletions

     thirdkind -f alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile.reconciled.xml -o alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile.reconciled.svg
     
This converts XML to SVG to be pushed to github for viewing.

     java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile --root --speciestag prefix --savepng --events
     
This minimizes duplications and deletions.

Iqtree was then used to calculate the number of bootstrap support and tree search for the alpha-adducin gene. To evaluate the node (branch) support for the alpha-adducin gene tree using the bootstrap, the following commands were used.

     iqtree -s alpha_adducin.blastp.detail.filtered.aligned_.fas -bb 1000 -nt 2 -m VT+F+R4 -t alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile -pre alpha_adducin.genes.ufboot
     
Lastly, gotree reroots the bootstrap support tree into the midpoint root tree. This calculates the bootstrap support and performs tree search of alpha-adducin using fasta and gene tree files.

     gotree reroot midpoint -i alpha_adducin.genes.ufboot.treefile -o alpha_adducin.genes.ufboot.treefile.midpoint.ufboot
     
This reroots the tree to midpoint.

## IV. Graph the predicted Pfam domains onto the alpha-adducin tree phylogeny

After the alpha-adducin tree has been reconciled and rerooted to the midpoint root tree, the predicted Pfam domains can be graphed back to the alpha-adducin tree phylogeny. Iprscan5 analyzes the alpha-adducin gene sequences to the protein signature databases. 

     cd ~/labs/lab8-L05-Evelyn/myfamily

The sequence file ``alpha_adducin.blastp.detail.filtered.fas`` should also be there. 

     iprscan5   --email evelyn.zheng@stonybrook.edu  --multifasta --useSeqId --sequence   alpha_adducin.blastp.detail.filtered.fas

This command runs Iprscan5. The genes were then sorted according to the Pfam database.

     cat ~/labs/lab8-L05-Evelyn/myfamily/*.tsv.tsv > ~/labs/lab8-L05-Evelyn/alpha_adducin.domains.all.tsv

     grep Pfam ~/labs/lab8-L05-Evelyn/alpha_adducin.domains.all.tsv >  ~/labs/lab8-L05-Evelyn/alpha_adducin.domains.pfam.tsv

The following commands concatenate the files into a single file, then filters it by focusing only on the domains defined by the Pfam database.

     grep Pfam ~/labs/lab8-L05-Evelyn/alpha_adducin.domains.all.tsv >  ~/labs/lab8-L05-Evelyn/alpha_adducin.domains.pfam.tsv

This fixes the output by re-arranging the interproscan output. To locate the newick filefor the phylogeny associated with these sequences, choose the alpha-adducin rooted IQTREE gene tree from lab6 or lab7.

     cd ../lab7-L05-Evelyn/ cp alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile ~/labs/lab8-L05-Evelyn/ cd ../lab8-L05-Evelyn/
     
The midpoint treefile was used from lab 7. The name of the file is ``alpha_adducin.blastp.detail.filtered.aligned_.fas.midpoint.treefile``.

Finally, Evolview (an online web-service) plotted the alpha-adducin gene tree with the Pfam protein domains. To do that, go to the website: https://www.evolgenius.info/evolview/

Upload the tree by clicking on the File folder in the upper left hand corner.


To load the annotation file into Evolview, first update the annotation file by clicking on the upload tree datasets icon. This icon is the "up arrow" on the left hand side of the screen. Then, download it as a png or pdf and upload it to the github repository.

Lastly, to upload everything from the terminal into the github repository, type in the following:
     
     git add .
     git commit -a -m "Adding all new data files I generated in AWS to the repository."
     git pull --no-edit

Then, type:

     git push

Thank you for reading through my methods! The rest of the paper could be found in my pdf file of my Final Paper.
