
Example Analyses
=================

    This section includes multiple example stepwise analyses of data from a
    qCLASH experiment described in [Gay2018]_, with data acquired from
    the NCBI Gene Expression Omnibus (GEO) accession
    `GSE101978 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978>`_.

    Each analysis is implemented both using the Python3 API, and as a sequence of shell
    executable commands in a bash script. The Python API implementations are generally
    significantly more efficient as more steps can be performed on a single iteration
    over the input data.

    Each analysis performs quality control steps on the data by checking data integrity
    (:ref:`hyb_check`) and removing artifactual
    ribosomal- and mitochondrial-RNA hybrids (:ref:`hyb_filter`).
    Further filtration may be performed, and then each described analysis is carried out.

        ====================================== ===================================================
        Pipeline                               Description
        ====================================== ===================================================
        :ref:`Example Type-miRNA Analysis`     Quantify the sequence and miRNA types in a hyb file
        :ref:`Example Target Analysis`         Analyze targets of a set of miRNAs from a single
                                               experiment
        :ref:`Example Grouped Target Analysis` Analyze and plot targets of a set of miRNAs from
                                               pooled experimental replicates
        :ref:`Example Fold Analysis`           Analyze and plot predicted miRNA
                                               folding patterns in
                                               miRNA-containing hybrids
        ====================================== ===================================================


    Further details on each respective example analysis can be found in each section.

.. toctree::
   :maxdepth: 4

   /example_01_link
   /example_02_link
   /example_03_link
   /example_04_link

