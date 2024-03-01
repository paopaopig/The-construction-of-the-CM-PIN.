Description of the project file: This project was compiled using python3.9.

Instructions for use: First, the refinement network CM-PIN of the given protein-protein interaction network can be obtained through 'main.py',
and then, the number of essential proteins identified by node ranking methods can be verified using the py files of 9-10 as follows.

1. All datasets used in the paper are in the 'datasets' folder.

2.The 'main.py' is the main function for building CM-PIN, which is used to output the refined network CM-PIN corresponding to the input network.
Among them, the 'refined_by_max_con.py' is used to extract the maximum connected subgraph of the given network; 
the 'fast_unfolding.py' is used to divide the modules of the protein-protein interaction network;
the 'calculate_pcc.py' ,'calculate_nsl.py' and 'calculate_topo.py' are used to calculate the three scores of each modules.

3. The 'network_based_centralities.py' is used to implement network-based centrality methods.

4. The 'PEC_method.py' is used to implement the PEC method.

5. The 'WDC_method.py' is used to implement the WDC method.

6. The 'construction_of_dpin.py' is used to build the DPIN.

7. The 'construction_of_rdpin.py' is used to build the RDPIN.

8. The 'identification_ep_on_spin_dpin_rdpin.py'  is used to identify the number of essential proteins for node ranking methods on SPIN,DPIN, and RDPIN.

9. The 'identification_of_essential_proteins.py' is used to identify the number of essential proteins for the network-based centrality methods on CMPIN.

10. The 'identification_ep_pec_wdc_on_CMPIN.py' is used to identify the number of essential proteins for the PEC and the WDC methods on CMPIN.
