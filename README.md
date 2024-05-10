# Simplicial_vector_autoregressive_model
In this study, we address the scalability challenge of the vector autoregressive (VAR) model in modeling dynamic processes with numerous time series by leveraging data relations as inductive priors. We propose two simplicial VAR models: one for modeling time series on a single simplicial level (e.g., edge flows) and another for jointly modeling multiple time series across different levels, capturing their spatio-temporal interdependencies. These models employ simplicial convolutional filters to effectively capture structure-aware dependencies in a multiresolution manner.  Through experiments on synthetic and real-world networks, including water distribution networks and collaborative agent scenarios, we demonstrate the effectiveness of our approach, achieving competitive signal modeling accuracy with significantly fewer parameters compared to state-of-the-art alternatives.

# Experiment C. Collaborative agents
After downloading the repository, run the file 'Experiment_C_collaborative_agents.m'. This will enable you to compare the normalized mean squared error of SC-VAR with state-of-the-art competitors.

The experiment utilizes real dataset from https://datasets.simula.no/alfheim/, containing time-varying positional data of football players. This data is utilized to construct a simplicial complex (SC), where the average position of players is calculated, and the SC is constructed based on the distances between these average positions.In this experiment, the vertex signal represents speed, the edge signal corresponds to distance, and the triangle signal signifies area.

<img width="725" alt="Screenshot 2024-05-03 at 00 12 11" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/92a4bd2d-53fe-4da8-96f5-0a44c1e1d562">

# Experiment B. Water distribution network
After downloading the repository, run the file 'Experiment_B_Cherryhill.m'. This will enable you to compare the normalized mean squared error of SC-VAR with state-of-the-art competitors.

The experiment utilizes data generated through EPANET software for a real water distribution network. The network's structure is given, with vertex signal representing pressure and edge signal indicating flow rate.

<img width="886" alt="Screenshot 2024-05-08 at 22 52 19" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/1c8c7c6c-3f76-4411-86d7-1ea2088c8e94">

# Experiment A. Synthatic data
After downloading the repository:
(1) Run the file 'Experiment_A_Synthatic_1.m': This will enable you to compare the normalized mean squared error of SC-VAR with state-of-the-art competitors for time-varying model. 
(2) Run the file 'Experiment_A_Synthatic_2.m': This will enable you to compare the normalized mean squared error of SC-VAR with S-VAR for stationary model. 


<img width="644" alt="Screenshot 2024-05-09 at 01 52 22" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/210422fd-74c1-4853-a7ec-38a10ba0225f">

# Paramter comparison 
The SC-VAR and S-VAR models exhibit several orders of magnitude fewer parameters compared to TIRSO and RFNL-TIRSO, as illustrated in the graph (Note that y-axis is in log scale).

<img width="986" alt="Screenshot 2024-05-09 at 02 03 15" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/c83127ec-8047-4459-99fe-d0dcf712b839">




The MATLAB 2019b is employed to produce results for the experiments.
