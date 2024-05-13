This repository contains the code for the experiments outlined in the paper  " Simplicial vector autoregressive models" .
https://www.techrxiv.org/users/779840/articles/921287-simplicial-vector-autoregressive-models

# Simplicial_vector_autoregressive_model
The VAR model is widely used for dynamic processes, but handling numerous time series becomes challenging due to the explosion of parameters. To tackle this, we suggest using data relations as priors to cope with dimensionality while still effectively modeling time series. We explore using simplicial complexes as priors for time series on higher-level networks like edges and triangles. We introduce two types of simplicial VAR models: one for time series on a single simplicial level (e.g., edge flows), and another for modeling multiple time series across different simplicial levels, capturing their spatio-temporal interdependencies. These models employ simplicial convolutional filters to share parameters and capture structure-aware dependencies at various resolutions. Additionally, we develop a joint simplicial-temporal Fourier transform to analyze their spectral characteristics. For streaming signals, we design an online learning algorithm for simplicial VAR models, proving its convergence under reasonable assumptions. Finally, experiments on synthetic networks, water distribution networks, and collaborative agents confirm that our models achieve competitive accuracy with far fewer parameters than current methods.

# Experiment C. Collaborative agents
After downloading the repository, run the file 'Experiment_C_collaborative_agents.m'. This will enable you to compare the normalized mean squared error of SC-VAR with state-of-the-art competitors.

The experiment utilizes real dataset from https://datasets.simula.no/alfheim/, containing time-varying positional data of football players. This data is utilized to construct a simplicial complex (SC), where the average position of players is calculated, and the SC is constructed based on the distances between these average positions.In this experiment, the vertex signal represents speed, the edge signal corresponds to distance, and the triangle signal signifies area.


<img width="952" alt="Screenshot 2024-05-11 at 22 53 13" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/6b61c52e-5a09-4673-be97-5878de7a02eb">


# Experiment B. Water distribution network
After downloading the repository, run the file 'Experiment_B_Cherryhill.m'. This will enable you to compare the normalized mean squared error of SC-VAR with state-of-the-art competitors.

The experiment utilizes data generated through EPANET software for a real water distribution network. The network's structure is given, with vertex signal representing pressure and edge signal indicating flow rate.

<img width="841" alt="Screenshot 2024-05-08 at 22 52 19" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/58653726-dd2a-4913-aae6-358fe0672f00">


# Experiment A. Synthetic data
After downloading the repository:
(1) Run the file 'Experiment_A_Synthetic_1.m': This will enable you to compare the normalized mean squared error of SC-VAR with state-of-the-art competitors for time-varying data generation model. 
(2) Run the file 'Experiment_A_Synthetic_2.m': This will enable you to compare the normalized mean squared error of SC-VAR with S-VAR for stationary data generation model. 


<img width="644" alt="Screenshot 2024-05-09 at 01 52 22" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/210422fd-74c1-4853-a7ec-38a10ba0225f">

# Parameter comparison 
The SC-VAR and S-VAR models exhibit several orders of magnitude fewer parameters compared to TIRSO and RFNL-TIRSO, as illustrated in the graph (Note that y-axis is in log scale).

<img width="986" alt="Screenshot 2024-05-09 at 02 03 15" src="https://github.com/rohantmoney/Simplicial-vector-autoregressive-model/assets/61416415/c83127ec-8047-4459-99fe-d0dcf712b839">




    ****************** The aforementioned experiments were carried out using MATLAB 2019b ********************
