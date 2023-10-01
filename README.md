# BIOFUTURE
BIOFUTURE is a decision-support framework for sustainable Waste-to-X biorefinery supply chains, proposed and developed as part of my Master Thesis in Industrial Engineering and Management at FEUP. 

This package is a work in progress. BIOFUTURE will be in continuous development and all contributions are welcomed.

# Framework

BIOFUTURE (Biorefinery supply-chaIn Optimization Framework for sUsTainable Upcycling of oRganic wastE)assists in the analysis of the optimal design and planning of biorefinery supplychains in three main interactive stages that feed information into each other: (i) Selection of Supply Chain Network, (ii) Supply, Demand, and Price Forecasting; and (iii) Supply-Chain MultiObjective Optimization. Each stage is further characterized by its own step-by-step methodology, which requires decision-maker interaction with the open-source tool.

![framework](/assets/BIOFUTURE_Framework.png)

## 1. Selection of Supply Chain Network

To reduce the complexity of supply chain design optimization, a method for pre-filtering waste feedstocks, suppliers, products, and markets was developed. As an open source tool for the efficient evaluation of biorefining supply chains, this framework focuses on creating a repeatable systematic approach to leverage publicly and readily available data. Therefore, BIOFUTURE is directly integrated with Eurostat as it was considered the most reliable, rich, and frequently updated open-source database, which includes all the necessary data for such a framework. This stage of the methodology is specific to the exploration of Eurostat’s database and culminates in the creation of high quality datasets with feedstock and product data.

![framework](/assets/BIOFUTURE_Stage1.png)

## 2. Supply, Demand, and Price Forecasting

By integrating various data analysis, preprocessing, and modeling techniques, robust product demands and prices’ estimates are used as input for a multi-period biorefinery supply chain optimization model. The nature of the process implemented ensures that the methodology is adaptable to different products and requirements, staying up-to-date with the most recent available data in Eurostat, making it a versatile tool for other biorefining or supply chain studies. Models implemented are ARIMA and LSTM.

![framework](/assets/BIOFUTURE_Stage2.png)

## 3. Supply-Chain MultiObjective Optimization

Superstructure optimization is a powerful but computationally demanding task that can be used to select the optimal flowsheet among many alternatives within a single optimization. It includes all the alternative processes in the design space, as well as the feedstocks and main outputs. All the links between these elements also have to be defined. BIOFUTURE includes a superstructure optimization module implemented as a mixed-integer linear programming model (MILP) via the Pyomo modeling language, which allows the use of multiple solvers.

![framework](/assets/BIOFUTURE_Stage3.png)

