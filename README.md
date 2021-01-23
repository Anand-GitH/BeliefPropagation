# BeliefPropagation
Belief Propagation over a Bayesian Network

Belied Propagation: If there is a variable X in a cluster and it carries its own
belief of the variable X and even it applies to other clusters having X. There exists
a unique path between such clusters. If X is changed in one of the clusters we
can verify how it will affect in other clusters as there exists a path and each
cluster starts updating the new belief of X and act accordingly and we can verify
the new probability of the cluster and determines predictions.


Packages Used: gRain, igraph, ggm
Dataset: Coronary artery disease data
Methods: dSep, cptable, grain, compileCPT, extractCPT, setFinding, getFinding, querygrain,
simulate.grain

Goal: Construct the network, identify nodes that are d-separable, structure learning and parameter
learning and then belied propagation and then simulation of the model over a data to verify results and
to make predictions.

Newly created functions:
ffinddsepingraph- function which returns all possible d-separations in the given input graph

CAD dataset contains 14 columns but we will use just 6 columns to create a network and check belief
propagation


Outcome of 500 observations:
Probability of Heart Failure – 46%
Probability of CAD - 55%

Observation: Predicted data for 25 observations and 500 observations is nearly same at least for CAD.
With graphical modelling, structure and parameter learning and belief propagation we can predict the
results which is nearly to reality. Predicted results can be put into use to device action plan for public
safety measures or raising awareness of high cholesterol levels in females and national nutrition can
work out dietary plan.
