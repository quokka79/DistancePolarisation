# DistancePolarisation
Calculates a 'DP' value for coordinate data from two channels indicating juxtaposition from +1 (identical) to 0 (randomly arranged) to -1 (perfectly dispersed).

Uses only coordinate (x,y) data.

Distance Polatisation Histogram. Here it was applied to STED data, from which clusters were identified. The centroid coordinates for each cluster was used to determine if the clusters in each channel were associated with each other. The association can work in both directions, e.g. a red cluster might always be found near a green cluster but perhaps green clusters have no strong interest in the red clusters, hence there is a '1-vs-2' and '2-vs-1' measurement for each image.

![Distance Polarisation Histogram](Histogram%20of%20Distance%20Polarisation.png?raw=true "Distance Polarisation Histogram")

Test Data ('sample' in the histogram) ... are red & green associated with each other?

![Test data](Sample%20-%20Composite.png?raw=true "Sample Data - Overlay")

Test Data - cluster centroids:

![Test data - centroids](sample_data.png?raw=true "Sample Data - Cluster centroids")

Positive Control Data ('control' in the histogram) ... the same protein stained with two different reporters, i.e. a well juxtaposed channels

![Highly juxtaposed data](Positive%20Control%20-%20Composite.png?raw=true "Positive Control (dual labelled) - Overlay")

Positive Control Data - cluster centroids:

![Positive Control data - centroids](ctrl_data.png?raw=true "Positive Control Data - Cluster centroids")
