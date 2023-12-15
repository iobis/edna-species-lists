# edna-species-lists

This repository contains species lists for fish, mammals, and turtles for the marine World Heritage sites and feeds directly into the [eDNA Expeditions sample tracking website](https://samples.ednaexpeditions.org/). These lists are a combination of observations from OBIS as well as detections from eDNA sampling.

Species lists based on OBIS are generated in the [mwhs-obis-species](https://github.com/iobis/mwhs-obis-species) repository, eDNA results are generated in [edna-results](https://github.com/iobis/edna-results), and combined lists are generated using a Python script in this repository as well (folder [lists](lists)). The lists are available as CSV and JSON.


## Download data

```
rm -r edna-results
git clone --branch data --depth 1 git@github.com:iobis/edna-results.git
```
