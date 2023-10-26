# edna-species-lists

This repository contains species lists for fish, mammals, and turtles for the marine World Heritage sites and feeds directly into the [eDNA Expeditions sample tracking website](https://samples.ednaexpeditions.org/). These lists are a combination of observations from OBIS as well as detections from eDNA sampling.

Species lists based on OBIS are generated in the [mwhs-obis-species](https://github.com/iobis/mwhs-obis-species) repository, eDNA results are uploaded to this repository (folder [pipeline_results](pipeline_results)), and the combined lists are generated in this repository using a GitHub action (folder [lists](lists)). The lists are available as CSV and JSON.