# edna-species-lists

This repository contains species lists for the marine World Heritage sites and feeds directly into the [eDNA Expeditions sample tracking website](https://samples.ednaexpeditions.org/). These lists are a combination of observations from OBIS as well as detections from eDNA sampling.

Species lists based on OBIS are generated in the [mwhs-obis-species](https://github.com/iobis/mwhs-obis-species) repository, eDNA results are generated in [edna-results](https://github.com/iobis/edna-results), and combined lists are generated using a Python script in this repository as well (folders [lists](lists) and [lists_full](lists_full)). The lists are available as CSV and JSON.


## For maintainers

To download the eDNA dataset:

```
rm -r output.zip output
wget https://obis-edna-results.s3.amazonaws.com/output.zip
unzip output.zip
```
