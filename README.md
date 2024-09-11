# Dataset for the ICPE 2023 Data Challenge track

Information about the track:
[https://icpe2023.spec.org/tracks-and-submissions/data-challenge-track/](https://icpe2023.spec.org/tracks-and-submissions/data-challenge-track/)

The dataset contains performance measurements of JMH microbenchmarks from 30 Java open source projects. The list of projects, along with the revision at which the microbenchmarks were executed, can be found in [benchmarks_revision.csv](benchmarks_revision.csv).

The measurements are organized in time series available in the [timeseries](timeseries) folder. Morevover, the raw samples (JMH output) in JSON format can be found on [https://zenodo.org/record/5961018](https://zenodo.org/record/5961018) (~65GB when unpacked).

Questions about the dataset can be asked by opening issues on this repository, or by sendind an e-mail to icpe2023-data-challenge@easychair.org.

## Usage example

In the Python script [example_viz.py](example_viz.py) you can find an example of how to read the data and generate a simple plot for a random benchmark in the dataset. The script requires `pandas` and `matplotlib`:
```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
## For analysis data:

In /icpe/config.py

set the configs as your target running:
```Python
config = {
    "using_last_stable_segmet": True, # using when apply paper strategy, otherwise using the longest segment for classifying steady fork.
    "output_dir" : "./csv/output.csv",
    "figs" : "figs",
    "visualize": True
}
```
If the element "visualize" is set to True, then the figure fig.1 of the original paper. Otherwise, it will be ignore.

```bash
python run.py -i data --rlink "your revision/or result repo" --dlink "your data link" --action "your action" --range "range of data (file)"
```
For example:
```bash
python run.py -i data --rlink "./icpe/benchmarks_revision.csv" --dlink ".icpe/timeseries" --action "run" --range 0 1
```
If you action equal to debug, then it will run one by one. Someone can also config visual studio code debug configs as:
```sh
{
    "name": "Python: Debug1 run.py",
    "type": "python",
    "request": "launch",
    "program": "${workspaceFolder}/run.py", 
    "args": [
        "-i", "data",
        "--action", "debug",
    ],
    "console": "integratedTerminal"
    }
```

## For result analysis:
After getting the data analysis result, you can run the below code to visualize these result.
```
python run.py -i result -rlink "your revision/or result repo" dlink "your data link" --action "your action"
```
---
Example
1. Plot figures 6(a), 6(b) as the paper.
```sh
python run.py -i result --action rq1 --rlink ./csv/output.csv
```
2. Plot Figure 4 in the report (extra part).
```sh
python run.py -i result --action rq6 --rlink ./csv/output.csv ./csv/output_non_last_steady.csv
```
Find the results at /figs/ and /csv/.

The dataset was originally created for the paper:

Luca Traini, Vittorio Cortellessa, Daniele Di Pompeo, Michele Tucci  
**Towards effective assessment of steady state performance in Java software: Are we there yet?**  
Empirical Software Engineering (EMSE) - 28, 13 (2023)  
[https://doi.org/10.1007%2Fs10664-022-10247-x](https://doi.org/10.1007%2Fs10664-022-10247-x)  
[https://github.com/SEALABQualityGroup/steady-state](https://github.com/SEALABQualityGroup/steady-state)


