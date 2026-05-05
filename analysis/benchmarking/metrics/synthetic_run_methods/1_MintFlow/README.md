

To run MintFlow on synthetic data, please:
- Install MintFlow v0.3.0 (or the latest version).
- Download the following files from our **[Data & Weights Registry (Google Sheet)](https://docs.google.com/spreadsheets/d/1MZ575oXGSmGBPvi-qBCIgjDxeWn2tQmCxVbO9LWlhNs/edit?usp=sharing)**. 
- Correct the paths in `./run.sh` to where you have downloaded the files.
- In `./run.sh`, modify `conda activate` so your actual conda environment is activated.
- Run `./run.sh`.


The metric values will be created and saved in csv files in `./NonGit/OutputRuns/`.
