# OL-to-SR
*Source code for results presented in the paper "From order lifting to social ranking: recovering preferences from partial extensions" for ECAI25*

Launching by writing command
```
python3 main.py
```
once in the OL-to-SR repository, will start the tests. 

Parameters are set by default to 10 000 runs being made to test the exactness of the recovery, the exactness of the top item retrieved, the Kendall-Tau distance to the ground truth and the number of errors found. They will be made based on rankings over all coalitions in a population N ranging from 4 to 9. It is possible to change these default parameters in the Terminal. **Please keep in mind that tests may take a while to run, especially as the population size increases.**

Results from tests will be saved in two formats in a folder "out":

```bash
OL-to-SR
├── out
│   ├── data
│   ├── plots
├── launch.py
├── main.py
├── OL.py
├── SR.py
├── tools.py
├── README.md
```
Subfolder "data" stores raw data from run results in csv files; "plots" stores plotted results.

For legibility purposes, results are stored in files titled after the content of the results ("exact" for study of number of times each method recover the correct exact order; "win" for when methods correctly recover the top item; "KT" for Kendall-Tau distance studies; "error" for the study of number of errors). Note that, due to the complexity of the results for the studies of scenarios where rankings are over coalitions of size k, only raw data in csv files will be saved.
