# ProHeatNet_Sim
A python based simulation framework for the thermohydraulic simulation of prosumer-dominated heat networks.

# How to run the example
- Download / clone the latest release.
- Make sure our system meets the requirements (see below):
  - python 3.7 or higher installed
  - Gurobi optimizer version 9.0.1 or higher installed
  - suitable license for Gurobi
  - module "gurobipy" installed
  - IPython (Interactive Python) version 7.12.0 or higher installed
- open command shell of your operating system
- start ipython
  ```
  ipython
  ```
- navigate to the directory that you downloaded the repository to
  ```
  cd *your_directory*
  ```
- do some settings to automatically reload changes in the python code files
  ```
  %load_ext autoreload
  %autoreload 2
  ```
- run the main file
  ```
  run A_main.py
  ```
- the example set-up in directory *./set-up/example* will be used
- calculation starts
- results show in command window
- figures with visualizations of the results should show
- results are printed to textfile *results.txt* and saved as a python pickle in *results.pkl*
- export to Excel-file can be done by running *H_toExcel.py*
```
run H_toExcel.py
```

# How to modify the set-up
- go to directory *./set-up/example* and copy the three files
- paste the copied files in folder *./set-up/my_set-up*
- 

# Requirements
The code was tested under Linux and under Windows with the following specifications.
## Linux
- Operating system: Ubuntu 18.04.5 LTS
- Python version: Python 3.7.6
- Optimizer: Gurobi Optimizer version 9.0.1 build v9.0.1rc0 (linux64), Academic license - for non-commercial use only
- install module "gurobipy" (https://www.gurobi.com/documentation/9.0/quickstart_mac/the_grb_python_interface_f.html)
- IDE: none, instead: IPython (Interactive Python) version 7.12.0
## Windows
- operating system: Microsoft Windows Version 10.0.18363.1139
- Python version: Python 3.8.3
- Optimizer: Gurobi Optimizer version 9.1.0 build v9.1.0rc0 (win64)
- install module "gurobiypy" (https://www.gurobi.com/documentation/9.0/quickstart_mac/the_grb_python_interface_f.html)
- - IDE: none, instead: IPython (Interactive Python) version 7.16.1
