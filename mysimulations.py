import glob
import os
from datetime import datetime
import openpyxl as pyxl
import compare2excel


now = datetime.now() # current date and time
date_time = now.strftime("%Y%m%d_%H%M%S")
standarddir = os.getcwd()

# set base case 
basecase = "A"
path_config     =   "simulationcases/%s/setup/configuration.py" % basecase
path_param      =   "simulationcases/%s/setup/parametrization.py" % basecase
path_res = "simulationcases/%s/results/" % (basecase)
workbook_path = path_res + 'results_%s_%s.xlsx' % (basecase,date_time)


# read all scenario variations
os.chdir("simulationcases/%s/" % basecase)
scenarios = list(glob.glob("*.py"))
os.chdir(standarddir)
scenarios = sorted(scenarios)
scenarios = [entry.replace(".py", "") for entry in scenarios]


# initiate excel file
workbook = pyxl.Workbook(workbook_path)
workbook.save(workbook_path)
workbook = pyxl.load_workbook(workbook_path)
sheet = workbook.active

# calculate and save all scenario results
for i in range(len(scenarios)):
    scenario = scenarios[i]
    iteration = i+2
    setup_name = scenario
    path_scenario   =   "simulationcases/%s/%s.py" % (basecase, scenario)
    print(path_scenario)
    exec(open("A_main.py").read())
    compare2excel.excelcompare(setup_name , workbook_path, iteration)