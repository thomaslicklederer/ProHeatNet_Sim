import pickle
import math
import pprint
import xlsxwriter
from datetime import datetime


def excelexport(setup_name):
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all
    rights reserved'''
    # create datetime string
    now = datetime.now() # current date and time
    date_time = now.strftime("%Y%m%d_%H%M%S")


    # parameters
    case_name = '%s_%s' % (setup_name, date_time)

    # load data from pickle file
    myfile = open('results/results.pkl', 'rb')
    myresults = pickle.load(myfile)
    #pprint.pprint(myresults)
    myfile.close()
    PSM = myresults['PSM']
    v_vec = myresults['v_vec']
    e_vec = myresults['e_vec']
    topology = myresults['topology']
    scenario = myresults['scenario']
    solutions = myresults['solutions']

    # create excel file and worksheet
    workbook = xlsxwriter.Workbook('results/%s.xlsx' % case_name)
    worksheet = workbook.add_worksheet(case_name)
    worksheet.write(1, 0, case_name)

    # content
    row = 0
    col = 1

    # transferred heat
    for P in PSM:
        header = 'dotQ_%1.0f' % P
        entry = -solutions['Q_trnsf2'][P]
        
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        
        col +=1
        
    # sum of heat losses
    sum_Q_loss = 0
    for key in solutions['Q_loss'].keys():
        sum_Q_loss += solutions['Q_loss'][key]
    header = 'sum_dotQ_loss'
    entry = sum_Q_loss
    worksheet.write(0, col, header)
    worksheet.write(1, col, entry)
    col +=1

    # network (primary side) temperatures at prosumers
    for P in PSM:
        header = 'T_%1.0f_prim_hot' % P
        key= ('%1.0fh' % P, '%1.0fc' % P)
        entry = solutions['T'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1
        
        header = 'T_%1.0f_prim_cold' % P
        key= ('%1.0fc' % P, '%1.0fh' % P)
        entry = solutions['T'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1

    # prosumer (secondary side) temperatures at prosumers
    for P in PSM:
        header = 'T_%1.0f_sec_hot' % P
        key= 'PSM%1.0fh' % P
        entry = solutions['T'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1
        
        header = 'T_%1.0f_sec_cold' % P
        key= 'PSM%1.0fc' % P
        entry = solutions['T'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1

    # volume flows through prosumer substations at primary side
    for P in PSM:
        header = 'dotV_%1.0f_prim' % P
        key= ('%1.0fh' % P, '%1.0fc' % P)
        entry = solutions['dotV'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1
     
    # volume flows through prosumer substations at secondary side
    for P in PSM:
        header = 'dotV_%1.0f_sec' % P
        entry = solutions['dotV'][P]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1

    # pressure differences over substations primary side
    for P in PSM:
        header = 'Deltap_%1.0f_prim' % P
        key= ('%1.0fh' % P, '%1.0fc' % P)
        entry = solutions['Deltap'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1

    #############
    col += 1
    #############

    # volume flows in network
    for key in solutions['dotV'].keys():
        try:
            header = 'dotV_%s_%s' % (key[0], key[1])
        except:
            header = 'dotV_%s' % key
        entry = solutions['dotV'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1
        
    # pressure differences in network
    for key in solutions['Deltap'].keys():
        if len(key)>1:
            header = 'Deltap_%s_%s' % (key[0], key[1])
        else:
            header = 'Deltap_%s' % key
        entry = solutions['Deltap'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1
        
    # temperatures in network
    for key in solutions['T'].keys():
        if not 'PSM' in key:
            header = 'T_%s_%s' % (key[0], key[1])
        else:
            header = 'T_%s' % key
        entry = solutions['T'][key]
        worksheet.write(0, col, header)
        worksheet.write(1, col, entry)
        col += 1

    # close workbook
    workbook.close()