# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:02:25 2013

@author: team
"""

import hydro_wrapper

hydro_wrapper.change_param(12,9, 180, 1, 1, 1)
hydro_wrapper.change_data_time(12,12,14,2013,12,23,2013)
hydro_wrapper.timer(15600)
hydro_wrapper.runSELFE(12,10800)
