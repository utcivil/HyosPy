# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:07:23 2013

@author: team
"""

import oilspill_wrapper

oilspill_wrapper.mul_GNOME_inputs(12,2013,12,21,0,48)

oilspill_wrapper.run_mul_GNOME(12,666393,3076780,2013,12,21,48,900)

oilspill_wrapper.GNOME_GM_visualization(12)

oilspill_wrapper.GNOME_GE_animation(12,13,2013,12,21)

# oilspill_wrapper.rk4(12,900,900,48,666393,3076780)
