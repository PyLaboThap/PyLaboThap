from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc

case_study = "RecupTest"

#Exo ORC M&S
HTX = HexCstEffDisc()

"Simple test - CO2 HTHP GasCooler"

if case_study == 'GasCooler':

    HTX.set_inputs(
        fluid_C = 'Water',
        T_su_C = 273.15 + 15,
        m_dot_C = 0.1,
        P_su_C = 10e5,
    
        fluid_H = 'CO2',
        T_su_H = 450,
        m_dot_H = 0.16,
        P_su_H = 180*1e5,
    )
    
    HTX.set_parameters(**{
        'eta_max' : 0.95,
        'n_disc' : 100, 
        'Pinch_min' : 10,
        'DP_c' : 50*1e3,
        'DP_h' : 50*1e3,    
    })
    
    HTX.solve()
    HTX.plot_disc()
    
    # fig = HTX.plot_Ts(choose_HX_side='H')

elif case_study == "RecupHT":
    
    HTX.set_inputs(
        fluid_C = 'CO2',
        T_su_C = 307.7544676675842,
        m_dot_C = 30,
        P_su_C = 3991606,
    
        fluid_H = 'CO2',
        T_su_H = 325.403943370302,
        m_dot_H = 30,
        P_su_H = 150*1e5,
    )
    
    HTX.set_parameters(**{
        'eta_max' : 0.95,
        'n_disc' : 20, 
        'Pinch_min' : 0,
        'DP_c' : 0,
        'DP_h' : 0*1e3,    
    })
    
    HTX.solve()
    HTX.plot_disc()
    

elif case_study == "RecupTest":
    
    HTX.set_inputs(
        fluid_C = 'CO2',
        T_su_C = 306.51,
        m_dot_C = 365.4,
        P_su_C = 17830200,
    
        fluid_H = 'CO2',
        T_su_H = 324.403943370302,
        m_dot_H = 365.4,
        P_su_H = 5965838,
    )
    
    HTX.set_parameters(**{
        'eta_max' : 0.95,
        'n_disc' : 100, 
        'Pinch_min' : 0,
        'DP_c' : 0,
        'DP_h' : 0*1e3,    
    })
    
    HTX.solve()
    HTX.plot_disc()
    
