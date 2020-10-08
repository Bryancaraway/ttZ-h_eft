'''
Store unc normalizations for 
Signal and background
and tt_bb

dict format
year:process:uncertainty:normalization
'''

shape_unc_norms = { 
    '2016':
    {
        'ttH':
        {
            'mu_r_Up'       :1.154172219,
            'mu_r_Down'     :0.7947879344,
            'mu_f_Up'       :1.019241013,
            'mu_f_Down'     :0.9749255808,
            'pdfWeight_Up'  :0.9785013492,
            'pdfWeight_Down':1.022464567
        },
        'ttZ':
        {
            'mu_r_Up'       :1.08959876,
            'mu_r_Down'     :0.9270304533,
            'mu_f_Up'       :1.022914196,
            'mu_f_Down'     :0.9733106715,
            'pdfWeight_Up'  :0.9639621111,
            'pdfWeight_Down':1.038837103,
        },
        'ttHbb':
        {
            'mu_r_Up'       :1.223879752,
            'mu_r_Down'     :0.708285456,
            'mu_f_Up'       :1.019004895,
            'mu_f_Down'     :0.9754186276,
            'pdfWeight_Up'  :0.9886965412,
            'pdfWeight_Down':1.011564906
        },
        'ttZbb':
        {
            'mu_r_Up'       :1.089579691,
            'mu_r_Down'     :1.078785105,
            'mu_f_Up'       :0.9773626304,
            'mu_f_Down'     :1.027684208,
            'pdfWeight_Up'  :1.037516441,
            'pdfWeight_Down':0.9624835586
        },
        'Vjets':
        {
            'mu_r_Up'       :1.041775985,
            'mu_r_Down'     :0.9517037872,
            'mu_f_Up'       :1.004489954,
            'mu_f_Down'     :1.002978595,
            'pdfWeight_Up'  :0.9451599757,
            'pdfWeight_Down':1.061595861
        },
        'ttX':
        {
            'mu_r_Up'       :1.054942922,
            'mu_r_Down'     :0.9481289124,
            'mu_f_Up'       :0.9980799314,
            'mu_f_Down'     :0.998019411,
            'pdfWeight_Up'  :0.9730232084,
            'pdfWeight_Down':1.028515294
        },
        'other':
        {
            'mu_r_Up'       :1.024665248,
            'mu_r_Down'     :0.972323301,
            'mu_f_Up'       :0.9871896069,
            'mu_f_Down'     :1.016601997,
            'pdfWeight_Up'  :0.9642494407,
            'pdfWeight_Down':1.038503611
        },
        'TTBar':
        {
            'mu_r_Up'       :1.11037998,
            'mu_r_Down'     :0.9034428481,
            'mu_f_Up'       :1.019278563,
            'mu_f_Down'     :0.9759970997,
            'pdfWeight_Up'  :0.9916673881,
            'pdfWeight_Down':1.00847383
        },
        'tt_bb':
        {
            'mu_r_Up'       :1.36878335,
            'mu_r_Down'     :0.7585017229,
            'mu_f_Up'       :1.020980265,
            'mu_f_Down'     :0.9708575015,
            'pdfWeight_Up'  :0.9674025093,
            'pdfWeight_Down':1.034746148,
            'ISR_Up'        :1.001090801,
            'ISR_Down'      :0.9975892668,
            'FSR_Up'        :0.8628071886,
            'FSR_Down'      :0.9104763627
            
        },
        'tt_2b':
        {
            'mu_r_Up'       :1.36878335,
            'mu_r_Down'     :0.7585017229,
            'mu_f_Up'       :1.020980265,
            'mu_f_Down'     :0.9708575015,
            'pdfWeight_Up'  :0.9674025093,
            'pdfWeight_Down':1.034746148,
            'ISR_Up'        :1.001090801,
            'ISR_Down'      :0.9975892668,
            'FSR_Up'        :0.8628071886,
            'FSR_Down'      :0.9104763627
        },
    },
    '2017':
    {
        'ttH':
        {
            'mu_r_Up'       :1.153983255,
            'mu_r_Down'     :0.7949417097,
            'mu_f_Up'       :1.019145214,
            'mu_f_Down'     :0.975039758,
            'pdfWeight_Up'  :0.9779578227,
            'pdfWeight_Down':1.023058706
        },
        'ttZ':
        {
            'mu_r_Up'       :1.089343077,
            'mu_r_Down'     :0.9272916579,
            'mu_f_Up'       :0.9272916579,
            'mu_f_Down'     :0.9739511566,
            'pdfWeight_Up'  :0.9971435622,
            'pdfWeight_Down':1.00287285
        },
        'new_ttZbb': # not verified
        {
            'mu_r_Up'       :1.122440087,
            'mu_r_Down'     :0.9270244634,
            'mu_f_Up'       :1.052569601,
            'mu_f_Down'     :0.9741160454,
            'pdfWeight_Up'  :0.9895149904,
            'pdfWeight_Down':1.01070959
        },
        'ttHbb':
        {
            'mu_r_Up'       :1.223901083,
            'mu_r_Down'     :0.7083289977,
            'mu_f_Up'       :1.018935143,
            'mu_f_Down'     :0.9755043454,
            'pdfWeight_Up'  :0.9876652389,
            'pdfWeight_Down':1.01264675
        },
        'ttZbb':
        {
            'mu_r_Up'       :1.089403551,
            'mu_r_Down'     :0.9271305916,
            'mu_f_Up'       :1.022520642,
            'mu_f_Down'     :0.9737183687,
            'pdfWeight_Up'  :0.998941416,
            'pdfWeight_Down':1.00106083
        },
        'Vjets':
        {
            'mu_r_Up'       :1.034439556,
            'mu_r_Down'     :0.9608115094,
            'mu_f_Up'       :0.9867753929,
            'mu_f_Down'     :1.029381929,
            'pdfWeight_Up'  :0.9880826901,
            'pdfWeight_Down':1.01220829
        },
        'ttX':
        {
            'mu_r_Up'       :1.067569182,
            'mu_r_Down'     :0.936076902,
            'mu_f_Up'       :0.9966101886,
            'mu_f_Down'     :0.998655213,
            'pdfWeight_Up'  :0.9985193727,
            'pdfWeight_Down':1.001485025
        },
        'other':
        {
            'mu_r_Up'       :1.023090294,
            'mu_r_Down'     :0.9740709749,
            'mu_f_Up'       :0.9892341041,
            'mu_f_Down'     :1.013858917,
            'pdfWeight_Up'  :0.9904350087,
            'pdfWeight_Down':1.009751538
        },
        'TTBar':
        {
            'mu_r_Up'       :1.110371976,
            'mu_r_Down'     :0.9034522385,
            'mu_f_Up'       :1.01928524,
            'mu_f_Down'     :0.9759878037,
            'pdfWeight_Up'  :0.9912449624,
            'pdfWeight_Down':1.008911071
        },
        'tt_bb':
        {
            'mu_r_Up'       :1.368920511,
            'mu_r_Down'     :0.7582487489,
            'mu_f_Up'       :1.020871941,
            'mu_f_Down'     :0.9709472992,
            'pdfWeight_Up'  :0.9675287834,
            'pdfWeight_Down':1.034726436,
            'ISR_Up'        :1.001079023,
            'ISR_Down'      :0.9977401484,
            'FSR_Up'        :0.9592251866,
            'FSR_Down'      :0.938298306
        },
        'tt_2b':
        {
            'mu_r_Up'       :1.368920511,
            'mu_r_Down'     :0.7582487489,
            'mu_f_Up'       :1.020871941,
            'mu_f_Down'     :0.9709472992,
            'pdfWeight_Up'  :0.9675287834,
            'pdfWeight_Down':1.034726436,
            'ISR_Up'        :1.001079023,
            'ISR_Down'      :0.9977401484,
            'FSR_Up'        :0.9592251866,
            'FSR_Down'      :0.938298306
        },
    },
    '2018':
    {
        'ttH':
        {
            'mu_r_Up'       :1.223951931,
            'mu_r_Down'     :0.7083202143,
            'mu_f_Up'       :1.019029157,
            'mu_f_Down'     :0.9753779636,
            'pdfWeight_Up'  :0.9880984117,
            'pdfWeight_Down':1.012191792
        },
        'ttZ':
        {
            'mu_r_Up'       :1.089358994,
            'mu_r_Down'     :0.9272679036,
            'mu_f_Up'       :1.022246366,
            'mu_f_Down'     :0.9740186372,
            'pdfWeight_Up'  :0.9900751773,
            'pdfWeight_Down':1.010125817
        },
        'ttHbb':
        {
            'mu_r_Up'       :1.22396623,
            'mu_r_Down'     :0.7083418477,
            'mu_f_Up'       :1.019060665,
            'mu_f_Down'     :0.9753604582,
            'pdfWeight_Up'  :0.9880306323,
            'pdfWeight_Down':1.012262927
        },
        'ttZbb':
        {
            'mu_r_Up'       :1.08941883,
            'mu_r_Down'     :0.9271044569,
            'mu_f_Up'       :1.022422238,
            'mu_f_Down'     :0.9738510912,
            'pdfWeight_Up'  :0.9895278558,
            'pdfWeight_Down':1.010696168
        },
        'new_ttZbb': # not verified
        {
            'mu_r_Up'       :1.08941883,
            'mu_r_Down'     :0.9271044569,
            'mu_f_Up'       :1.022422238,
            'mu_f_Down'     :0.9738510912,
            'pdfWeight_Up'  :0.9895278558,
            'pdfWeight_Down':1.010696168
        },
        'Vjets':
        {
            'mu_r_Up'       :1.034406517,
            'mu_r_Down'     :0.9608491273,
            'mu_f_Up'       :0.98674442,
            'mu_f_Down'     :1.029418647,
            'pdfWeight_Up'  :0.9880942359,
            'pdfWeight_Down':1.012196174
        },
        'ttX':
        {
            'mu_r_Up'       :1.057770602,
            'mu_r_Down'     :0.9517853584,
            'mu_f_Up'       :0.9935263339,
            'mu_f_Down'     :0.9935263339,
            'pdfWeight_Up'  :0.9935285463,
            'pdfWeight_Down':1.006556311
        },
        'other':
        {
            'mu_r_Up'       :1.021577146,
            'mu_r_Down'     :0.9754442119,
            'mu_f_Up'       :0.9883471812,
            'mu_f_Down'     :1.015073056,
            'pdfWeight_Up'  :0.9901528784,
            'pdfWeight_Down':1.010044949
        },
        'TTBar':
        {
            'mu_r_Up'       :1.110374282,
            'mu_r_Down'     :0.9034496276,
            'mu_f_Up'       :1.01928302,
            'mu_f_Down'     :0.975990606,
            'pdfWeight_Up'  :0.990834563,
            'pdfWeight_Down':1.009336585
        },
        'tt_bb':
        {
            'mu_r_Up'       :1.368222738,
            'mu_r_Down'     :0.7581770089,
            'mu_f_Up'       :1.020928553,
            'mu_f_Down'     :0.9708995389,
            'pdfWeight_Up'  :0.9673913873,
            'pdfWeight_Down':1.034883626,
            'ISR_Up'        :1.001132127,
            'ISR_Down'      :0.9976717585,
            'FSR_Up'        :0.934261331,
            'FSR_Down'      :0.9521073108
        },
        'tt_2b':
        {
            'mu_r_Up'       :1.368222738,
            'mu_r_Down'     :0.7581770089,
            'mu_f_Up'       :1.020928553,
            'mu_f_Down'     :0.9708995389,
            'pdfWeight_Up'  :0.9673913873,
            'pdfWeight_Down':1.034883626,
            'ISR_Up'        :1.001132127,
            'ISR_Down'      :0.9976717585,
            'FSR_Up'        :0.934261331,
            'FSR_Down'      :0.9521073108
        },
    }
}

rate_unc = {
    'QCD_scale':
    {
        'Vjets':[1.012,.980],
        'ttX'  :[1.071,.938],
        'other':[1.012,.992]
    },
    'pdf':
    {
        'Vjets':1.012,
        'ttX'  :1.007,
        'other':1.01
    }
}
    
